function Child = Mutation_SOR_Overlapping(Child, A, ...
                                           IndicesInteractionProtein, NumInteractionProtein, ...
                                           Pm, params)
if nargin < 6, params = struct(); end
if ~isfield(params,'MinSize'), params.MinSize = 3; end
if ~isfield(params,'pAddTriangle'), params.pAddTriangle = 0.6; end
if ~isfield(params,'pPrune'), params.pPrune = 0.4; end
if ~isfield(params,'MaxComplexAddsPerNode'), params.MaxComplexAddsPerNode = 2; end
if ~isfield(params,'MaxMovesPerNode'), params.MaxMovesPerNode = 1; end

A = sparse(A);
N = size(A,1);

Comms = Child.CommsNodes.CommsNodes;
K = numel(Comms);

% ---------- Build membership map: complexes containing each node ----------
proteinComplexMap = cell(N,1);
for c = 1:K
    members = unique(Comms{c}(:))';
    members = members(members>=1 & members<=N);
    Comms{c} = members;
    for ii = 1:numel(members)
        u = members(ii);
        proteinComplexMap{u}(end+1) = c;
    end
end

% Helper: fast membership check using logical mask per complex (rebuild lazily)
% For simplicity, we build masks once (cost O(sum |C|))
Cmask = false(N, K);
for c = 1:K
    Cmask(Comms{c}, c) = true;
end

% ---------- Mutation over nodes ----------
for u = 1:N
    if NumInteractionProtein(u) <= 0 || rand > Pm
        continue;
    end

    currCs = proteinComplexMap{u};
    if isempty(currCs)
        continue;
    end

    % ========= (A) Boundary-Repair Move =========
    % Evaluate if u is "boundary-heavy" in any of its current complexes
    moveCount = 0;
    for idxC = 1:numel(currCs)
        if moveCount >= params.MaxMovesPerNode, break; end

        c = currCs(idxC);
        membersMask = Cmask(:,c);

        % Compute in/out edge counts of u relative to complex c
        kin = 0; kout = 0;
        for t = 1:NumInteractionProtein(u)
            v = IndicesInteractionProtein(u,t);
            if v < 1 || v > N, continue; end
            if membersMask(v)
                kin  = kin  + A(u,v);
            else
                kout = kout + A(u,v);
            end
        end

        % If u leaks more than it connects inside -> try relocate/add elsewhere
        if kin <= kout
            bestC = [];
            bestScore = -Inf;

            % Try candidate complexes from neighbors' complexes (more focused than all K)
            cand = [];
            for t = 1:NumInteractionProtein(u)
                v = IndicesInteractionProtein(u,t);
                if v>=1 && v<=N
                    cand = [cand, proteinComplexMap{v}]; %#ok<AGROW>
                end
            end
            cand = unique(cand);
            cand(cand==c) = [];

            for cc = cand
                % Skip if already in that complex
                if Cmask(u,cc), continue; end

                mask2 = Cmask(:,cc);

                kin2 = 0; kout2 = 0;
                for t = 1:NumInteractionProtein(u)
                    v = IndicesInteractionProtein(u,t);
                    if v < 1 || v > N, continue; end
                    if mask2(v)
                        kin2  = kin2  + A(u,v);
                    else
                        kout2 = kout2 + A(u,v);
                    end
                end

                % Triangle-closure gain approximation:
                % triangles gained ~ number of common neighbors of u with nodes in complex
                % Approx via: count neighbors v in complex where u connects AND v has another neighbor in complex
                triGain = 0;
                for t = 1:NumInteractionProtein(u)
                    v = IndicesInteractionProtein(u,t);
                    if v<1 || v>N || ~mask2(v) || A(u,v)==0, continue; end
                    % count v's neighbors inside complex
                    vn = IndicesInteractionProtein(v, 1:NumInteractionProtein(v));
                    vn = vn(vn>=1 & vn<=N);
                    triGain = triGain + sum(mask2(vn));
                end

                % Score encourages high kin2, low kout2, high triGain
                score = (kin2 - kout2) + 0.1*triGain;

                if score > bestScore
                    bestScore = score;
                    bestC = cc;
                end
            end

            if ~isempty(bestC)
                % Add u to best complex (overlap allowed)
                Comms{bestC}(end+1) = u;
                Cmask(u,bestC) = true;
                proteinComplexMap{u}(end+1) = bestC;
                moveCount = moveCount + 1;
            end
        end
    end

    % ========= (B) Triangle-Closure Add =========
    if rand <= params.pAddTriangle
        adds = 0;

        % Consider complexes of neighbors as targets
        cand = [];
        for t = 1:NumInteractionProtein(u)
            v = IndicesInteractionProtein(u,t);
            if v>=1 && v<=N
                cand = [cand, proteinComplexMap{v}]; %#ok<AGROW>
            end
        end
        cand = unique(cand);

        for cc = cand
            if adds >= params.MaxComplexAddsPerNode, break; end
            if Cmask(u,cc), continue; end  % already member

            mask = Cmask(:,cc);

            % Triangle closure condition: u has at least 2 neighbors inside cc
            nin = 0;
            for t = 1:NumInteractionProtein(u)
                v = IndicesInteractionProtein(u,t);
                if v>=1 && v<=N && mask(v) && A(u,v)==1
                    nin = nin + 1;
                    if nin >= 2, break; end
                end
            end

            if nin >= 2
                Comms{cc}(end+1) = u;
                Cmask(u,cc) = true;
                proteinComplexMap{u}(end+1) = cc;
                adds = adds + 1;
            end
        end
    end

    % ========= (C) Prune Triangle-Poor Boundary Nodes =========
    if rand <= params.pPrune
        % For each complex containing u, remove u if it is boundary-heavy AND triangle-poor
        currCs2 = unique(proteinComplexMap{u});
        for idxC = 1:numel(currCs2)
            c = currCs2(idxC);
            members = unique(Comms{c}(:))';
            if numel(members) <= params.MinSize, continue; end

            mask = Cmask(:,c);

            kin = 0; kout = 0; triSupport = 0;
            for t = 1:NumInteractionProtein(u)
                v = IndicesInteractionProtein(u,t);
                if v<1 || v>N, continue; end

                if mask(v)
                    kin = kin + A(u,v);

                    % tri-support approx: v's neighbors inside complex
                    vn = IndicesInteractionProtein(v, 1:NumInteractionProtein(v));
                    vn = vn(vn>=1 & vn<=N);
                    triSupport = triSupport + sum(mask(vn));
                else
                    kout = kout + A(u,v);
                end
            end

            % Remove if leaks more than connects AND has weak triangle support
            if (kin <= kout) && (triSupport <= 1)
                % remove u from this complex
                Comms{c} = members(members ~= u);
                Cmask(u,c) = false;
                % also remove from map
                proteinComplexMap{u} = proteinComplexMap{u}(proteinComplexMap{u} ~= c);
            end
        end
    end

end

% Final clean: unique members in each complex
for c = 1:K
    Comms{c} = unique(Comms{c}(:))';
end

Child.CommsNodes.CommsNodes = Comms;
end
