function [Population] = ComputeOurFitnessCollectionOverlapping(A, N, ...
                                                      IndicesInteractionProtein, NumInteractionProtein, ...
                                                      Population)

% ComputeOurFitnessCollectionOverlapping (NEW)
% Two-objective topology-based evaluation for MOEA/D (Tchebycheff, Minimization)
%
% Obj1 (min): Avg Conductance across detected complexes
% Obj2 (min): 1 - Avg Triangle Density across detected complexes
%
% Notes:
% - A must be NxN binary adjacency (0/1), undirected. Diagonal will be zeroed.
% - Distance / IndicesInteractionProtein / NumInteractionProtein are kept for signature compatibility,
%   but not used in this NEW objective set.
%
% Output stored in:
%   Population.Obj = [Obj1, Obj2]
%   Population.Obj1_Conductance, Population.Obj2_TriLack, Population.Kv

% ---- Hyperparameters ----
MinSize = 3;
epsVal  = 1e-12;

% ---- Prepare graph ----
A = sparse(A);
A(1:N+1:end) = 0;

deg = full(sum(A,2));      % degree of each node
totalVol = sum(deg);       % = 2M

% ---- Read complexes ----
CmplxID = Population.CommsNodes.CommsNodes;
K = numel(CmplxID);

condSum = 0;
triSum  = 0;
Kv = 0;

for c = 1:K
    members = unique(CmplxID{c}(:))';
    s = numel(members);

    if s < MinSize
        continue;
    end

    Kv = Kv + 1;

    % Subgraph
    A_sub = A(members, members);

    % Internal edges (binary undirected)
    m_in = nnz(triu(A_sub, 1));

    % vol(C)
    volC = sum(deg(members));

    % cut(C) = vol(C) - 2*m_in
    cutC = volC - 2*m_in;

    % Conductance
    denom = min(volC, totalVol - volC) + epsVal;
    phi = cutC / denom;
    condSum = condSum + phi;

    % Triangle density
    % triangles = trace(A_sub^3)/6
    triangles = trace(A_sub * A_sub * A_sub) / 6;
    triDen = triangles / (nchoosek(s,3) + epsVal);
    triSum = triSum + triDen;
end

% ---- Aggregate objectives (minimization) ----
if Kv == 0
    Obj1 = 1;     % worst
    Obj2 = 1;     % worst (1 - 0)
    avgTri = 0;
else
    Obj1 = condSum / Kv;      % minimize
    avgTri = triSum / Kv;     % maximize but we convert
    Obj2 = 1 - avgTri;        % minimize
end

Population.Obj = [Obj1, Obj2];
Population.Obj1_Conductance = Obj1;
Population.Obj2_TriLack = Obj2;
Population.Kv = Kv;
Population.avgTri = avgTri;





% function [ Population ] = ComputeOurFitnessCollectionOverlapping(N, m, Distance, ...
%                                                       IndicesInteractionProtein,NumInteractionProtein,...
%                                                       Population)
%                                                   
% 
% CmplxID = Population.CommsNodes.CommsNodes;
% k = length(CmplxID);
% 
% % Preallocate
% Node_i_in = zeros(1, N);
% Node_i_out = zeros(1, N);
% Cluster_k_in = zeros(1, k);
% Cluster_k_out = zeros(1, k);
% Cluster_k_Volume = zeros(1, k);
% Cluster_k_Cardinality = zeros(1, k);
% Cluster_k_Degree = zeros(1, k);
% Cluster_k_Strong_Neighbor = zeros(1, k);
% Cluster_k_Weak_Neighbor = zeros(1, k);
% Cluster_k_IntraCmplxScore = zeros(1, k);
% Cluster_k_InterCmplxScore = zeros(1, k);
% 
% % Create a map of each node to all its complexes
% NodeToComplexes = cell(1, N);
% for c = 1:k
%     for p = CmplxID{c}
%         NodeToComplexes{p}(end + 1) = c;
%     end
% end
% 
% % Compute values for each complex
% for c = 1:k
%     members = CmplxID{c};
%     Cluster_k_Cardinality(c) = length(members);
% 
%     for i = members
%         Cluster_k_Degree(c) = Cluster_k_Degree(c) + NumInteractionProtein(i);
%         Node_i_in(i) = 0;
%         Node_i_out(i) = 0;
% 
%         for j = 1:NumInteractionProtein(i)
%             neighbor = IndicesInteractionProtein(i, j);
%             if ismember(neighbor, members)
%                 % Internal link
%                 Cluster_k_Volume(c) = Cluster_k_Volume(c) + 1;
%                 Cluster_k_in(c) = Cluster_k_in(c) + 1;
%                 Node_i_in(i) = Node_i_in(i) + 1;
%             else
%                 % External link
%                 Cluster_k_out(c) = Cluster_k_out(c) + 1;
%                 Node_i_out(i) = Node_i_out(i) + 1;
%             end
%         end
% 
%         % Intra/Inter evaluation
%         if NumInteractionProtein(i) > 0
%             if Node_i_in(i) >= Node_i_out(i)
%                 Cluster_k_Strong_Neighbor(c) = Cluster_k_Strong_Neighbor(c) + (Node_i_in(i) / NumInteractionProtein(i));
%             else
%                 Cluster_k_Weak_Neighbor(c) = Cluster_k_Weak_Neighbor(c) + (Node_i_in(i) / NumInteractionProtein(i));
%             end
%         end
%     end
% end
% 
% % Intra and Inter-complex scores
% Cluster_k_IntraCmplxScore = ((Cluster_k_Volume / 2) + Cluster_k_Strong_Neighbor) ./ Cluster_k_Cardinality;
% IntraCmplxScore = sum(Cluster_k_IntraCmplxScore);
% 
% Cluster_k_InterCmplxScore = (Cluster_k_out ./ Cluster_k_Degree) + Cluster_k_Weak_Neighbor;
% InterCmplxScore = sum(Cluster_k_InterCmplxScore) * k;
% 
% Population.Intra_CmplxScore = (N ^ 2) - IntraCmplxScore;
% if k == 1
%     InterCmplxScore = (N ^ 2);
% end
% Population.Inter_CmplxScore = InterCmplxScore;
% 
% 
% TotalIntra_CC = 0;
% DistantProteinCCInEachCluster = zeros(1, k);  % For Inter-CC needs
% DistantProteinInEachCluster = zeros(1, k);    % The most distant protein in each cluster
% 
% for CurrentCluster = 1 : k
%     ProteinsInCurrentCluster = CmplxID{CurrentCluster};
%     Intra_CC = 0;
%     
%     % Handle edge case: cluster is empty or singleton
%     if length(ProteinsInCurrentCluster) <= 1
%         DistantProteinCCInEachCluster(CurrentCluster) = 0;
%         DistantProteinInEachCluster(CurrentCluster) = ProteinsInCurrentCluster(1);
%         continue;
%     end
% 
%     MaxCC = -inf; % for finding the most distant node
%     
%     for p = 1 : length(ProteinsInCurrentCluster)
%         CurrentProtein = ProteinsInCurrentCluster(p);
%         CurrentProteinCC = 0;
% 
%         for q = 1 : length(ProteinsInCurrentCluster)
%             if p ~= q
%                 NextProtein = ProteinsInCurrentCluster(q);
%                 CurrentProteinCC = CurrentProteinCC + Distance(CurrentProtein, NextProtein);
%             end
%         end
% 
%         Intra_CC = Intra_CC + CurrentProteinCC;
% 
%         if CurrentProteinCC > MaxCC
%             MaxCC = CurrentProteinCC;
%             DistantProteinCCInEachCluster(CurrentCluster) = MaxCC;
%             DistantProteinInEachCluster(CurrentCluster) = CurrentProtein;
%         end
%     end
% 
%     % Average closeness centrality (based on distance)
%     TotalIntra_CC = TotalIntra_CC + (Intra_CC / (length(ProteinsInCurrentCluster)^2));
% end
% 
% Population.Intra_CmplxScore = Population.Intra_CmplxScore + (TotalIntra_CC / k);
