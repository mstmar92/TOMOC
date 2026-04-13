function [Individual] = ComputeCmplxDecodingOverlapping(Individual, ProteinLabelPairs)

% 1) Extraction process of the communities of edges

chromosome = Individual.Chromosome; % Alleles of the individual
numEdges = length(chromosome); % Total number of edges

% Step 1: Build adjacency list (undirected)
adjList = cell(1, numEdges);
for i = 1:numEdges
    j = chromosome(i);
    if (j > 0)
    % Add bi-directional connection
    adjList{i} = [adjList{i}, j];
    adjList{j} = [adjList{j}, i];
    end
end

% Step 2: Initialize
visited = false(1, numEdges);
commEdges = {}; % List of communities
stack = [];

% Step 3: DFS traversal for each unvisited edge
for i = 1:numEdges
    if ~visited(i)
        comm = []; % New community
        stack = [i]; % Start DFS from this edge
        visited(i) = true;
        
        while ~isempty(stack)
            current = stack(end);
            stack(end) = []; % Pop
            comm = [comm, current];
            
            for neighbor = adjList{current}
                if ~visited(neighbor)
                    visited(neighbor) = true;
                    stack(end+1) = neighbor; % Push
                end
            end
        end
        
        commEdges{end + 1} = sort(comm); % Store sorted community
    end
end

% % Store the communities in the output structure
% Individual.Communities = commEdges;


% function commsNodes = TransformEdgeToNodeCommunities(commsEdges, ProteinLabelPairs)
% 2) Transform communities of edges into overlapping communities of nodes
% commsEdges: Cell array where each element is a set of edges in a community
% ProteinLabelPairs: Nx2 matrix where each row represents an edge with source and target nodes

    commsNodes = {}; % Initialize the communities of nodes

    for i = 1:length(commEdges)
        commE = commEdges{i}; % Extract the current community of edges
        commN = []; % Initialize the current community of nodes

        for j = 1:length(commE)
            edgeId = commE(j); % Get the edge ID

            % Get the source and target nodes for the edge
            idNodeSource = ProteinLabelPairs(edgeId, 1);
            idNodeTarget = ProteinLabelPairs(edgeId, 2);

            % Add the source and target nodes to the node community
            commN = unique([commN, idNodeSource, idNodeTarget]);
        end

        % Add the current node community to the list of communities
        commsNodes{end+1} = commN;
    end
    Individual.CommsNodes = commsNodes;


