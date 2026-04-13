function [Population] = CreatePopulationOverlapping(proteinPairs, ...
                                                    PopulationSize)


% Number of edges
numEdges = size(proteinPairs, 1);

% Initialize the output array
adjacentEdges = zeros(numEdges, numEdges - 1); % Each edge can have at most numEdges - 1 adjacent edges

% Iterate through each edge
for i = 1:numEdges
    % Get the nodes of the current edge
    nodeA = proteinPairs(i, 1);
    nodeB = proteinPairs(i, 2);

    % Find adjacent edges
    adjacent = [];
    for j = 1:numEdges
        if i ~= j % Avoid self-comparison
            % Check if the edges share a node
            if any(ismember(proteinPairs(j, :), [nodeA, nodeB]))
                adjacent = [adjacent, j];
            end
        end
    end

    % Store the indices of adjacent edges
    adjacentEdges(i, 1:length(adjacent)) = adjacent;

end


for IndividualCounter = 1 : PopulationSize
    for EdgeCounter = 1 : numEdges
        if (sum(adjacentEdges(EdgeCounter,:) > 0))
            Population(IndividualCounter).Chromosome(EdgeCounter) = adjacentEdges(EdgeCounter,randi(length(find(adjacentEdges(EdgeCounter,:) > 0))));
        end
    end;
end;

