%% generateGraphBinom.m
% This function generates a graph where the distribution in the number of
% edges connected to each node follows a binomial distribution.

function connections = generateGraphBinom(numNodes, mu)
    edgeOptions = sum(1:numNodes-1);
    edgeLocations = zeros(edgeOptions,2);
    ind = 1;
    
    % Randomly choose which potential edges become actual edges such that
    % the mean number connected to each node is mu
    for i = 1:numNodes-1
        for j = i+1:numNodes
            edgeLocations(ind, :) = [i, j];
            ind = ind+1;
        end
    end
    
    edgeOptions = length(edgeLocations(:,1));
            
    edgeChoices = randperm(edgeOptions, mu*numNodes/2);
    
    % Generate a symmetric matrix that represents connections in the graph 
    connections = zeros(numNodes, numNodes);
    for i = 1:length(edgeChoices)
        ind = edgeChoices(i);
        connections(edgeLocations(ind,1), edgeLocations(ind,2)) = 1;
        connections(edgeLocations(ind,2), edgeLocations(ind,1)) = 1;
    end
    
    
end
