%% generateGraphUnif.m
% This function generates a graph where the distribution in the number of
% edges connected to each node follows a uniform distribution.

function connections = generateGraphUnif(numNodes, a, b)
    connections = zeros(numNodes, numNodes);
    
    numCon = randi(b-a+1, numNodes, 1) + a - 1;
    [~, sortInd] = sort(numCon);
    errors = [];
    
    % Generate connections for each row 
    for i = 1:numNodes
        rowNum = sortInd(i);
        
        numConsAlready = sum(connections(rowNum,:));
        
        % Only generate new connections if the row doesn't have an
        % appropriate number of edges yet
        numNewCons = numCon(rowNum) - numConsAlready;
        if numNewCons > 0
            potentialConnections = zeros(numNodes,1);
            for j = 1:numNodes
                if connections(rowNum, j) == 0 && j ~= rowNum && sum(connections(j,:)) < numCon(j)
                    potentialConnections(j) = 1;
                end
            end

            % Randomly create new edges
            conInds = find(potentialConnections == 1);
            if length(conInds) >= numNewCons
                
                newCons = randperm(length(conInds), numNewCons);
                for j = 1:numNewCons
                    connections(rowNum, conInds(newCons(j))) = 1;
                    connections(conInds(newCons(j)), rowNum) = 1;
                end
            else
                for j = 1:length(conInds)
                    connections(rowNum, conInds(j)) = 1;
                    connections(conInds(j), rowNum) = 1;
                end
                errors = [errors; rowNum, numNewCons-length(conInds)];
            end
        end
    end
    
    % Correct any rows that weren't able to meet the uniform distibution
    % requirement
    if ~isempty(errors)
        for i = 1:length(errors(:,1))
            rowNum = errors(i,1);
            numNeeded = a - sum(connections(rowNum,:));
            if numNeeded > 0
                potentialConnections = zeros(numNodes,1);
                for j = 1:numNodes
                    if connections(rowNum, j) == 0 && j ~= rowNum && sum(connections(j,:)) < b
                        potentialConnections(j) = 1;
                    end
                end
                conInds = find(potentialConnections == 1);

                newCons = randperm(length(conInds), numNeeded);
                for j = 1:numNeeded
                    connections(rowNum, conInds(newCons(j))) = 1;
                    connections(conInds(newCons(j)), rowNum) = 1;
                end
            end
        end
    end
end