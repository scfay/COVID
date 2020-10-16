%% infectionIBM.m
% This function simulates the spread of an infectious disease where the
% underlying network of interactions is developed using an individual-based
% model (as described in the manuscript "Simple control for complex
% pandemics")

%% Infection simulator for an individual based model
function result = infectionIBM(diseaseParams, testingParams, modelParams, networkParams)

     % Relevant paramters as defined in
     % main_simpleControlForComplexPandemics.m
    Tend = modelParams.Tend;
    numIters = modelParams.numIters;
    numInf0 = modelParams.numInf0;
    
    
    rho = diseaseParams.rho;
    d = diseaseParams.d;
    
    c = testingParams.c;
    nu = testingParams.nu;
    
    % Define underlying network distribution
    numNodes = networkParams.numNodes;
    if networkParams.type == 1
        generateGraph = @() generateGraphBinom(numNodes, networkParams.mu);
        mu = networkParams.mu;
        p = mu/numNodes;
        sigma2 = mu*(1-p);
    else
        a = networkParams.a;
        b = networkParams.b;
        generateGraph = @() generateGraphUnif(numNodes, a, b);
        mu = (b + a)/2;
        sigma2 = ((b-a+1)^2 - 1)/12;
    end
    
    % Define epidemiological parameters
    R0 = rho*mu*d;
    Reff = R0*(1-c*nu)/(1+nu*(d-1));
    r = nu + 1/d - nu/d;
    lambda = 1 + r*(Reff-1);
    
    infEndRes = zeros(numIters,1);
    infTotalRes = zeros(numIters,1);
    
    % Run stochastic infection simulation
    for i = 1:numIters
        statusTrack = zeros(numNodes, Tend);
           
        connections = generateGraph();           
        infInds = randperm(numNodes, numInf0);
        statusTrack(infInds,1) = 1;
            
        % Simulate infections through desired time window
        for k = 2:Tend
            statusTrack(:,k) = statusTrack(:, k-1);
            for ii = 1:numNodes
                if statusTrack(ii,k-1) == 1
                    
                    % Infection
                    for jj = 1:numNodes
                        if connections(ii,jj) == 1 && statusTrack(jj,k-1) == 0 
                            if rho > rand(1)
                                statusTrack(jj, k) = 1;
                            end
                        end
                    end
                    
                    % Recovery
                    recov = (1/d > rand(1));
                    
                    % Detection
                    tested = (nu > rand(1));
                    
                    if tested || recov
                        statusTrack(ii, k) = 2;
                        % Contact tracing
                        if tested
                            for jj = 1:numNodes
                                if connections(ii,jj) == 1 && statusTrack(jj, k-1) == 1 
                                    if c > rand(1)
                                        statusTrack(jj, k) = 2;
                                    end
                                end    
                            end
                        end
                    end

                end
            end
        end
        % Store data from trial
        infEndRes(i) = sum(statusTrack(:,Tend) == 1);
        infTotalRes(i) = numNodes - sum(statusTrack(:,Tend) == 0);
    end
    
    % Store summary data and parameters from trials
    meanInfEnd = mean(infEndRes);
    varInfEnd = var(infEndRes);
    
    meanInfTotal = mean(infTotalRes);
    varInfTotal = var(infTotalRes);
    
    result.meanInfEnd = meanInfEnd;
    result.varInfEnd = varInfEnd;
    result.meanInfTotal = meanInfTotal;
    result.varInfTotal = varInfTotal;
    result.Reff = Reff;
    result.lambda = lambda;
    result.mu = mu;
    result.sigma2 = sigma2;
    result.R0 = R0;
    result.networkParams = networkParams;
    result.modelParams = modelParams;
    result.diseaseParams = diseaseParams;
end
