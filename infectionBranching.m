%% infectionBranching.m
% This function simulates the spread of an infectious disease where the
% underlying network of interactions is developed using a branching model
% (as described in the manuscript "Simple control for complex pandemics").

%% Infection simulator for a branching model
function result = infectionBranching(diseaseParams, testingParams, modelParams, networkParams)
    
    % Relevant paramters as defined in
    % main_simpleControlForComplexPandemics.m
    Tend = modelParams.Tend;
    numIters = modelParams.numIters;
    numInf0 = modelParams.numInf0;
    
    rho = diseaseParams.rho;
    d = diseaseParams.d;
    
    c = testingParams.c;
    nu = testingParams.nu;
    
    % Define network interaction distribution
    numNodes = networkParams.numNodes;
    if networkParams.type == 1
        mu = networkParams.mu;
        
        nBinom = numNodes-1;
        pBinom = mu/nBinom;
        sigma2 = nBinom*pBinom*(1-pBinom);
        
        generateEdges = @() binornd(nBinom, pBinom);    
    else
        a = networkParams.a;
        b = networkParams.b;
        generateEdges = @() randi(b+1-a) + (a-1); 
        mu = (b+a)/2;
        sigma2 = ((b-a+1)^2 - 1)/12;
    end
    
    
    Imax = numNodes;
    
    % Define epidemiological parameters
    R0 = rho*mu*d;
    Reff = R0*(1-c*nu)/(1+nu*(d-1));
    r = nu + 1/d - nu/d;
    lambda = 1 + r*(Reff-1);
    
    infEndRes = zeros(numIters,1);
    infTotalRes = zeros(numIters,1);

    % Run stochastic infection simulation
    for i = 1:numIters
        I = zeros(Tend, 1);
        InewAll = zeros(Tend,1);
        
        I(1) = numInf0;
        InewAll(1) = numInf0;
        
        % Simulate infections through desired time window
        for t = 2:Tend
            Iold = I(t-1); 
            Inew = 0;
            if Iold > 0 && Iold < Imax
                for j = 1:Iold
                    new = 0;
                    Xjinf = generateEdges();
                    % Infection
                    for k = 1:Xjinf
                        if rho > rand(1)
                            new = new + 1;
                        end
                    end
                    
                    % Recovery
                    Zj = ((1/d) > rand(1));
                    
                    % Detection
                    Vj = (nu > rand(1));

                    new = new - Zj - Vj + Zj*Vj;
                    if Vj == 1 
                        %Xjcont = generateEdges();
                        % Contact tracing
                        for k = 1:Xjinf
                            if c*rho > rand(1)
                                new = new - 1;
                            end
                        end
                    end
                    Inew = Inew + new;
                end
                
                % Update to next time step
                I(t) = Iold + Inew;
                InewAll(t) = Inew;
                
            elseif Iold >= Imax
                I(t) = Imax;
                InewAll(t) = 0;
            else
                I(t) = 0;
                InewAll(t) = 0;
            end
        end
        
        % Store data from trial
        infEndRes(i) = I(end);
        infTotalRes(i) = sum(InewAll);
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
    result.nu = nu;
    result.c = c;
    result.networkParams = networkParams;
    result.modelParams = modelParams;
    result.diseaseParams = diseaseParams;
end



