%% main_simpleControlForComplexPandemics.m
% Main code file for "Simple control for complex pandemics" by Sarah C.
% Fay, Dalton J. Jones, Munther A. Dahleh, A.E. Hosoi
%
% The purpose of the code is to simulate the spread of a disease outbreak
% in a large susceptible population. Random viral testing and detection as
% well as contact tracing are implement. Two strategies for modeling
% interactions are included: 
%
% infectionBranching.m ---simulates a branching model (as described in the
% manuscript)
%
% infectionIBM.m ---simulates an individual-based model (as described in the manuscript
%% Parameter definitions
clear all
close all


% Model Parameters
modelParams.Tend = 6;           % number of days to simulate [days]
modelParams.numIters = 100;     % number of trials at each combination of parameters
modelParams.numInf0 = 10;       % size of original outbreak [people]
    

% Disease Parameters
diseaseParams.rho = 0.025;      % transmissibility []
diseaseParams.d = 14;           % mean infectious period [days]
    
% Testing Parameters
testingParams.c = 0.5;          % fraction of contacts traced [%]

% Variables to test
numNuTests = 4;                 % number of different daily testing fractions to try     
muTests = [10, 20, 30];         % set of mean number of daily contacts to try
nodeTests = [500 1000 2000];    % set of # nodes in IBM to try

% Plotting information
lw = 1.5;
varMax = 40;
shift = 0;
symbsType1 = ["^", "s", "p"];


% Data storage location
allResults = {};
count = 1;


%% Network as a binomial distribution
networkParams.type = 1;
for m = 1:length(muTests)
    networkParams.mu = muTests(m);
    
    for s = 1:length(nodeTests)
        networkParams.numNodes = nodeTests(s);
        
        for v = 1:numNuTests
            nu1 = rand(1);
            nu2 = rand(1);
            
            % Run simulation for branching model
            testingParams.nu = nu1;
            result1 = infectionBranching(diseaseParams, testingParams, modelParams, networkParams); 
            colorVal1 = (result1.sigma2+shift)/(varMax+shift)*ones(1,3);
            allResults{count} = result1;
            count = count + 1;
            
            % Run simulation for individual-based model
            testingParams.nu = nu2;
            result2 = infectionIBM(diseaseParams, testingParams, modelParams, networkParams);
            colorVal2 = (result2.sigma2+shift)/(varMax+shift)*ones(1,3);
            allResults{count} = result2;
            count = count + 1;
            
            
            % Plot mean number of infectious people vs. daily testing
            % fraction
            figure(1)
            hold on
            plot(nu1, result1.meanInfEnd, symbsType1(m), 'MarkerSize', 8, 'MarkerFaceColor', colorVal1, 'MarkerEdgeColor', colorVal1, 'LineWidth', lw)
            plot(nu2, result2.meanInfEnd, symbsType1(m), 'MarkerSize', 8, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', colorVal2, 'LineWidth', lw)
            
            % Plot mean number of infectious people vs. lambda
            figure(2)
            hold on
            plot(result1.lambda, result1.meanInfEnd, symbsType1(m), 'MarkerSize', 8, 'MarkerFaceColor', colorVal1, 'MarkerEdgeColor', colorVal1, 'LineWidth', lw)
            plot(result2.lambda, result2.meanInfEnd, symbsType1(m), 'MarkerSize', 8, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', colorVal2, 'LineWidth', lw)
             
        end
    end
end




%% Network as a uniform distribution

as = [muTests-1; muTests-4; muTests-8];
bs = [muTests+1; muTests+4; muTests+8];
networkParams.numNodes = 1000;
symbsType2 = ["v", "d", "h"];

networkParams.type = 2;
for m = 1:length(as(:,1))
    for s = 1:length(as(1,:))
        networkParams.a = as(s,m);
        networkParams.b = bs(s,m);
        
        for v = 1:numNuTests             
            nu1 = rand(1);
            nu2 = rand(1);
            
            % Run infection simulation for branching model
            testingParams.nu = nu1;
            result1 = infectionBranching(diseaseParams, testingParams, modelParams, networkParams); 
            colorVal1 = (result1.sigma2+shift)/(varMax+shift)*ones(1,3);
            allResults{count} = result1;
            count = count + 1;
            
            % Run infection simulation for individual-based model
            testingParams.nu = nu2;
            result2 = infectionIBM(diseaseParams, testingParams, modelParams, networkParams);
            colorVal2 = (result2.sigma2+shift)/(varMax+shift)*ones(1,3);
            allResults{count} = result2;
            count = count + 1;
            
            % Plot mean number of infectious people vs. daily testing
            % fraction
            figure(1)
            hold on
            plot(nu1, result1.meanInfEnd, symbsType2(m), 'MarkerSize', 8, 'MarkerFaceColor', colorVal1, 'MarkerEdgeColor', colorVal1, 'LineWidth', lw)
            plot(nu2, result2.meanInfEnd, symbsType2(m), 'MarkerSize', 8, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', colorVal2, 'LineWidth', lw)
            
            % Plot mean number of infectious people vs. lambda
            figure(2)
            hold on
            plot(result1.lambda, result1.meanInfEnd, symbsType2(m), 'MarkerSize', 8, 'MarkerFaceColor', colorVal1, 'MarkerEdgeColor', colorVal1, 'LineWidth', lw)
            plot(result2.lambda, result2.meanInfEnd, symbsType2(m), 'MarkerSize', 8, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', colorVal2, 'LineWidth', lw)
            
        end
    end
end

%% Data storage and figure formatting
% Save data
save('combinedInfections.mat')


% Label and format plots
figure(1)
hold on
xlabel('$$\nu =$$ daily testing fraction', 'Interpreter', 'latex', 'Color', 'k')
ylabel('$$\#$$ infectious people after 5 days', 'Interpreter', 'latex', 'Color', 'k')
xlim([0 1])
ylim([0 120])
set(gca, 'FontSize', 20)
set(gca, 'TickLabelInterpreter', 'latex')


figure(2)
hold on
xlabel('$$\lambda =$$ growth/decay parameter', 'Interpreter', 'latex', 'Color', 'k')
xlim([0 2])
ylim([0 120])
set(gca, 'FontSize', 20)
set(gca, 'TickLabelInterpreter', 'latex')


% Plot theory line and stability region
lambdaTest = 0:0.01:2;
resTest = modelParams.numInf0*lambdaTest.^(modelParams.Tend-1);
plot(lambdaTest, resTest, '-r', 'LineWidth', 1.5)

facesf = [1 2 3 4];
vertsv = [0 0; 1 0; 1 120; 0 120];

patch('Faces', facesf, 'Vertices', vertsv, 'FaceColor', 'black', 'FaceAlpha', 0.2, 'EdgeColor', 'black', 'EdgeAlpha', 0)
plot([0 2],[10 10], '--k', 'LineWidth', 1.5)

% Finish formatting plot
ax = gca;
ax.YTickLabel = [];
ax.YAxisLocation = 'right';
