%% Test MCMC
% lclemens@uci.edu
% Creates simple two dimensional, two parameter function to test MCMC
% convergence and plots posterior distributions for parameters and
% correlations

%% Construct test function and test data
testFunction = @(p,t) ((t(:,1)-p(1)).^2 + (t(:,2)-p(2)).^2)/2;

testParams = [3 -1.5]; % true parameters
testDomain = [linspace(0,100,1000)',linspace(0,100,1000)']; 
testData = testFunction(testParams,testDomain);

% plot test function
figure; clf; hold on;
plot(testDomain,testData);

%% Run MCMC
% set initial parameters to start search from
initialparams = [3.1 -1.7];

% run and time MCMC
tic
% fitted parameters = MCMC(function, xdata, ydata, start parameters, lower bounds, upper bounds)
params_fit = simulateMCMC(testFunction,testDomain,testData,initialparams,[-100 -100],[100 100]);
toc

%% Plot posteriors and correlations
% VERY SLOW: calls subfunction once for every iteration for posterior
% construction
vParticles();

%% Plot data with fit
figure; clf; hold on;
plot(testDomain, testData,'.b');
plot(testDomain, testFunction(params_fit,testDomain),'r');
legend('Data','Fit');

