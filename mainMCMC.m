%% Markov Chain Monte Carlo algorithm
% Developed for use with SPR MPDPDE Fitting algorithm


%% Test MCMC

testFunction = @(p,t) ((t-p(1)).^2 + (t-p(2)).^2)/2;

testParams = [3 -1.5];
testDomain = linspace(0,100,1000);
testData = testFunction(testParams,testDomain);

figure; clf; hold on;
plot(testDomain,testData);


initialparams = [1 -4];
params_fit = simulateMCMC(testFunction,testDomain,testData,initialparams,[-100 -100],[100 100]);

%% Initialize parameters

% Open data

% Set initial parameters



%% Call MCMC





%% Call function to display posterior distributions