function parameters = simulateMCMC(fun,tdata,ydata,initialparams,lb,ub)
%% simulateMCMC accepts input:
% function (function to call MCMC on - accepts parameters and tdata)
% tdata
% ydata - data to compare function output
% initialparams - parameters to start MCMC chain at
% lb - lower bounds for parameters
% ub - upper bounds for parameters

%% Set constants
NTMAX = 1e9;
NTCHECK = 3000;
ntNextStationarityCheck = 2*NTCHECK;
KSCRITICAL = 0.01;
NTADAPT = 200;


%% Initialize parameters
% initialize number of parameters
N_params = length(initialparams);

% set current parameters
parameters = initialparams;

% set lower, upper bounds based on input or defaults
if(exist('lb','var'))
    lowerbounds = lb;
else
    disp('No lower bounds specified.');
    lowerbounds = -Inf*ones(N_params,1);
end

if(exist('ub','var'))
    upperbounds = ub;
else
    disp('No upper bounds specified.');
    upperbounds = Inf*ones(N_params,1);
end

disp('lower bounds');
disp(lowerbounds);
disp('upper bounds');
disp(upperbounds);


% initialize step size for each parameter
dChi = ones(N_params,1);
proposals = zeros(N_params,1);
accepts = zeros(N_params,1);

% intialize stuff for convergence check
ntNextStationarityCheck = 3*NTCHECK;

% Histograms over lower/upper bound range
% function to bin parameters for posterior distribution
NBINS = 500;

binSize = (upperbounds - lowerbounds)./NBINS;

paramHistCounts = zeros(N_params,NBINS);
paramHistCountsPrevious = zeros(N_params,NBINS);


convergedTF=0;
nt=1;
E=Inf;

%% Find initial energy
% find initial fit
yCurrent = fun(parameters,tdata);
% find inital sum of squared residuals (SSR)
SSR = sum((yCurrent-ydata).^2);
E = SSR;

%% Propose new parameters

while(~convergedTF && nt < NTMAX)
    
    
    % adapt step size, but only until NTCHECK?
    if(mod(nt,NTADAPT)==0)
        rate = accepts./proposals;
        
        accepts(:) = 0;
        proposals(:) = 0;
        
        if(nt<NTCHECK)
            for iParam = 1:N_params
                if(rate(iParam) > 0.5 || rate(iParam) < 0.5)
                    dChi(iParam) = dChi(iParam)*rate(iParam)/0.44;
                    % possibly include max/min step sizes
                end
            end
        end
        
    end
%         /********* OUTLINE OF ALGORITHM *******************/
%         // 1. We create a new proposal configuration and
%         // then decide:
%         // 2. whether the proposal configuration satisfies any constraints (e.g., floor?), and
%         // 3. whether the proposal configuration is accepted by Metropolis test (always accept for a freely-jointed chain)
%         // After the configuration update, we:
%         // 4. collect any data of interest (e.g., ree; whether there is steric occlusion), and write it to file
%         // 5. test for stationarity (convergence)
%         // 6. increment the iteration and repeat

    % choose parameter to propose
    pPropose = ceil(N_params*rand);
    
    % pick step size
    dChiHere = dChi(pPropose);
    proposals(pPropose) = proposals(pPropose)+1;
    
    constraintSatisfiedTF = 0;
    
    while(~constraintSatisfiedTF) % no actual constraints right now

        % set initial configuration to current configuration
        paramsPropose = parameters;
        
        % update chosen parameter
        paramsPropose(pPropose) = paramsPropose(pPropose) + dChiHere*(2*rand-1);
        
        % Test constraints: Lower, upper bounds
        if(paramsPropose(pPropose) < lowerbounds(pPropose))
            paramsPropose(pPropose) = lowerbounds(pPropose);
        end
        if(paramsPropose(pPropose) > upperbounds(pPropose))
            paramsPropose(pPropose) = upperbounds(pPropose);
        end
       
        constraintSatisfiedTF=1;
    end
    
    %% Metropolis test
    
    % call function with proposed parameters
    yPropose = fun(paramsPropose,tdata);
    
    % calculate sum of square residuals (SSR)
    SSR = sum((yPropose-ydata).^2);
    
    % Compare new SSR to previous SSR
    ENew = SSR; % Calculate energy here
    if( rand < exp(E-ENew))
        % update energy
        E = ENew;
        
        % set parameters to proposed parameters
        parameters = paramsPropose;
        % set 
        yCurrent = yPropose;
        
        % increment accepted proposals
        accepts(pPropose) = accepts(pPropose)+1;
        
    end
    
    %% Collect data
    
    % collection of all parameter values
    parameters_all(:,nt) = parameters;


    % find bins for the current parameters
    binCurrent(:) = ceil( ( parameters(:) + abs(lowerbounds(:)) ) ./ binSize(:) );
    %disp(binCurrent);

    % Params histograms over bounds
    for p=1:N_params
        paramHistCounts(p,binCurrent(p)) = paramHistCounts(p,binCurrent(p))+1;
    end
    
    %% Test Stationarity
    
    % first append bins
    if(nt==NTCHECK)
        paramHistCountsPrevious = paramHistCountsPrevious + paramHistCounts;
        paramHistCounts = zeros(p,NBINS);
    end
    

    % test stationarity
    if(nt==ntNextStationarityCheck)

        % intialize
        ksStatistic = zeros(N_params,1);
        cdf1 = zeros(N_params,NBINS);
        cdf2 = zeros(N_params,NBINS);

        % compares first half to second half of data for each parameter
        for p=1:N_params
            % compute cumulative distribution functions
            cdf1(p,:) = cumsum(paramHistCountsPrevious(p,:))./(nt/2);
            cdf2(p,:) = cumsum(paramHistCounts(p,:))./(nt/2);
            % compute ksStatistic
            ksStatistic(p) = max(abs(cdf1(p,:)-cdf2(p,:)));
        end

        convergedTF = 1; % assume converged
        if(any(ksStatistic(:) >= KSCRITICAL))
            convergedTF = 0; % set to 0 if any of the parameters are not converged
        end
        
        % append bins
        paramHistCountsPrevious = paramHistCountsPrevious + paramHistCounts;
        paramHistCounts = zeros(p,NBINS);
        
        % set next stationarity checkpoint
        ntNextStationarityCheck = 2*ntNextStationarityCheck;

    end



    %% Iterate
    nt = nt+1;
end


    
disp(nt);
disp(parameters);

figure; clf; hold on;
plot(1:1:NBINS,paramHistCountsPrevious);
        
        
        
        
        
        
        
        









end


