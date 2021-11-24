%%% Objective function
function f = ObjFun_BetaWithoutAge(t_actual,params,data,options,priors,yinit,unknowns)
%%% t_actual = actual time (double)
%%% params = model parameters (struct)
%%% data = COVID-19 reports (double)
%%% options = options for the ODE solver
%%% priors = a priori set of parameters (double)
%%% yinit = initial condition for the ODE system (double)
%%% unknowns = set of unknown variables (double)

Number = params.NumberOfAgeClasses;
params.beta = unknowns(1);
tspan = [t_actual(1),t_actual(end)];

N = params.N;


sigma = params.sigma;

%%% Model prediction evaluation
[~,y]=ode45(@(t,y)seir_death_age_beta_b2(t,y, params),tspan,yinit,options);
NewInfections = sigma*y(end,Number+1:2*Number)*N;

% Gaussian Misfit or Likelihood
% f = (N*NewInfections-data(:,1));
% f = [f;10*(N*NewDeaths-data(:,2))];
% f = [f;(unknowns-priors)'];

% log-Poisson Misfit of Likelihood
Stirling = 0.5*log(2*pi*data(:,1)) + data(:,1).*log(data(:,1)) - data(:,1);
f = data(:,1).*log(NewInfections) - NewInfections - Stirling;
f = [f;1E-3*(unknowns-priors)'];
f(isnan(f)~=0)=zeros;