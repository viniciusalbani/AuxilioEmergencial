%%% Evaluation of the reproduction number using the next-generation matrix
%%% technique
function R0 = basic_reproduction_rate_beta(S,params,beta,t)
%%% S = susceptible individuals
%%% params = set of model parameters
%%% beta = transmission parameter
%%% t = current time

factorD = params.factorDeath;

Number = params.NumberOfAgeClasses;
mu = factorD(t).*params.Death;
nu = ones-mu;
sigma = params.sigma;

f = [0,beta*S;0,0];
v = diag([sigma,nu+mu]);
v = v - [zeros(1,2);[sigma,0]];
aux = f/v;
eigen = eig(aux);
R0 = max(eigen); % Time-dependent reproduction number