%%% SEIR-type model solver
function  yprime = seir_death_age_beta_b2(t,y, params)
%%% t = current time (double)
%%% y = current values for the ODE system solution (double)
%%% params  = set of parameters of the ODE system (struct)

factorD = params.factorDeath; %%% daily death proportion

Number = params.NumberOfAgeClasses; 
beta = params.beta;
mu = factorD(t).*params.Death;
nu = ones-mu;
sigma = params.sigma;

S = y(1);
E = y(2);
I = y(3);
%      R = y(4);
%      D = y(5);

yprime = zeros(length(y),1);
for jj = 1:Number
yprime(1) = -S*(beta*I);
yprime(2) = S*(beta*I) - sigma*E;
yprime(3) = sigma*E-(nu+mu)*I;              
yprime(4) = nu*I;
yprime(5) = mu*I;
end