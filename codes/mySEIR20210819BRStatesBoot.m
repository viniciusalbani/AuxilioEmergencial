clear all; clc; close all; format long e;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% We estimate the parameters of a SEIR-type epidemiological model by
%%%% using a maximum a posteriori estimator. All the estimation procedures
%%%% are carried out by LSQNONLIN, although the likelihood function (or
%%%% data misfit) is log-Poisson. The model parameters are estimated from
%%%% the daily records of infections and deaths.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Setup for the ODE solver and the least-square function
tol = 1.e-6;  % ode solver tolerance
% set 'Stats','on' to get more info
% options = odeset('AbsTol', tol,'RelTol',tol,'MaxOrder',5,'Stats','on');

% note: set Refine switch to avoid interpolation
options = odeset('AbsTol', tol,'RelTol',tol,'MaxOrder',5,'Stats',...
                                                         'off','Refine',1);
options2 = [];%optimset('MaxFunEvals',10000,'MaxIter',7000,'TolFun',...
      %                                                 1e-30,'TolX',1e-30);
options3 = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Experimental Data

load('dataBrStates_20210511.mat')
States = ['AC';'AL';'AM';'AP';'BA';'CE';'DF';'ES';'GO';'MA';...
'MG';'MS';'MT';'PA';'PB';'PE';'PI';'PR';'RJ';'RN';'RO';'RR';'RS';'SC';...
'SE';'SP';'TO'];

S = 1:27;

NSamples = 200;

CasesBoot=zeros(length(dates),NSamples,length(S));
DeathsBoot = zeros(length(dates),NSamples,length(S));
R0StatesBoot = zeros(length(dates),NSamples,length(S));
BETAStatesBoot = zeros(length(dates),NSamples,length(S));
MEANR = zeros(16,length(S));
BETAStates = zeros(length(dates),length(S));
R = zeros(length(dates),length(S));
NCases = zeros(length(dates),length(S));
NDeaths = zeros(length(dates),length(S));


for zz = 1:length(S)
tic;
params = [];
data = [Cases2(:,S(zz)),Deaths2(:,S(zz))];
data = abs(data);
t_span2 = t_span(data(:,1)>0);
data = data(data(:,1)>0,:);
%%% We shall delete the last 10 days.
t_actual = 0:size(data,1);

%%%% Smoothing the data - averaging every 7e consecutive days:

%%%% Total Population:
N = Population(S(zz));      

%%%% Population proportion on each age range:
Proportion = ones;
PropInfections = sum(data(:,1))/N;
PropDeath = sum(data(:,2))/(0.4*sum(data(:,1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial parameters

I_0 = data(1,1);      % Potential initial infected and mild at t=0
E_0 = data(1,1);       % initial exposed
R_0 = 0;       % initial recovered 
D_0 = 0;       % initial dead

%--------------------------------------------------------------------------

%  params is a structure used to pass parameters to the
%   ODE solver

S_0 = N-(I_0+E_0+R_0+D_0);    % Suceptible pop.,  excluding initial infected 
params.N = N;  % N = total population

NumberOfAgeClasses = 1;  % total age ranges

yinit(1:NumberOfAgeClasses) = S_0*Proportion;
yinit(NumberOfAgeClasses+1:2*NumberOfAgeClasses) = E_0*Proportion;
yinit(2*NumberOfAgeClasses+1:3*NumberOfAgeClasses) = I_0*Proportion;
yinit(3*NumberOfAgeClasses+1:4*NumberOfAgeClasses) = R_0*Proportion;
yinit(4*NumberOfAgeClasses+1:5*NumberOfAgeClasses) = D_0*Proportion;
yinit = yinit/N;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Model Parameters

params.sigma = 1./5.1;   % inverse of mean incubation time
params.NumberOfAgeClasses = NumberOfAgeClasses;


nu = 1./14; % Mean time until recovery
gamma = PropDeath; % Mean time until death
p  = (1-gamma); % Proportion of individuals that will recover
params.p = p; % In Mild conditions

%--------------------------------------------------------------------------
%%% RATES

Recovery = (nu+gamma)*p; % Recovery Rate
params.Recovery = Recovery;

%%% Death Rate
Death = gamma;
params.Death = Death; % in ICU individuals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Correcting Hospitalization and Death Rates
Death = ones(size(t_actual));
Death(3:end) = min(1,data(2:end,2)./data(1:end-1,1))/gamma;
params.factorDeath = @(t)interp1(t_actual,Death,t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% A priori value for the transmission parameter
R_zero = 1.4*2.8/0.1782;
gamma = 1/18;

%--------------------------------------------------------------------------
beta = 2.2911*R_zero*gamma;
params.beta = beta; % Transmission parameter

paramsOld = params;
yinitOld = yinit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Estimating a time-dependent transmission parameter (beta)
BETA = zeros(length(t_actual),1);
BETA(1) = beta; 
unknowns0 = beta;
priors = unknowns0;
yinit2 = yinit;
yb = zeros(length(t_actual),length(yinit));
yb(1,:) = yinit;
R0 = zeros(1,length(t_actual)-1);
LB = 1E-3;
UB = 10;
for jj = 1:length(t_actual)-1

OF = @(unknowns)ObjFun_BetaWithoutAge(t_actual(jj:jj+1),params,data(jj,:),options,priors,yinit2,unknowns);
unknowns = lsqnonlin(OF,unknowns0,LB,UB,options2);
params.beta= unknowns(1);
priors = unknowns;
unknowns0 = unknowns;
BETA(jj+1,:) = unknowns;
[~,y2] = ode45(@(t,y)seir_death_age_beta_b2(t,y, params),t_actual(jj:jj+1),yinit2,options);
yinit2 = y2(end,:);
yb(jj+1,:) = yinit2;
R0(jj) = basic_reproduction_rate_beta(Proportion,params,BETA(jj+1),t_actual(jj+1));
end

factorD = zeros(length(t_actual),1);
factorDeath = params.factorDeath;
for jj = 1:length(t_actual)
factorD(jj) = factorDeath(t_actual(jj));
end  

%% Final Number of Cases for each Age Range
sigma = params.sigma;
NewCases = sigma*yb(:,2)*N;
NewDeaths = params.Death*factorD.*yb(:,3)*N;

NCases(1:length(NewCases),zz) = NewCases;
NDeaths(1:length(NewDeaths),zz) = NewDeaths;
R(1:length(R0),zz)=R0';
BETAStates(1:length(BETA),zz) = BETA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate scenarios by bootstrapping
Bootstrapping_20210819;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CasesBoot(1:size(NewCasesBoot,1),:,zz)=NewCasesBoot;
DeathsBoot(1:size(NewDeathsBoot,1),:,zz)=NewDeathsBoot;
R0StatesBoot(1:size(R0Boot,2),:,zz)=R0Boot';
BETAStatesBoot(1:size(BETABoot,2),:,zz)=BETABoot';

elt = toc;
disp(['Elapsed time: ',num2str(elt),' seconds'])
end
save data_BRStates_20210820;
close all; clc;
% S = [1,2,5,10,16,19,23,24,25,26];
X = -10:1:0;
IncrementCases = zeros(size(NCasesCI));
IncrementDeaths = zeros(size(NCasesCI));
IncrementCasesPerc = zeros(3,size(Increment,2),length(S));
IncrementDeathsPerc = zeros(3,size(Increment,2),length(S));
for jj=1:length(S)
figure
hold on
grid on
box on
title(['Cases in ',States(S(jj),:)])
plot(X,100*NCases(jj,:)/NCases(jj,end)-100,'r')
plot(X,100*NCasesCI(1,:,jj)/NCases(jj,end)-100,'k')
plot(X,100*NCasesCI(2,:,jj)/NCases(jj,end)-100,'k')
xlabel('Reduction in Social Isolation Index')
ylabel('Increment in Cases (%)')
hold off
figure
hold on
grid on
box on
title(['Deaths in ',States(S(jj),:)])
plot(X,100*NDeaths(jj,:)/NDeaths(jj,end)-100,'-sr')
plot(X,100*NDeathsCI(1,:,jj)/NDeaths(jj,end)-100,'k')
plot(X,100*NDeathsCI(2,:,jj)/NDeaths(jj,end)-100,'k')
xlabel('Reduction in Social Isolation Index')
ylabel('Increment in Deaths (%)')
hold off
IncrementCasesPerc(1,:,jj) = 100*NCasesCI(1,:,jj)/NCases(jj,end)-100;
IncrementCasesPerc(2,:,jj) = 100*NCasesCI(2,:,jj)/NCases(jj,end)-100;
IncrementCasesPerc(3,:,jj) = 100*NCases(jj,:)/NCases(jj,end)-100;
IncrementDeathsPerc(1,:,jj) = 100*NDeathsCI(1,:,jj)/NDeaths(jj,end)-100;
IncrementDeathsPerc(2,:,jj) = 100*NDeathsCI(2,:,jj)/NDeaths(jj,end)-100;
IncrementDeathsPerc(3,:,jj) = 100*NDeaths(jj,:)/NDeaths(jj,end)-100;
IncrementCases(:,:,jj) = round(NCasesCI(:,:,jj)-NCases(jj,end));
IncrementDeaths(:,:,jj) = round(NDeathsCI(:,:,jj)-NDeaths(jj,end));
disp(['State: ',States(S(jj),:)])
% disp(num2str([X(1),X(6),X(10)]))
disp(['Cases; ',num2str(IncrementCasesPerc(3,1,jj)),'% (',num2str(IncrementCasesPerc(1,1,jj)),'%--',num2str(IncrementCasesPerc(2,1,jj)),'%);',...
    num2str(IncrementCasesPerc(3,6,jj)),'% (',num2str(IncrementCasesPerc(1,6,jj)),'%--',num2str(IncrementCasesPerc(2,6,jj)),'%);',...
    num2str(IncrementCasesPerc(3,10,jj)),'% (',num2str(IncrementCasesPerc(1,10,jj)),'%--',num2str(IncrementCasesPerc(2,10,jj)),'%);'])
disp(['Deaths; ',num2str(IncrementDeathsPerc(3,1,jj)),'% (',num2str(IncrementDeathsPerc(1,1,jj)),'%--',num2str(IncrementDeathsPerc(2,1,jj)),'%);',...
    num2str(IncrementDeathsPerc(3,6,jj)),'% (',num2str(IncrementDeathsPerc(1,6,jj)),'%--',num2str(IncrementDeathsPerc(2,6,jj)),'%);',...
    num2str(IncrementDeathsPerc(3,10,jj)),'% (',num2str(IncrementDeathsPerc(1,10,jj)),'%--',num2str(IncrementDeathsPerc(2,10,jj)),'%);'])
end
%     