clear all; clc; close all; format long e;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% We estimate the parameters of a SEIR-type epidemiological model by
%%%% using a maximum a posteriori estimator. All the estimation procedures
%%%% are carried out by LSQNONLIN, although the likelihood function (or
%%%% data misfit) is log-Poisson. The model parameters are estimated from
%%%% the daily records of infections and deaths.
%%%%-----------------------------------------------------------------------
%%%% This code generates the figures in the supplementary materials. The
%%%% reproduction number is evaluated using the next-generation matrix
%%%% technique
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Setup for the ODE solver and the least-square function
tol = 1.e-6;  % ode solver tolerance
options = odeset('AbsTol', tol,'RelTol',tol,'MaxOrder',5,'Stats',...
                                                         'off','Refine',1);
options2 = [];
options3 = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Experimental Data

load('dataBrStates_20210511.mat')
States = ['AC';'AL';'AM';'AP';'BA';'CE';'DF';'ES';'GO';'MA';'MG';'MS';...
'MT';'PA';'PB';'PE';'PI';'PR';'RJ';'RN';'RO';'RR';'RS';'SC';'SE';'SP';'TO'];

CASES=zeros(length(dates),27);
DEATHS = zeros(length(dates),27);
R = zeros(length(dates),27);
MEANR = zeros(16,27);
BETAStates = zeros(length(dates),27);

for zz = 1:size(States,1)
params = [];
%%% Selecting the datasets for each State: daily numbers of cases and time-
%%% interval,
data = [Cases2(:,zz),Deaths2(:,zz)];
data = abs(data);
t_span2 = t_span(data(:,1)>0); % calendar time
data = data(data(:,1)>0,:); % nonzero COVID-19 reports
t_actual = 0:size(data,1); % time
%%%% We use 7-day moving averaged data:

%%%% Total Population:
N = Population(zz);      

Proportion = ones; % Change if considering age-ranges or spatial-distribution
PropInfections = sum(data(:,1))/N;
PropDeath = sum(data(:,2))/(0.4*sum(data(:,1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial parameters

I_0 = data(1,1);  % Potential initial infected and mild at t=0
E_0 = data(1,1);  % initial exposed
R_0 = 0;       % initial recovered 
D_0 = 0;       % initial dead

%--------------------------------------------------------------------------

%  params is a structure used to pass parameters to the
%   ODE solver

S_0 = N-(I_0+E_0+R_0+D_0); % Suceptible pop.,  excluding initial infected 
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
gamma = PropDeath; % observed death proportion
p  = (1-gamma); % observed proportion of individuals that will recover
params.p = p; % In Mild conditions

%--------------------------------------------------------------------------
%%% RATES

Recovery = (nu+gamma)*p; % Recovery Rate
params.Recovery = Recovery;

%%% Death Rate
Death = gamma;
params.Death = Death; % in ICU individuals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Correcting Death Rate using daily reports
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

paramsOld = params; % saving the initial set of parameters
yinitOld = yinit;   % saving the initial condition

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Estimating the time-dependent transmission parameter (beta) 
BETA = zeros(length(t_actual),1); % time-dependent transmission parameter
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
%%% Defining the objective function:
OF = @(unknowns)ObjFun_BetaWithoutAge(t_actual(jj:jj+1),params,...
                                data(jj,:),options,priors,yinit2,unknowns);
%%% Minimization:
unknowns = lsqnonlin(OF,unknowns0,LB,UB,options2);
params.beta= unknowns(1);
priors = unknowns;
unknowns0 = unknowns;
BETA(jj+1,:) = unknowns;

%%% Generating model predictions with the calibrated parameters:
[~,y2] = ode45(@(t,y)seir_death_age_beta_b2(t,y, params),...
                                         t_actual(jj:jj+1),yinit2,options);
yinit2 = y2(end,:);
yb(jj+1,:) = yinit2;
disp(num2str(unknowns))

%%% Evaluating the time-dependent reproduction number:
R0(jj) = basic_reproduction_rate_beta(Proportion,params,BETA(jj+1),...
                                                           t_actual(jj+1));
end

%%% Evaluating the model predictions of numbers of infections and deaths: 
factorD = zeros(length(t_actual),1);
factorDeath = params.factorDeath;
for jj = 1:length(t_actual)
factorD(jj) = factorDeath(t_actual(jj));
end  

%%% Final Number of Cases for each Age Range
sigma = params.sigma;
NewCases = sigma*yb(:,2)'*N;
NewDeaths = params.Death*factorD.*yb(:,3);

disp('  Infections   Deaths  ')
disp(num2str(round([sum(NewCases),sum(NewDeaths)])))

%%% Smoothing the reproduction number to see its behavior:
R0t = R0;
len = 3;
for jj = 1+len:length(t_actual)-(1+len) 
R0t(jj) = mean(R0(jj-len:jj+len));
end

H = [100 100 1000 400];
datas = [datetime(2020,02:12,01),datetime(2021,1:6,01)];

%%% plotting Results:
figure
hold on
grid on
box on
title('Daily Infections')
bar(t_span2,data(:,1),'FaceColor',[0 0.75 0.75],'EdgeColor','none','FaceAlpha',0.5)
plot(t_span2,NewCases(2:end),'r','LineWidth',1)
legend('Reported','Estimated','Location','NorthWest')
ylabel('Number of Individuals')
xlim([datetime(2020,03,01),datetime(2021,06,01)])
xticks(datas);
xtickformat('MMM')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
hold off
saveas(gcf,['Infections',States(zz,:),'.fig']);
print('-dpng',['Infections',States(zz,:)]);



figure
hold on
grid on
box on
title('Daily Deaths')
bar(t_span2,data(:,2),'FaceColor',[0 0.75 0.75],'EdgeColor','none','FaceAlpha',0.5)
plot(t_span2,NewDeaths(2:end)*N,'r','LineWidth',1)
legend('Reported','Estimated','Location','NorthWest')
ylabel('Number of Individuals')
xlim([datetime(2020,03,01),datetime(2021,06,01)])
xticks(datas);
xtickformat('MMM')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
hold off
saveas(gcf,['Deaths',States(zz,:),'.fig']);
print('-dpng',['Deaths',States(zz,:)]);

%%% R0

figure
hold on
grid on
box on
title('Reproduction Number') % Without smoothing
plot(t_span2,R0,'b','LineWidth',1)
plot(t_span2,ones(size(t_span2)),'k')
ylabel('Basic Reproduction Rate')
xlim([datetime(2020,03,01),datetime(2021,06,01)])
ylim([0,3])
xticks(datas);
xtickformat('MMM')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
hold off
saveas(gcf,['R0',States(zz,:),'.fig']);
print('-dpng',['R0',States(zz,:)]);


%%% Evaluating monthly medians
MeanR0 = zeros(1,length(datas)-1);
for ii = 1:length(datas)-1
if ii < length(datas)-1
MeanR0(ii) = median(R0((t_span2>=datas(ii))&(t_span2<datas(ii+1))));
else
MeanR0(ii) = median(R0(t_span2>=datas(ii)));
end
end

figure
hold on
grid on
box on
title('Effective Reproduction Number')
bar(datas(1:end-1)+15,MeanR0,'r','FaceAlpha',0.3)
plot(t_span2,R0t,'r','LineWidth',2)
plot(t_span2,ones(size(t_span2)),'k')
ylabel('Reproduction Number')
xlim([datetime(2020,03,01),datetime(2021,06,01)])
ylim([0,3])
xticks(datas);
xtickformat('MMM')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
hold off
saveas(gcf,['R0smooth',States(zz,:),'.fig']);
print('-dpng',['R0smooth',States(zz,:)]);

CASES(1:length(NewCases),zz)=NewCases;
DEATHS(1:length(NewDeaths),zz)=NewDeaths;
R(1:length(R0),zz)=R0';
MEANR(:,zz) = MeanR0';
BETAStates(1:length(BETA),zz) = BETA;
end