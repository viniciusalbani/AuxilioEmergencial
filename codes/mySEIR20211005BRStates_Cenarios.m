clear all; clc; close all; format long e;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% We estimate the parameters of a SEIR-type epidemiological model by
%%%% using a maximum a posteriori estimator. All the estimation procedures
%%%% are carried out by LSQNONLIN, although the likelihood function (or
%%%% data misfit) is log-Poisson. The model parameters are estimated from
%%%% the daily records of infections and deaths.
%%%%-----------------------------------------------------------------------
%%%% This script generate scenarios in the case of absence of Auxilio
%%%% Emergencial.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Setup for the ODE solver and the least-square function
tol = 1.e-6;  % ode solver tolerance
options = odeset('AbsTol', tol,'RelTol',tol,'MaxOrder',5,'Stats',...
                                                         'off','Refine',1);
options2 = [];
options3 = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Experimental Data
%%% Loading data:
load('data_CORR_20211004.mat')
NSamples = 200;

%%%% Variables Concerning Infections and deaths
CASES2=zeros(length(dates),NSamples,size(States,1));
DEATHS2 = zeros(length(dates),NSamples,size(States,1));
MEDCASES = zeros(length(dates),size(States,1));
MEDDEATHS = zeros(length(dates),size(States,1));
CICASES = zeros(2,length(dates),size(States,1));
CIDEATHS = zeros(2,length(dates),size(States,1));

Reduction = (1:10)/100;% [1,5,10]/100; %
NCases = zeros(NSamples,size(States,1),length(Reduction));
MEDNCases = zeros(size(States,1),length(Reduction));
CINCases = zeros(2,size(States,1),length(Reduction));
NDeaths = zeros(NSamples,size(States,1),length(Reduction));
MEDNDeaths = zeros(size(States,1),length(Reduction));
CINDeaths = zeros(2,size(States,1),length(Reduction));

%%%% Variables concerning Social Isolation and the transmission Parameter
Mobility = zeros(413,27);
Mobility2 = Mobility;
H = [100 100 1000 400];
X = zeros(413,NSamples,size(States,1));
Y = zeros(413,NSamples,size(States,1));


for ss = 1:length(Reduction)
disp(['Reduction in SII: ',num2str(Reduction(ss))])
for zz = 1:size(States,1)

data = [Cases2(:,zz),Deaths2(:,zz)];
data = abs(data);
t_span2 = t_span(data(:,1)>0);
t_span3 = t_span2(1):dates(end);
data = data(data(:,1)>0,:);

t_actual = 0:size(data,1);

Mobility(:,zz) = auxE.data((zz-1)*413+1:zz*413,1);

%%%% Total Population:
N = Population(zz);      

%%%% Population proportion on each age range:
Proportion = ones;
PropInfections = sum(data(:,1))/N;
PropDeath = sum(data(:,2))/(0.4*sum(data(:,1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Initial parameters

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Estimating a time-dependent transmission parameter (beta) from day 1
%%%% until 19:
CASES3=zeros(length(t_actual),NSamples);
DEATHS3 = zeros(length(t_actual),NSamples);

parfor ll=1:NSamples
cases = CasesBoot(2:length(t_span2)+1,ll,zz);
deaths = DeathsBoot(2:length(t_span2)+1,ll,zz);
BETA = BETAStatesBoot(1:length(t_span2)+1,ll,zz);
BETA2 = BETA;

%%%% Smoothing data
Mobility2 = Mobility(:,zz);

%%%%%%
MobilityB = interp1(dates,Mobility2,t_span3)';
BETA3 = BETA2(DELAY(ll,zz):length(t_span3)+DELAY(ll,zz)-1);
[AA,idx] = sort(MobilityB);
BB = BETA3(idx);

yinit2 = yinit;
yb = zeros(length(t_actual),length(yinit));
yb(1,:) = yinit;
Beta = BETA;
for jj = 1:length(t_actual)-1
params2 = params;
if t_span(jj)>=datetime(2020,04,01) && t_span(jj)<datetime(2020,09,01)
Beta(jj+1) = interp1(AA,BB,Mobility2(jj)-Reduction(ss));
end
Beta(isnan(Beta)==1)=BETA(isnan(Beta)==1);
params2.beta = Beta(jj+1);
[~,y2] = ode45(@(t,y)seir_death_age_beta_b2(t,y, params2),t_actual(jj:jj+1),yinit2,options);
yinit2 = y2(end,:);
yb(jj+1,:) = yinit2;
end
factorD = zeros(length(t_actual),1);
factorDeath = params.factorDeath;
for jj = 1:length(t_actual)
factorD(jj) = factorDeath(t_actual(jj));
end  

%%%%% Final Number of Cases for each Age Range
sigma = params.sigma;
NewCases = sigma*yb(:,2)*N;
NewDeaths = params.Death*factorD.*yb(:,3)*N;
CASES3(:,ll)=NewCases;
DEATHS3(:,ll)=NewDeaths;

aux = interp1(t_span2,NewCases(2:end),t_span2(t_span2<datetime(2020,09,01)));
aux = max(aux,interp1(t_span2,data(:,1),t_span2(t_span2<datetime(2020,09,01))));
NCases(ll,zz,ss) = sum(aux);
aux = interp1(t_span2,NewDeaths(2:end),t_span2(t_span2<datetime(2020,09,01)));
aux = max(aux,interp1(t_span2,data(:,2),t_span2(t_span2<datetime(2020,09,01))));
NDeaths(ll,zz,ss) = sum(aux);
end

CASES2(1:length(t_actual),:,zz)=CASES3;
DEATHS2(1:length(t_actual),:,zz)=DEATHS3;

%%%%% Confidence Intervals Evaluation:
aux = sort(CASES2(1:length(t_actual),:,zz)');
MEDCASES(1:length(t_actual),zz) = median(aux)';
aux2 = round(0.05*NSamples);
aux = aux(aux2+1:end-aux2,:);
CICASES(:,1:length(t_actual),zz) = [min(aux);max(aux)];

aux = sort(DEATHS2(1:length(t_actual),:,zz)');
MEDDEATHS(1:length(t_actual),zz) = median(aux)';
aux2 = round(0.05*NSamples);
aux = aux(aux2+1:end-aux2,:);
CIDEATHS(:,1:length(t_actual),zz) = [min(aux);max(aux)];

aux = interp1(t_span2,data(:,1),t_span2(t_span2<datetime(2020,09,01)));
NC = sum(aux);
aux = interp1(t_span2,data(:,2),t_span2(t_span2<datetime(2020,09,01)));
ND = sum(aux);

aux = sort(NCases(:,zz,ss));
aux = aux(aux>=NC);
if isempty(aux)==1
    aux = sort(NCases(:,zz,ss));
end
MEDNCases(zz,ss) = max(NC,median(aux));
aux2 = round(0.05*length(aux));
aux = aux(aux2+1:end-aux2,:);
if isempty(aux)==1
    aux = sort(NCases(:,zz,ss));
end
CINCases(:,zz,ss) = max(NC,[min(aux);max(aux)]);

aux = sort(NDeaths(:,zz,ss));
aux = aux(aux>=ND);
if isempty(aux)==1
    aux = sort(NCases(:,zz,ss));
end
MEDNDeaths(zz,ss) = max(ND,median(aux));
aux2 = round(0.05*length(aux));
aux = aux(aux2+1:end-aux2,:);
if isempty(aux)==1
    aux = sort(NDeaths(:,zz,ss));
end
CINDeaths(:,zz,ss) = max(ND,[min(aux);max(aux)]);


disp([States(S(zz),:),';',num2str(round(NC)),';',num2str(round(MEDNCases(zz,ss))),...
';',num2str(round(CINCases(1,zz,ss))),';',num2str(round(CINCases(2,zz,ss))),';',...
    num2str(round(ND)),';',num2str(round(MEDNDeaths(zz,ss))),';',...
num2str(round(CINDeaths(1,zz,ss))),';',num2str(round(CINDeaths(2,zz,ss)))])
end
end
close all;clc;
H = [100,100,800,400];
NO = [1,3,4,14,21,22,27];
NE = [2,5,6,10,15,16,17,20,25];
CW = [7,9,12,13];
SE = [8,11,19,26];
SO = [18,23,24];

CasesObs = [527198,1125794,421417,1327970,393347,3795726];
DeathsObs = [13337,34641,8915,53711,8691,119295];

Med = zeros(6,length(Reduction));
CI1 = zeros(6,length(Reduction));
CI2 = zeros(6,length(Reduction));
for ss=1:length(Reduction)
Med(1,ss) = sum(MEDNCases(NO,ss));
Med(2,ss) = sum(MEDNCases(NE,ss));
Med(3,ss) = sum(MEDNCases(CW,ss));
Med(4,ss) = sum(MEDNCases(SE,ss));
Med(5,ss) = sum(MEDNCases(SO,ss));
Med(6,ss) = sum(MEDNCases(:,ss));

CI1(1,ss) = sum(CINCases(1,NO,ss));
CI1(2,ss) = sum(CINCases(1,NE,ss));
CI1(3,ss) = sum(CINCases(1,CW,ss));
CI1(4,ss) = sum(CINCases(1,SE,ss));
CI1(5,ss) = sum(CINCases(1,SO,ss));
CI1(6,ss) = sum(CINCases(1,:,ss));

CI2(1,ss) = sum(CINCases(2,NO,ss));
CI2(2,ss) = sum(CINCases(2,NE,ss));
CI2(3,ss) = sum(CINCases(2,CW,ss));
CI2(4,ss) = sum(CINCases(2,SE,ss));
CI2(5,ss) = sum(CINCases(2,SO,ss));
CI2(6,ss) = sum(CINCases(2,:,ss));
end


for jj=1:6

figure
hold on
grid on
box on
title('Cases')
area(100*Reduction,100*(CI2(jj,:)/CasesObs(jj)-ones),'linestyle',':',...
    'FaceColor',[255,160,122]/255,'FaceAlpha',0.5);
plot(100*Reduction,100*(Med(jj,:)/CasesObs(jj)-ones),'r','LineWidth',2)
area(100*Reduction,100*(CI1(jj,:)/CasesObs(jj)-ones),'linestyle',':','FaceColor',[1,1,1]);
xlim([1,10])
legend('90% CI','Median')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
ylabel('Increment in Cases (%)')
xlabel('Points Reduced in the Social Isolation Index')
hold off
saveas(gcf,'CasesIncrement.fig');
print('-dpng','CasesIncrement');


figure
hold on
grid on
box on
title('Cases')
area(100*Reduction,100*(CI2(jj,:)/CasesObs(jj)-ones),'linestyle',':',...
    'FaceColor',[255,160,122]/255,'FaceAlpha',0.5);%[51,236,255]/255);
plot(100*Reduction,100*(Med(jj,:)/CasesObs(jj)-ones),'r','LineWidth',2)
area(100*Reduction,100*(CI1(jj,:)/CasesObs(jj)-ones),'linestyle',':','FaceColor',[1,1,1]);
xlim([1,10])
legend('90% CI','Median')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
set ( gca, 'xdir', 'reverse' )
ylabel('Increment in Cases (%)')
xlabel('Points Reduced in the Social Isolation Index')
hold off
saveas(gcf,'CasesIncrementB.fig');
print('-dpng','CasesIncrementB');

disp(num2str([Reduction',100*(Med(jj,:)'/CasesObs(jj)-ones),100*(CI1(jj,:)'/CasesObs(jj)-ones),100*(CI2(jj,:)'/CasesObs(jj)-ones),]))
A = [(Med(jj,:)'/CasesObs(jj)-ones),(CI1(jj,:)'/CasesObs(jj)-ones),(CI2(jj,:)'/CasesObs(jj)-ones)];
end






for ss=1:length(Reduction)
Med(1,ss) = sum(MEDNDeaths(NO,ss));
Med(2,ss) = sum(MEDNDeaths(NE,ss));
Med(3,ss) = sum(MEDNDeaths(CW,ss));
Med(4,ss) = sum(MEDNDeaths(SE,ss));
Med(5,ss) = sum(MEDNDeaths(SO,ss));
Med(6,ss) = sum(MEDNDeaths(:,ss));

CI1(1,ss) = sum(CINDeaths(1,NO,ss));
CI1(2,ss) = sum(CINDeaths(1,NE,ss));
CI1(3,ss) = sum(CINDeaths(1,CW,ss));
CI1(4,ss) = sum(CINDeaths(1,SE,ss));
CI1(5,ss) = sum(CINDeaths(1,SO,ss));
CI1(6,ss) = sum(CINDeaths(1,:,ss));

CI2(1,ss) = sum(CINDeaths(2,NO,ss));
CI2(2,ss) = sum(CINDeaths(2,NE,ss));
CI2(3,ss) = sum(CINDeaths(2,CW,ss));
CI2(4,ss) = sum(CINDeaths(2,SE,ss));
CI2(5,ss) = sum(CINDeaths(2,SO,ss));
CI2(6,ss) = sum(CINDeaths(2,:,ss));
end


for jj=6:6

figure
hold on
grid on
box on
title('Deaths')
area(100*Reduction,100*(CI2(jj,:)/DeathsObs(jj)-ones),'linestyle',':','FaceColor',[255,160,122]/255,'FaceAlpha',0.5);%[51,236,255]/255);
plot(100*Reduction,100*(Med(jj,:)/DeathsObs(jj)-ones),'r','LineWidth',2)
area(100*Reduction,100*(CI1(jj,:)/DeathsObs(jj)-ones),'linestyle',':','FaceColor',[1,1,1]);
xlim([1,10])
legend('90% CI','Median')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
set ( gca, 'xdir', 'reverse' )
ylabel('Increment in Cases (%)')
xlabel('Points Reduced in the Social Isolation Index')
hold off
saveas(gcf,'DeathsIncrement.fig');
print('-dpng','DeathsIncrement');

figure
hold on
grid on
box on
title('Deaths')
area(100*Reduction,100*(CI2(jj,:)/DeathsObs(jj)-ones),'linestyle',':','FaceColor',[255,160,122]/255,'FaceAlpha',0.5);%[51,236,255]/255);
plot(100*Reduction,100*(Med(jj,:)/DeathsObs(jj)-ones),'r','LineWidth',2)
area(100*Reduction,100*(CI1(jj,:)/DeathsObs(jj)-ones),'linestyle',':','FaceColor',[1,1,1]);
xlim([1,10])
legend('90% CI','Median')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
ylabel('Increment in Cases (%)')
xlabel('Points Reduced in the Social Isolation Index')
hold off
saveas(gcf,'DeathsIncrementB.fig');
print('-dpng','DeathsIncrementB');

disp(num2str([Reduction',100*(Med(jj,:)'/DeathsObs(jj)-ones),100*(CI1(jj,:)'/DeathsObs(jj)-ones),100*(CI2(jj,:)'/DeathsObs(jj)-ones),]))
B = [(Med(jj,:)/DeathsObs(jj)-ones)',(CI1(jj,:)'/DeathsObs(jj)-ones),(CI2(jj,:)'/DeathsObs(jj)-ones)];
end