clear all; close all; clc;

%%% Evaluation of Correlation

load('data_BRStates_20210820.mat')
auxE = importdata('IIS_UF(1).csv');
dates = string(auxE.textdata(2:414,1));
dates = datetime(dates,'InputFormat','yyyy-MM-dd');
Mobility = zeros(413,27);
Mobility2 = Mobility;
DELAY = zeros(200,size(States,1));
H = [100 100 1000 400];
CORR = zeros(413,200,size(States,1));
LEN = zeros(200,size(States,1));
for jj = 1:size(States,1)
data = [Cases2(:,jj),Deaths2(:,jj)];
data = abs(data);
t_span2 = t_span(data(:,1)>0);
t_span3 = t_span2(1):dates(end);
data = data(data(:,1)>0,:);

t_actual = 0:size(data,1);
Mobility(:,jj) = auxE.data((jj-1)*413+1:jj*413,1);
for ll=1:200
cases = CasesBoot(2:length(t_span2)+1,ll,jj);
deaths = DeathsBoot(2:length(t_span2)+1,ll,jj);
Rt = R0StatesBoot(1:length(t_span2),ll,jj);
Rt2 = Rt;


Mobility2(:,jj) = Mobility(:,jj);
len = 6;
for ii = 1+len:length(Mobility2)%-(1+len) 
Mobility2(ii,jj) = mean(Mobility(ii-len:ii,jj));
Rt2(ii) = mean(Rt(ii-len:ii));
end
Dt = 3:30;
MobilityB = interp1(dates,Mobility2(:,jj),t_span3)';
DT = 30;
CORR1 = zeros(length(Dt),length(t_span3)-DT);
for zz = 1:length(Dt)
RtB =zeros(length(Dt(zz)+1:Dt(zz)+length(t_span3)),1);
RtB(2:end) = diff(Rt2(Dt(zz)+1:Dt(zz)+length(t_span3)));
for ii=1:length(t_span3)-DT
aux1 = zeros(size(MobilityB(ii:ii+DT)));
aux1(2:end) = diff(MobilityB(ii:ii+DT));
aux2 = corr([aux1,RtB(ii:ii+DT)]);
CORR1(zz,ii) = aux2(2,1);
end
end
SIGN1 = sign(CORR1);
SIGN1 = abs(min(0,SIGN1));
AUX1 = sum(SIGN1,2);
[~,AUX] = max(AUX1);
DELAY(ll,jj) = Dt(AUX);
CORR(1:length(CORR1(AUX,:)),ll,jj)=CORR1(AUX,:);
LEN(ll,jj) = length(CORR1(AUX,:));
end


%%% CORR:
aux = sort(CORR(1:max(LEN(:,jj)),:,jj)');
MECORR = median(aux);
aux2 = round(0.25*NSamples);
aux = aux(aux2+1:end-aux2,:);
CICORR = [min(aux);max(aux)];

%%% DELAY
aux = sort(DELAY(:,jj));
MEDELAY = median(aux);
aux2 = round(0.25*NSamples);
aux = aux(aux2+1:end-aux2,:);
CIDELAY = [min(aux);max(aux)];
disp([States(jj,:),num2str([MEDELAY,CIDELAY'])])

figure
hold on
grid on
box on
title('Correlation')
plot(t_span3(MEDELAY+1:max(LEN(:,jj))+MEDELAY),CICORR(2,:),'k');%[51,236,255]/255);
plot(t_span3(MEDELAY+1:max(LEN(:,jj))+MEDELAY),MECORR,'k','LineWidth',2);%[51,236,255]/255);
plot(t_span3(MEDELAY+1:max(LEN(:,jj))+MEDELAY),CICORR(1,:),'k');%[51,236,255]/255);
legend('50% CI','Median')
xlim([datetime(2020,03,01),datetime(2021,06,01)])
xticks(datas);
xtickformat('MMM')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
hold off
saveas(gcf,['Correlation',States(jj,:),'.fig']);
print('-dpng',['Correlation',States(jj,:)]);
end

A = median(Mobility(61:335,:))';
B = median(Mobility(336:end,:))';

aux = sort(Mobility(61:335,:));
aux2 = round(0.125*(335-60));
aux = aux(aux2+1:end-aux2,:);
CIMobility2020 = [min(aux);max(aux)];

aux = sort(Mobility(336:end,:));
aux2 = round(0.125*size(Mobility(336:end,:),1));
aux = aux(aux2+1:end-aux2,:);
CIMobility2021 = [min(aux);max(aux)];

AUX = [A,CIMobility2020',B,CIMobility2021'];
