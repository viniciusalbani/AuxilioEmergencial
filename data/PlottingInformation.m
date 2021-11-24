% close all;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importing the number of individuals receiving the Auxilio Emergencial

aux = importdata('dados_beneficiario.txt');
ProportionPop = aux.data(:,12:end);
datesCD = string(aux.textdata(1,13:end));
datesCD = datetime(datesCD,'InputFormat','dd/MM/yyyy');

States = ['AC';'AL';'AM';'AP';'BA';'CE';'DF';'ES';'GO';'MA';...
'MG';'MS';'MT';'PA';'PB';'PE';'PI';'PR';'RJ';'RN';'RO';'RR';'RS';'SC';'SE';'SP';'TO'];
% % States = string(States);
H = [100 100 1000 400];
datas = [datetime(2020,02:12,01),datetime(2021,1:6,01)];

%%% Population Proportion receiving the Auxilio Emergencial:
for zz = 1:size(States,1)
figure
hold on
grid on
box on
% title(States(zz,:))
title('Auxilio Emergencial')
bar(datesCD+15,100*ProportionPop(zz,:),'FaceAlpha',0.5)%,'FaceColor',[0 0.75 0.75],'EdgeColor','none')
ylabel('Population (%)')
xlim([datetime(2020,03,01),datetime(2021,06,01)])
ylim([0,40])
xticks(datas);
xtickformat('MMM')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
hold off
saveas(gcf,['AuxilioPago',States(zz,:),'.fig']);
print('-dpng',['AuxilioPago',States(zz,:)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Importing the Social Isolation Index

close all;
aux = importdata('IIS_UF(1).csv');
dates = string(aux.textdata(2:414,1));
dates = datetime(dates,'InputFormat','yyyy-MM-dd');
Mobility = zeros(413,27);
Mobility2 = Mobility;
MeanMobility = zeros(length(datas)-3,27);

for jj = 1:size(States,1)
Mobility(:,jj) = aux.data((jj-1)*413+1:jj*413,1);
Mobility2(:,jj) = Mobility(:,jj);
len = 6;
for ii = 1+len:length(Mobility2)%-(1+len) 
Mobility2(ii,jj) = mean(Mobility(ii-len:ii,jj));
end
for ii = 1:length(datas)-3
if ii < length(datas)-3
MeanMobility(ii,jj) = median(Mobility((dates>=datas(ii))&(dates<datas(ii+1)),jj));
else
MeanMobility(ii,jj) = median(Mobility(dates>=datas(ii),jj));
end
end

figure
hold on
grid on
box on
% title(States(jj,:))
title('Social Isolation Index')
% bar(dates,100*Mobility(:,jj),'FaceColor',[0 0.75 0.75],'EdgeColor','none')
bar(datas(1:end-3)+15,100*MeanMobility(:,jj),'r','FaceAlpha',0.8)%'FaceColor',[0 0.75 0.75])
area(dates,100*Mobility2(:,jj),'FaceColor',[255,160,122]/255,'FaceAlpha',0.5);%[51,236,255]/255);
legend('Monthly Median','Daily Data')
ylabel('Population (%)')
xlim([datetime(2020,03,01),datetime(2021,06,01)])
ylim([25,65])
xticks(datas);
xtickformat('MMM')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
hold off
saveas(gcf,['Mobility',States(jj,:),'.fig']);
print('-dpng',['Mobility',States(jj,:)]);
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


%AUX = [median(Mobility(61:335,:))',median(Mobility(335:end,:))'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 
% close all;
aux = importdata('dados_beneficiario2.csv');
ProportionPop = aux.data(:,1:end);
datesCD = string(aux.textdata(1,2:end));
datesCD = datetime(datesCD,'InputFormat','dd/MM/yyyy');

%%% Amount paid by the Auxilio Emergencial (Mean Value):
H = [100 100 1000 400];
datas = [datetime(2020,02:12,01),datetime(2021,1:6,01)];
for zz = 1:size(States,1)
figure
hold on
grid on
box on
% title(States(zz,:))
title('Auxilio Emergencial')
bar(datesCD+15,[ProportionPop(zz,1:end-1),0],'FaceAlpha',0.5,'FaceColor',[0 0.75 0.75])%,'FaceColor',[0 0.75 0.75],'EdgeColor','none')
ylabel('Mean Value (BRL)')
xlim([datetime(2020,03,01),datetime(2021,06,01)])
ylim([0,1000])
xticks(datas);
xtickformat('MMM')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
hold off
saveas(gcf,['AuxilioPagoMean',States(zz,:),'.fig']);
print('-dpng',['AuxilioPagoMean',States(zz,:)]);
end