%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Reading Daily Cases and Deaths in Texas

aux = importdata('HIST_PAINEL_COVIDBR_11mai2021.csv');
States = ['AC';'AL';'AM';'AP';'BA';'CE';'DF';'ES';'GO';'MA';...
'MG';'MS';'MT';'PA';'PB';'PE';'PI';'PR';'RJ';'RN';'RO';'RR';'RS';'SC';'SE';'SP';'TO'];
States = string(States);
dailyCasesDeaths = aux.data(:,[4,6]);
datesCD = string(aux.textdata(2:end,8));
datesCD = datetime(datesCD,'InputFormat','yyyy-MM-dd');

Population = [894470
3351543
4207714
861773
14930634
9187103
3055149
4064052
7113540
7114598
21292666
2809394
3526220
8690745
4039277
9616621
3281480
11516840
17366189
3534165
1796460
631181
11422973
7252502
2318822
46289333
1590248];
for ll = 1:27
ii = ones;
for jj = 1:size(datesCD,1)
if string(aux.textdata(jj+1,2))==States(ll) && string(aux.textdata(jj+1,3)) == string(aux.textdata(2,3))
Cases(ii,ll) = dailyCasesDeaths(jj,1);    
Deaths(ii,ll) = dailyCasesDeaths(jj,2);
if ll == 1
dates(ii) = datesCD(jj);
end
ii = ii+1;
end
end
end
Cases = Cases(1:length(dates),:);
Deaths = Deaths(1:length(dates),:);
t_span = dates;
Cases2 = Cases;
Deaths2 = Deaths;
for jj=7:size(Cases2,1)
for ii = 1:size(Cases2,2)
Cases2(jj,ii) = mean(Cases(jj-6:jj,ii));
Deaths2(jj,ii) = mean(Deaths(jj-6:jj,ii));
end
end

save('dataBrStates_20210511.mat')

CasesBy100K = 1E5*(sum(Cases(1:311,:))./Population')';
CasesBy100K = [CasesBy100K,1E5*(sum(Cases(312:end,:))./Population')'];
DeathsBy100K = 1E5*(sum(Deaths(1:311,:))./Population')';
DeathsBy100K = [DeathsBy100K,1E5*(sum(Deaths(312:end,:))./Population')'];