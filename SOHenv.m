function envdata = SOHenv(envlist)

%Extracts external voltages
[Dat] = importdata(envlist);
envdata.time = Dat.textdata(:,2);
envdata.Vol_A1 = Dat.data(~isnan(Dat.data(:,1)),1);
envdata.Vol_A2 = Dat.data(~isnan(Dat.data(:,2)),2);
envdata.Vol_A3 = Dat.data(~isnan(Dat.data(:,3)),3);
envdata.Ext_SOH1= Dat.data(~isnan(Dat.data(:,4)),4);
envdata.Ext_SOH2= Dat.data(~isnan(Dat.data(:,5)),5);
envdata.Ext_SOH3= Dat.data(~isnan(Dat.data(:,6)),6);

%assignin('base','SOHenv',env);
end