function instdata=SOHinst(instlist)

%Extracts the current, voltages and temp from the instrument
[Dat] = importdata(instlist);
instdata.time = Dat.textdata(:,2);
instdata.Supply_Voltage = Dat.data(:,1);
instdata.Current = Dat.data(:,2);
instdata.Temp = Dat.data(:,3);

%assignin('base','inst',inst);
end