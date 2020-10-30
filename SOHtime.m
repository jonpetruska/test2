function timedata=SOHtime(timelist)


[Dat] = importdata(timelist);

x = length(Dat.textdata(:,1))-1;
endtime = str2num(cell2mat(Dat.textdata(x,1)));
starttime = str2num(cell2mat(Dat.textdata(2,1)));

timedata.Longsecs = Dat.textdata(2:end,1);
timedata.UTC = Dat.textdata(2:end,2);

Loc = Dat.textdata(2:end,3);

%
%Extract Latitude data (if S then lat is negative)
%This section extracts Longitude (Negative if west)
%also got elevation
Nrec = length(Loc);
lat = zeros(Nrec,1);
lon = zeros(Nrec,1);
elev = zeros(Nrec,1);
%If error occurs in next line then csv_1 is not time data, check csv_2 
%Rename switch csv_1 and csv_2, Currently working on try/catch feature
for i = 1:Nrec 
    A = textscan(Loc{i},'%f%1s %f%1s %f%1s');
    lat(i) = A{1};
    if strcmp(A{2}{1},'S'), lat(i,1) = -lat(i,1); end
    lon(i) = A{3};
    if strcmp(A{4},'W'), lon(i,1) = -lon(i,1); end
   elev(i)= A{5};
end






%Coordinates used in Work
timedata.Lat = lat;
timedata.Lon = lon;
timedata.Elev = elev;


timedata.GNSS = Dat.textdata(2:end,4);

%%% Creates a double with entries of number of satellites
nSatellites = Dat.textdata(2:end,5);
timedata.Satellites = zeros(Nrec,1); 
for i = 1:Nrec
    A = textscan(nSatellites{i},'%f');
    timedata.Satellites(i) = A{1};
end 
%Takes the timing status and produces a double with entires 1 or -1, where -1
%is TIMING NOT OK and 1=TIMING OK
%For some reson gives time.Timing_Status in 1 row, we want 1 column.
%Also not sure what the not OK output is, might not convert correctly with
%textscan
Status = Dat.textdata(2:end,6);
timedata.Timing_Status = zeros(Nrec,1);
for i = 1:Nrec
    B = textscan(Status{i},'%s%s');
    if strcmp(B{2}, 'OK'), timedata.Timing_Status(i) = 1;
    else timedata.Timing_Status(i) = (0);
    end
    
end 


%Get some numerical answers out of the phase lock;
%Creates a double out of a cell. The entries in Phase_Lock have values
%0=Free running, 1=Coarse Lock, 2=Fine Lock
Lock= Dat.textdata(2:end,7);
timedata.Phase_Lock = zeros(Nrec,1);
for i = 1:Nrec
    A = textscan(Lock{i},'%s%s');
    if strcmp(A{1},'Free'), timedata.Phase_Lock(i) = 0; end
    if strcmp(A{1},'Coarse'), timedata.Phase_Lock(i) = 1;end 
    if strcmp(A{1},'Fine'), timedata.Phase_Lock(i) = 2; end
end 
        


timedata.Time_uncertainty = Dat.data(:,1);
timedata.DAC_Count = Dat.data(:,2);
timedata.Timing_Quality = Dat.data(:,3);
timedata.Timing_Error = Dat.data(:,4);

%assignin('base','SOHtime',time);

end