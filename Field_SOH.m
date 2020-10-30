clear all
close all
%% Code to analyze State of health Data from a Centaur Data Logger. 

%% Evan Marschall with Collaboration\Assistance from Zach Eilon and Jon Petruska (UCSB)
%% Give top level directory for data card
% Set the desired file path, end with /station name (ex. /SR.G1C3)
topdir = '/Volumes/EMPORT/STATE_OF_HEALTH/SR.G3H1'; % no need for final slash after station name

delsac = false; % option to delete SAC files that are made in this dir
profile on
addpath = 'Functions';
%% work out when there is data

%%% years
temp = dir(topdir);
%
isyr = false(length(temp),1); %Creates a logical structure for years.
for ii = 1:length(temp) %Runs through names from top directory
        if strcmp(temp(ii).name(1),'.'), continue; end %skips name with'.'
        if length(temp(ii).name)~=4, continue; end %skips names not length 4
        isyr(ii) = true; % if both ifs satisfied than entry is 1.
        %indicates that the entry ii corresponds to a year
end

tempy = {temp(isyr).name}; % creates a cell for unique yrs from logical isyr
Ny = length(tempy); % number of years
for iy = 1:Ny %initiates creating a data structure (empty) for each year
    SOH(1).yrs(iy,1) = str2num(tempy{iy});
    SOH(1).(['yyyy',tempy{iy}]) = struct();
end

%%% months
for iy = 1:Ny
    ystr = num2str(SOH.yrs(iy));
    
    temp = dir([topdir,'/',ystr]); %files within the year folder for year iy
    
    ismo = false(length(temp),1); % Creates a logical structure for months
    for ii = 1:length(temp) % loops over number of files in year iy
            if strcmp(temp(ii).name(1),'.'), continue; end %skips name with'.'
            if length(temp(ii).name)>2, continue; end % skips name longer than 2
            ismo(ii) = true; % Both ifs satisfied than ii ismo = 1
    end
    tempm = {temp(ismo).name}; %cell of unique months for year iy
    Nm = length(tempm);
    SOH.(['yyyy',ystr])(1).Nm = Nm;
    for im = 1:Nm
        SOH.(['yyyy',ystr]).mos(im,1) = str2num(tempm{im});
        SOH.(['yyyy',ystr]).(['mm',tempm{im}]) = struct();
    end

end

%%% Days
for iy=1:Ny
    ystr = num2str(SOH.yrs(iy),'%04.f');
    Nm = SOH.(['yyyy',ystr]).Nm;
for im = 1:Nm
    mstr = num2str(SOH.(['yyyy',ystr]).mos(im),'%02.f');
    temp = dir([topdir,'/',ystr,'/',mstr,'/','soh']);
    isCSV = false(length(temp),1);
    for ii = 1:length(temp) % loops over number of files in year iy
            if strcmp(temp(ii).name(1),'.'), continue; end %skips name with'.'
            isCSV(ii) = true; 
    end
    
    for ii = 1:length(temp) % loops over number of files in year iy
            if strcmp(temp(ii).name(1),'.'), continue; end %skips name with'.'
    
    tempn = {temp(isCSV).name};
    tempf = {temp(isCSV).folder};
    Nf = length(tempn);
    SOH.(['yyyy',ystr]).(['mm',mstr])(1).Nf = Nf;
    SOH.(['yyyy',ystr]).(['mm',mstr])(1).csvfiles = {tempn{:}}';
    SOH.(['yyyy',ystr]).(['mm',mstr])(1).fpath = {tempf{:}}';
  
    end
    
    % Days
   temp = dir([topdir,'/',ystr,'/',mstr,'/','soh']);
   isdy = false(length(temp),1);
   
   for ii = 1:length(temp) % loops over number of files in year iy
            if strcmp(temp(ii).name(1),'.'), continue; end %skips name with'.'
            if temp(ii).date(1:2) == temp(ii-1).date(1:2), continue;end
            isdy(ii) = true;
            
    end
    tempd = temp(isdy);
    days30 = [4 6 9 11];
    days31 = [1 3 5 7 8 10 12];
    if ismember((SOH.(['yyyy',ystr]).mos(im)),days30) == 1
        numdays = 30;
    elseif ismember((SOH.(['yyyy',ystr]).mos(im)),days31) == 1
        numdays = 31;
    elseif SOH.(['yyyy',ystr]).mos(im) == 2
        if isleap(SOH.yrs(iy))==1,numdays = 29;
        else, numdays = 29; 
        end
    end
    for id = 1:numdays
        ymdir = [topdir,'/',ystr,'/',mstr,'/','soh'];
        dstr = sprintf('%02.0f',id);
        SOH.(['yyyy',ystr]).(['mm',mstr])(1).days(id,1) = id;
        SOH.(['yyyy',ystr]).(['mm',mstr]).(['dd',dstr]) = struct();
        
        %% get rid of any hidden
        
        

        temp_csv = dir([ymdir,'/*',ystr,mstr,dstr,'*','.csv']);
            kill = [];for jj = 1:length(temp_csv), if strcmp(temp_csv(jj).name(1),'.'), kill(jj) = 1; else kill(jj) = 0; end; end
            temp_csv = temp_csv(~kill);
            csv = [temp_csv.folder,'/',temp_csv.name];
        
        temp_csv1 = dir([ymdir,'/*',ystr,mstr,dstr,'*','.csv_1']);
            kill1 = [];for jj = 1:length(temp_csv1), if strcmp(temp_csv1(jj).name(1),'.'), kill1(jj) = 1; else kill1(jj) = 0; end; end
            temp_csv1 = temp_csv1(~kill1);
            csv1 = [temp_csv1.folder,'/',temp_csv1.name];
        
        
        temp_csv2 = dir([ymdir,'/*',ystr,mstr,dstr,'*','.csv_2']);
            kill2 = [];for jj = 1:length(temp_csv2), if strcmp(temp_csv2(jj).name(1),'.'), kill2(jj) = 1; else kill2(jj) = 0; end; end
            temp_csv2 = temp_csv2(~kill2);
            csv2 = [temp_csv2.folder,'/',temp_csv2.name];
        
        if isempty(csv)==0
            try 
                csv_test = importdata(csv);%gets data from csv
                [m1,n1] = size(csv_test.textdata);%size cell 
                if n1 >= 9
                    timelist = csv;
                elseif  n1 == 5
                    instlist = csv;
                else
                    envlist = csv;
                end
            catch
                fprintf('Error')
            end
        end
        if isempty(csv1)==0
        try 
            csv_test1 = importdata(csv1);%gets data from csv
            [m2,n2] = size(csv_test1.textdata);%size cell 
        if n2 >= 9,timelist = csv1;
        elseif  n2 == 5, instlist = csv1;
        else, envlist = csv1;
        end
        catch
            fprintf('Error')
        end
        end
        try
        csv_test2 = importdata(csv2);%gets data from csv
        [m3,n3] = size(csv_test2.textdata);%size cell 
        if n3 >= 9,timelist = csv2;
        elseif  n3 == 5, instlist = csv2;
        else, envlist = csv2;
        end
        catch
            fprintf('Error')
        end
        
        
        %inst
        if exist('instlist')
        instdata = SOHinst(instlist);
            SOH.(['yyyy',ystr]).(['mm',mstr]).(['dd',dstr]).Supply_Voltage = instdata.Supply_Voltage; %Supply voltage
            SOH.(['yyyy',ystr]).(['mm',mstr]).(['dd',dstr]).Current = instdata.Current; %Current
            SOH.(['yyyy',ystr]).(['mm',mstr]).(['dd',dstr]).Temp = instdata.Temp; %Temperature
        else, continue
        end
        
        
        
        %times
        if exist('timelist') 
            timedata = SOHtime(timelist);
            SOH.(['yyyy',ystr]).(['mm',mstr]).(['dd',dstr]).timeUTC = timedata.UTC;%UTC from time data
            SOH.(['yyyy',ystr]).(['mm',mstr]).(['dd',dstr]).Lat = timedata.Lat; %latitudes
            SOH.(['yyyy',ystr]).(['mm',mstr]).(['dd',dstr]).Lon = timedata.Lon; %Longitudes
            SOH.(['yyyy',ystr]).(['mm',mstr]).(['dd',dstr]).Elev = timedata.Elev; %Elevation
            SOH.(['yyyy',ystr]).(['mm',mstr]).(['dd',dstr]).GNSS = timedata.GNSS; %GNSS
            SOH.(['yyyy',ystr]).(['mm',mstr]).(['dd',dstr]).Sat = timedata.Satellites; %Number satellites
            SOH.(['yyyy',ystr]).(['mm',mstr]).(['dd',dstr]).Time_Status = timedata.Timing_Status; %Timing Status
            SOH.(['yyyy',ystr]).(['mm',mstr]).(['dd',dstr]).Lock = timedata.Phase_Lock; %Phase Locking
            SOH.(['yyyy',ystr]).(['mm',mstr]).(['dd',dstr]).Time_Uncertainty = timedata.Time_uncertainty; %Uncertainty in Time
            SOH.(['yyyy',ystr]).(['mm',mstr]).(['dd',dstr]).Time_Quality = timedata.Timing_Quality; %time Quality
            SOH.(['yyyy',ystr]).(['mm',mstr]).(['dd',dstr]).Time_Error = timedata.Timing_Error; %Error in timing
        else, continue 
        end
        
        
        %env
        if exist('envlist')
        envdata = SOHenv(envlist);
            SOH.(['yyyy',ystr]).(['mm',mstr]).(['dd',dstr]).mass1 = envdata.Vol_A1; %mass 1
            SOH.(['yyyy',ystr]).(['mm',mstr]).(['dd',dstr]).mass2 = envdata.Vol_A2; %mass 2
            SOH.(['yyyy',ystr]).(['mm',mstr]).(['dd',dstr]).mass3 = envdata. Vol_A3; %mass 3
            SOH.(['yyyy',ystr]).(['mm',mstr]).(['dd',dstr]).Ext1 = envdata.Ext_SOH1; 
            SOH.(['yyyy',ystr]).(['mm',mstr]).(['dd',dstr]).Ext2 = envdata.Ext_SOH2; 
            SOH.(['yyyy',ystr]).(['mm',mstr]).(['dd',dstr]).Ext3 = envdata.Ext_SOH3;
        else, continue
        end
        
                % Clear all to ensure no double counting
                clear env
                clear inst
                clear times
                clear envlist
                clear timelist
                clear instlist
                

end % loop on days
    clear numdays
end % loop on months
end % loop on years
profile viewer


%% Concat
%Produces a large structure with all the data for the deployment
ALL_SOH = Concat_SOH(SOH);
ALL_SOH.name = fliplr(strtok(fliplr(topdir),'/'));  


% Create a field of starttimes

%% PLOT
figure1 = Plot_ALL_Day(ALL_SOH)
set(findall(gcf,'-property','FontSize'),'FontSize',14)
%Plot_allcoord(ALL_SOH);
%Plot_allinst(ALL_SOH);
%Plot_alltime(ALL_SOH);
%Plot_allGPS(ALL_SOH);