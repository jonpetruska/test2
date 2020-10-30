clear all
close all

%% Give top level directory for data card
% This will probably be just /Volumes/Untitled
topdir = '/Volumes/G1C3'; % no need for final slash

delsac = false; % option to delete SAC files that are made in this dir

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
    data(1).yrs(iy,1) = str2num(tempy{iy});
    data(1).(['yyyy',tempy{iy}]) = struct();
end

%%% months
for iy = 1:Ny
    ystr = num2str(data.yrs(iy));
    
    temp = dir([topdir,'/',ystr]); %files within the year folder for year iy
    
    ismo = false(length(temp),1); % Creates a logical structure for months
    for ii = 1:length(temp) % loops over number of files in year iy
            if strcmp(temp(ii).name(1),'.'), continue; end %skips name with'.'
            if length(temp(ii).name)>2, continue; end % skips name longer than 2
            ismo(ii) = true; % Both ifs satisfied than ii ismo = 1
    end
    tempm = {temp(ismo).name}; %cell of unique months for year iy
    Nm = length(tempm);
    data.(['yyyy',ystr])(1).Nm = Nm;
    for im = 1:Nm
        data.(['yyyy',ystr])(1).mos(im,1) = str2num(tempm{im});
        data.(['yyyy',ystr])(1).(['mm',tempm{im}]) = struct();
    end

end

%%% SOH
for iy=1:Ny
    ystr = num2str(data.yrs(iy),'%04.f');
    Nm = data.(['yyyy',ystr])(1).Nm;
for im = 1:Nm
    mstr = num2str(data.(['yyyy',ystr]).mos(im),'%02.f');
    
    temp = dir([topdir,'/',ystr,'/',mstr]);
    isSOH = false(length(temp),1);
    for ii = 1:length(temp) % loops over number of files in year iy
            if strcmp(temp(ii).name(1),'.'), continue; end %skips name with'.'
            if length(temp(ii).name)~=3, continue; end % skips name longer than 2
            isSOH(ii) = true; % Both ifs satisfied than ii ismo = 1
    end
end
end
%%% days
for iy = 1:Ny
    ystr = num2str(data.yrs(iy),'%04.f');
    Nm = data.(['yyyy',ystr])(1).Nm;
for im = 1:Nm
    mstr = num2str(data.(['yyyy',ystr]).mos(im),'%02.f');
    
    temp = dir([topdir,'/',ystr,'/',mstr]);
    
    isda = false(length(temp),1);
    for ii = 1:length(temp)
            if strcmp(temp(ii).name(1),'.'), continue; end
            if ~any(regexp(temp(ii).name,'miniseed')), continue; end
            isda(ii) = true;
    end
    tempd = {temp(isda).name};
    tempp = {temp(isda).folder};
    Nd = length(tempd);
    data.(['yyyy',ystr])(1).(['mm',mstr])(1).Nd = Nd;
    data.(['yyyy',ystr])(1).(['mm',mstr])(1).datfiles = {tempd{:}}';
    data.(['yyyy',ystr])(1).(['mm',mstr])(1).fpath = {tempp{:}}';
    % parse miniseed file names
    mseedfiles_info = parse_miniseed_filename(tempd)';
    data.(['yyyy',ystr])(1).(['mm',mstr])(1).dom...
        = [mseedfiles_info.day]';
    data.(['yyyy',ystr])(1).(['mm',mstr])(1).serialstart...
        = [mseedfiles_info.serialstart]';
    % Get station info
    if iy== 1 && im == 1 
        data.sta = mseedfiles_info(end).sta;
        data.nwk = mseedfiles_info(end).nwk;
        data.centaur = mseedfiles_info(1).centaur_SN;
    end
    % get SOH files
    data.(['yyyy',ystr])(1).(['mm',mstr])(1).soh =...
        dir([topdir,'/',ystr,'/',mstr,'/soh/',data.nwk,'.',data.sta,'*.csv*']);
    
end; clear im
end; clear iy