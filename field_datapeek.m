clear all
close all

%% Give top level directory for data card
% This will probably be just /Volumes/Untitled
%topdir = '/Volumes/G4H1'; % no need for final slash
topdir = '/Users/jonpetruska/Desktop/Huddle_Tests/H4_Data_SOH/HTH1';
%Huddle3SOH/H3.H3H1.centaur-3_4311_20181129_000000.miniseed
%topdir = '/Users/jonpetruska/Desktop/Huddle_Tests/Huddle3';

delsac = false; % option to delete SAC files that are made in this dir

%% give station coordinates
slat = 35.69534;
slon = -120.04059;

%% Define event to look at
%elat = -117.59933;
%elon = 35.7695;
%edep = 8;
%evtimestr = '2019-07-04 17:33:49'; %time in local time
%filtfs = [0.01 2]



elon = -127.8761;
elat = 43.5436;
edep = 5.4;
evtimestr = '2019-02-14 15:07:58'; %UTC
%evtimestr = '2019-08-28 09:07:58'; %local
filtfs = [0.05 1];%teleseismic

viewwind = [-100 2000];% default window to view in seconds from event time
bigwind = [-500 3500]; % window to subset in seconds from event time

%% preamble
%addpath('datprocfunctions');
spd = 24*3600;
data = struct();

%% work out when there is data
% years
temp = dir([topdir,'/*.miniseed']);

isyr = false(length(temp),1);
for ii = 1:length(temp)
        if strcmp(temp(ii).name(1),'.'), continue; end
        if length(temp(ii).name)~=4, continue; end
        isyr(ii) = true;
end
tempy = {temp(isyr).name};
Ny = length(tempy);
for iy = 1:Ny
    data(1).yrs(iy,1) = str2num(tempy{iy});
    data(1).(['yyyy',tempy{iy}]) = struct();
end


% months
for iy = 1:Ny
    ystr = num2str(data.yrs(iy));
    
    temp = dir([topdir,'/',ystr]);
    
    ismo = false(length(temp),1);
    for ii = 1:length(temp)
            if strcmp(temp(ii).name(1),'.'), continue; end
            if length(temp(ii).name)>2, continue; end
            ismo(ii) = true;
    end
    tempm = {temp(ismo).name};
    Nm = length(tempm);
    data.(['yyyy',ystr])(1).Nm = Nm;
    for im = 1:Nm
        data.(['yyyy',ystr])(1).mos(im,1) = str2num(tempm{im});
        data.(['yyyy',ystr])(1).(['mm',tempm{im}]) = struct();
    end

end
    
% days
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

%% parse into allfiles
Ndf = 0; % number of data files
dfs = {}; % names of data files
dps = {}; % paths to data files
serialstarts = []; % serial start times for each file
for iy = 1:Ny
    ystr = num2str(data.yrs(iy),'%04.f');

    for im = 1:data.(['yyyy',ystr])(1).Nm
        mstr = num2str(data.(['yyyy',ystr]).mos(im),'%02.f');
        modat = data.(['yyyy',ystr])(1).(['mm',mstr]);
        Ndf = Ndf + length(modat.datfiles);
        dfs = cat(1,dfs,modat.datfiles);
        dps = cat(1,dps,modat.fpath);
        serialstarts = cat(1,serialstarts,modat.serialstart);
    end; clear im
end; clear iy

%% Predict EQ times
[gcarc,seaz] = distance(slat,slon,elat,elon);
addpath /Users/jonpetruska/Documents/MATLAB/seizmo-master/mattaup/
tpt = tauptime('deg',gcarc,'dep',edep);
evtime = datenum(evtimestr);
arrtime = tpt(1).time/spd + evtime;
t0 = arrtime + bigwind(1)/spd; % start of window to request 
t1 = arrtime + bigwind(2)/spd; % end of window to request (surface waves gone at least 180deg) 
% find data file(s)
fid0 = find(serialstarts<t0,1,'last');
fid1 = find(serialstarts<t1,1,'last');
fid2read = fid0:fid1;
for idf = 1:length(fid2read)
    dfile = [dps{fid2read(idf)},'/',dfs{fid2read(idf)}];
    % Use mseed2sac to make a sac file in the current directory
    system(sprintf('/usr/local/bin/mseed2sac -O %s',dfile));
    sacfiles = dir([data.nwk,'.',data.sta,'*.SAC']);
    sfileinfo = parse_sac_filename({sacfiles.name}');
    
    for icsf = 1:length(sfileinfo)
        % Uses sacpc2mat to extract data + metadata from SAC file
        [SACdata(icsf),seisdat{icsf},~]=sacpc2mat(sacfiles(icsf).name);
        if delsac
            fprintf('Deleting SAC file %s\n',sacfiles(icsf).name);
            delete(sacfiles(icsf).name);
        end
        [ST(icsf),fstarttime(icsf)] = parse_timing(SACdata(icsf).event);
        dt = round(SACdata(icsf).times.delta,4);
        tt(:,icsf) = fstarttime(icsf) + [0:(SACdata(icsf).trcLen-1)]'*dt./spd;
    end
    % find which chan is which (search second char in chan name + chstr)
    icn = find(contains({sfileinfo.chan},[sfileinfo(1).chan(2),'N'])); 
    ice = find(contains({sfileinfo.chan},[sfileinfo(1).chan(2),'E']));
    icz = find(contains({sfileinfo.chan},[sfileinfo(1).chan(2),'Z']));
    % insert data as ENZ
    dat = [seisdat{ice},seisdat{icn},seisdat{icz}];
   

end

%% isolate data around event
inwin = ((tt(:,1)-evtime)*spd >= bigwind(1)) & ((tt(:,1)-evtime)*spd < bigwind(2));
ttw = tt(inwin,:);
datw = dat(inwin,:);

%% Get radial and transverse chans
datZRT = zne2zrt( datw(:,[3,2,1]), seaz );
datw = [datw,datZRT(:,2:3)];

% %% spectrogram
% figure(86);
% spectrogram(datw(:,3),'yaxis')


%% plot event
cont=0;
ampf = 1;
replot = 1;
chans = {'E','N','Z','R','T'};
plotchans = [1,2,3];
npoles = 2;

figure(88); clf; set(gcf,'pos',[58 77 1259 701]); 
ax1 = axes(gcf,'pos',[0.05 0.07 0.9 0.85]); hold on
title(sprintf('Event on %s,   \\Delta=%.1f?,   seaz = %.0f?',evtimestr,gcarc,seaz),'fontsize',22)

%% CLEAN DATA + window to window time
pretime = (evtime-tt(1,1))*spd;
cp = struct('samprate',1./dt,'pretime',-bigwind(1),...
        'prex',-bigwind(1),'postx',bigwind(2),...
        'taperx',0.02,'fhi',filtfs(2),'flo',filtfs(1),'npoles',npoles,'norm',1);
[ datwf,~,~,~,~,ttws,~ ] = data_clean( datw,cp );

while cont==0

    normf = max(max(abs(datwf)));
    
    % plot data
    if replot
        try delete(hb); end
        for ic = 1:3
            hb(ic,1) = plot(ax1,ttws,ampf*datwf(:,plotchans(ic))/normf + ic,'b');
        end
        % plot predicted arrivals
        for ip = 1:length(tpt)
                plot(ax1,tpt(ip).time*[1 1],[0 4],':','linewidth',0.5);
                text(ax1,tpt(ip).time + 2,0.2,tpt(ip).phase,'fontsize',14);
        end
    end
    replot = 0;
    set(ax1,'box','on','fontsize',16,'xlim',viewwind,'ylim',[0 4],'fontsize',15);
    % channel names
    set(ax1.YAxis,'TickValues',[0:0.5:4],...
        'TickLabels',{'','',chans{plotchans(1)},'',chans{plotchans(2)},'',chans{plotchans(3)},'',''},...
        'fontsize',22,'fontweight','bold')

    %% Interact with figure
    [x,y,b] = ginput(1);

        %% ================ PROCESS THE KEYBOARD INPUT ================
        switch b 

            case 'x'  % select window
                x1=x;
                hp(1) = plot(ax1,[x1,x1],axlim(ax1,[3,4]),'--r','Linewidth',1.5);
                [x2,y2] = ginput(1);
                hp(2) = plot(ax1,[x2,x2],axlim(ax1,[3,4]),'--r','Linewidth',1.5);
                viewwind = sort(round_level([x1;x2],dt));
                drawnow; pause(0.3)
                delete(hp);

            case 'f' % change filter

                flt = inputdlg({'Enter  max  period','Enter  min  period','Enter  N  filter  poles'},...
                               'Filter details',1,{num2str(1/filtfs(1)),num2str(1/filtfs(2)),num2str(npoles)});
                flt = str2num(char(flt));
                filtfs(1) = 1/flt(1); % set high-pass (low-f corner)
                filtfs(2) = 1/flt(2); % set low-pass (high-f corner)
                npoles = flt(3);      % set N poles

                cp = struct('samprate',1./dt,'pretime',-bigwind(1),...
                        'prex',-bigwind(1),'postx',bigwind(2),...
                        'taperx',0.02,'fhi',filtfs(2),'flo',filtfs(1),'npoles',npoles,'norm',1);
                [ datwf,~,~,~,~,ttws,~ ] = data_clean( datw,cp );
                
                replot=1;

            case 'r'  % toggle ZNE-ZRT (from top down

                if isequal(plotchans,[1,2,3])
                    plotchans = [5,4,3];
                elseif isequal(plotchans,[5,4,3])
                    plotchans = [1,2,3];
                end
                replot=1;


            case '0'  % back to default window - all data

                viewwind = bigwind;

            case 'o'  % zoom out 20%

                viewwind = mean(viewwind) + diff(viewwind)*[-0.6 0.6];

            case 28 % left arrow
                viewwind = viewwind - 0.1*diff(viewwind); % move 10% back

            case 29  % right arrow
                viewwind = viewwind + 0.1*diff(viewwind); % move 10% forward

            case 30 % up arrow
                ampf = ampf*1.4; % increase amplitudes
                replot=1;

            case 31  % down arrow
                ampf = ampf/1.4; % decrease amplitudes
                replot=1;

            case 27  % ESC key
                break

        end
        % EXIT if ENTER selected
        if isempty(b), cont = 1; end
    
    end % continue to process 


%% ------------------------------------------------------------------------
%% Subfunctions



% time info
function [ST,starttime] = parse_timing(event)
    ST.year = event.nzyear;
    ST.jday = event.nzjday;
    [ST.month,ST.day] = calday(event.nzjday);
    ST.hour = event.nzhour;
    ST.min = event.nzmin;
    ST.sec = event.nzsec;
    ST.msec = event.nzmsec;
    ST.datestr = [num2str(ST.year),'/',num2str(ST.jday),' ',num2str(ST.hour,'%02d'),':',num2str(ST.min,'%02d'),':',num2str(ST.sec,'%02d')];
    starttime = serialmsecprec(datenum([ST.year,ST.month,ST.day,ST.hour,ST.min,ST.sec+ST.msec/1e3]));
end

function [ST] = timeinfo(tserial)
    v = datevec(tserial);
    ST.year = v(1);
    ST.jday = doy(v(1:3));
    ST.month = v(2);
    ST.day = v(3);
    ST.hour = v(4);
    ST.min = v(5);
    ST.sec = round(v(6),3);
    ST.datestr = datestr(tserial,'yyyy/mm/dd HH:MM:SS.FFF');
end


function tout = msecprec(tin)
    % ensure times in seconds are to millisecond precision
    tout = round(tin,3);
end

function tout = serialmsecprec(tin)
    % ensure times in days (e.g. serial times) are to millisecond precision
    spd = 24*60*60; % seconds per day for serial times
    tout = round(tin*spd,3)/spd;
end

