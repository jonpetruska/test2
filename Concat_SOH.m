function All_DeploySOH = Concat_SOH(SOH)
utc=[];
supvolt =[];
cur = [];
temp =[];
lat =[];
lon =[];
elev =[];
gnss =[];
sat =[];
timestat =[];
lock =[];
timeun =[];
timeerr =[];
timequal =[];
v1 =[];
v2 =[];
v3 =[];
ext1 =[];
ext2 =[];
ext3 =[];

badval = NaN;


Ny = length(SOH.yrs);
for iy=1:Ny
    ystr = num2str(SOH.yrs(iy),'%04.f');
    Nm = SOH.(['yyyy',ystr])(1).Nm;
    for im=1:Nm
        mstr = num2str(SOH.(['yyyy',ystr])(1).mos(im),'%02.f');
        intStr = SOH.(['yyyy',ystr])(1).(['mm',mstr]);        
        for id = 1:length(SOH.(['yyyy',ystr]).(['mm',mstr]).days)
            dstr = num2str(intStr.days(id),'%02.f');
            indfield = SOH.(['yyyy',ystr])(1).(['mm',mstr]).(['dd',dstr]);
            
            serstart = datenum([ystr mstr dstr '00' '00' '00'],'yyyymmddHHMMSS');
            
            try
            utc =    [utc       ; datenum(indfield.timeUTC(:))];
            sat =    [sat       ; indfield.Sat(:)];
            supvolt =[supvolt   ; indfield.Supply_Voltage(:)];
            cur =    [cur       ; indfield.Current(:)];
            temp =   [temp      ; indfield.Temp(:)];
            lat =    [lat       ; indfield.Lat(:)];
            lon =    [lon       ; indfield.Lon(:)];
            elev =   [elev      ; indfield.Elev(:)];
            gnss =   [gnss      ; indfield.GNSS(:)];
            timestat=[timestat  ; indfield.Time_Status(:)];
            lock =   [lock      ; indfield.Lock(:)];
            timeun = [timeun    ; indfield.Time_Uncertainty(:)];
            timeerr =[timeerr   ; indfield.Time_Error(:)];
            timequal=[timequal  ; indfield.Time_Quality(:)];
            v1 =     [v1        ; indfield.mass1(:)];
            v2 =     [v2        ; indfield.mass2(:)];
            v3 =     [v3        ; indfield.mass3(:)];
            ext1 =   [ext1      ; indfield.Ext1(:)];
            ext2 =   [ext2      ; indfield.Ext2(:)];
            ext3 =   [ext3      ; indfield.Ext3(:)];
            catch
                utc =    [utc       ; serstart];
                sat =    [sat       ; badval];
                supvolt =[supvolt   ; badval];
                cur =    [cur       ; badval];
                temp =   [temp      ; badval];
                lat =    [lat       ; badval];
                lon =    [lon       ; badval];
                elev =   [elev      ; badval];
                gnss =   [gnss      ; badval];
                timestat=[timestat  ; badval];
                lock =   [lock      ; badval];
                timeun = [timeun    ; badval];
                timeerr =[timeerr   ; badval];
                timequal=[timequal  ; badval];
                v1 =     [v1        ; badval];
                v2 =     [v2        ; badval];
                v3 =     [v3        ; badval];
                ext1 =   [ext1      ; badval];
                ext2 =   [ext2      ; badval];
                ext3 =   [ext3      ; badval];
            end
        end % loop on days
    end % loop on mos
end % loop on years
All_DeploySOH.Time = utc;
All_DeploySOH.Sat = sat;
All_DeploySOH.Supply_Voltage = supvolt;
All_DeploySOH.Current = cur;
All_DeploySOH.Temp = temp;
All_DeploySOH.Lat = lat;
All_DeploySOH.Lon = lon;
All_DeploySOH.Elev = elev;
All_DeploySOH.GNSS = gnss;
All_DeploySOH.Time_Status = timestat;
All_DeploySOH.Lock = lock;
All_DeploySOH.Time_Uncertainty = timeun;
All_DeploySOH.Time_Error = timeerr;
All_DeploySOH.Time_Quality = timequal;
All_DeploySOH.mass1 = v1;
All_DeploySOH.mass2 = v2;
All_DeploySOH.mass3 = v3;
All_DeploySOH.Ext1 = ext1;
All_DeploySOH.Ext2 = ext2;
All_DeploySOH.Ext3 = ext3;

end
