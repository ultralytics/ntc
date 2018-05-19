function C = loadDRS4xml(filename, pathname)
%pathname='/Users/glennjocher/Downloads/DRS4/autoTrigger/';
if nargin==0; filename='scintillator.air.coupling.5gsps.xml'; end %for debugging
fid = fopen([pathname filename]);

if ~exist('maxlines','var'); maxlines = 1E6; end %max number of lines to read

%SCAN FOR HEADER ROW ------------------------------------------------------
for j=1:3; fgetl(fid); end; %skip 4

i=0;
C=cell(10000,5);
timem=single((0:1023)'*.2); %monotonic 5 GSPS time vector
dateformat = 'yyyy/mm/dd HH:MM:SS.FFF';
while ~feof(fid)
    if i>1E4; break; end %max number of waveforms
    fgetl(fid);
    
    try
        c=regexp(fgetl(fid),'<Serial>(.+)</Serial>','tokens');  event=c{1};
    catch
        try
            c=regexp(fgetl(fid),'<Serial>(.+)</Serial>','tokens');  event=c{1};
        catch
            break
        end
    end
    c=regexp(fgetl(fid),'<Time>(.+)</Time>','tokens');  datetime=c{1}{1};      
    
    fgetl(fid);
    fgetl(fid);
    c=regexp(fgetl(fid),'<Board_(.+)>','tokens');  board=c{1};
    c=regexp(fgetl(fid),'<Trigger_Cell>(.+)</Trigger_Cell>','tokens');  triggercell=c{1};
    
    for CH=1:4
        c=regexp(fgetl(fid),'<Scaler\d>(.+)</Scaler','tokens');  if isempty(c); break; end; scalar=c{1};
        %scalar=0;
        c=regexp(fgetl(fid),'<CHN(.+)>','tokens');  channel=c{1};
        
        headers = str2double([event '0' board triggercell scalar channel]);
        
        B=textscan(fid,'%*s%*s%f%f%*s%*s',1024,'Delimiter','>,<,');
        time=B{1};  data=B{2};
        
        datam = interp1(time,data,timem); %monotonic
        i=i+1;  C{i,1}=headers';  C{i,2}=time; C{i,3}=data; C{i,4}=datam; C{i,5}=datetime;
        fgetl(fid);
        fgetl(fid);
    end
    fgetl(fid);
end
fclose(fid);
C=C(1:i,:);

datestr=C(:,5);  dates=datenum(datestr',dateformat);  for j=1:i; C{j,1}(2)=dates(j,:); end %optionally dates=now2unix(dates)

C=C(:,1:4);
