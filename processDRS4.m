function [] = processDRS4()
clc; f=''; tic
%p='/Users/glennjocher/Downloads/DATA/DRS4/Serge/11.25/';
p='/Users/glennjocher/Downloads/';
%f=uigetfile(p,'*.*');

[f, p] = uigetfile([p '/*.*'],'Select text file:',fcnlastfileloaded,'MultiSelect','on'); 
if ~iscell(f); f={f}; end

fileTypeFlag = 'IRS3D'; %'IRS3D' or 'DRS4'

for i=1:numel(f)
    close all;
    switch fileTypeFlag
        case 'DRS4'
            A=fcnloadtextfile([p f{i}],inf,0);  runDRS4(A)
        case 'IRS3D'
            A=fcnloadtextfile([p f{i}],1E7,3);  runIRSD(A)
    end
end
toc
end


function V=construct1IRSDevent(W,Wa,pid,Vx)
i=find(pid==22 | pid==83); %event number
pid=pid(i);  minWa=min(Wa);
Wa=Wa(i);
W=W(i);

nw=max(Wa)-min(Wa)+1;  nj=numel(pid);  nw=4;
V = zeros(nw*64,128,'int16');
k = sub2ind([64 nw 128],ones(nj,1),Wa-minWa+1,pid);  %DOUBLE BOOKEND
V(k+(0:63)) = Vx(:,i)';
V=single(V); %xa(xa==0)=nan;
V=V(:,[ 22 83 23 82]);

% %TIMEBASE AND MONOTONIC RESAMPLING
% A=webread('http://www.phys.hawaii.edu/~slavli/mtc-nist/dT/dTvalues_SCROD_MS215_carrier_2_ASIC_3_ch_1.txt');
% [G.ve,startrow,~,D] = fcnunique(double(G.x(:,2)));  endrow=startrow+D-1;  G.ei = [startrow, endrow];  %event start/end indices
% 
% nw=32; %number of windows
% RW=G.x(:,7); %reference window
% W=G.x(:,8); %window
% Wa=W-RW+1;  Wa(W<RW)=Wa(W<RW)+nw; %window adjusted
% t=(0:63)+(Wa*64); %(samples) absolute time
% t256=mod(t,256); %(samples) 256-vector time
% t256(:,1)

% %RESAMPLE DATA (CODE FROM DRS4)
% if nargin==3
%     timem=single((0:255)'/2.71); %monotonic 2.71 GSPS time vector (368 ps samples)
%     for i=1:ne
%        X(:,i)=interp1cfloat(t256(i,:),X(:,i),timem);
%     end
% end

end


function []=runIRSD(G)
f=G.filename;  p=G.pathname;
if ~strcmp(f(end-3:end),'.mat') && ~exist([p f '.mat'],'file'); A=G; save([p f '.mat'],'A','-v6'); delete([p f]); clear A; end
darkCounts = any(regexpi(f,'dark'));

G.E=G.x(:,2);  G.unixTime=G.double/1E6;  G.events=fcnunique(G.E);  %A.x=int16(A.x);
G.ve = G.events;  ne=numel(G.ve);
G.ne=ne;
G.promptFlag=true;

pfile = 'pedestals.exp_0099_run_0180.glenn.mat';
load(pfile); %load pedestals
G.pedestals=muvb;
G.pedestalOutliers=false(size(pedestalOutliers));
G.SRCCHi = MTCpixelID2SRCCH(1:1536);
load MTCoffsets.mat; G.offsets=map; %pixel offsets (1-1536);
H = design(fdesign.lowpass('Fp,Fst',0.01,.2), 'equiripple');  % Design an FIR equiripple filter
H.Arithmetic = class(single(1));  G.filter = H;
assignin('base','G',G)

dt = max(G.unixTime)-minnonzero(G.unixTime); %(s)
fprintf('\n%g events over %.1f s (%.1f Hz)\n%g missing unix times',G.ne,dt,G.ne/(dt),mean(G.unixTime==0));
[~,~,~,n]=fcnunique(G.x(:,2));  fprintf(', %g missing windows\n',sum(512-n));

%PRUNE
pid = NTCIRSDmap(G.x(:,3:6),'pixel');
i = ~any(pid==[106 98 108 99],2);  %ArrayX pruning list
%i=true(size(mtcpid));
G.x = G.x(i,:);  pid=pid(i);  mtcpid = MTCSRCCH2pixelID(G.x(:,3:6));  G.E=G.E(i);

nw=32;  RW=G.x(:,7);  W=G.x(:,8);  Wa=W-RW+1;  Wa(W<RW)=Wa(W<RW)+nw;
Vx=int16(G.x(:,10:73)');

%SUBTRACT PEDESTALS
Vx = subtractpedestals(G.pedestals,G.pedestalOutliers,mtcpid,W,Vx);

%CONSTRUCT WAVEFORMS
B=cell(G.ne,1);
a=[1; find(diff(G.E)); numel(G.E)];
for i=1:G.ne
    b=a(i):a(i+1);
    B{i}=construct1IRSDevent(W(b),Wa(b),pid(b),Vx(:,b));
end
i=~cellfun(@isempty,B); G.ne=sum(i);
B=B(i);

x=[B{:}];  zv=zeros(G.ne*4,1);
event = fcncol([1 1 1 1]'*(1:G.ne));
board = zv;
ch = repmat([1 0 2 3],[1 G.ne])';

try
    fid = fopen([p f(1:strfind(f,'.glenn')) 'event_steps']);
    M.x = textscan(fid,'%f',inf,'Headerlines',0);  M.x=M.x{:};  fclose(fid); %a little faster for small files
    spacing = 5; %(mm) laser steps
    motor=interp1([min(G.ve); M.x]+(0:numel(M.x))'*1E-3,(0:numel(M.x))*400*spacing/10,fcncol((G.ve*[1 1 1 1])'),'previous');
catch
    motor = zv;
end
temp=zv; %temperature
n=numel(event)/4; %numel(unique(event));

x=removeSpikes(x,2);
%x=removeSpikes(x,[4.5 1]);
%[~,x]=triggerCell2timebase(event,board,ch,tc,x);
%[~,x]=timebaseIRSD(G,pid,x);

%porch=mean(x(31:200,:),'omitnan');  
porch=nanporch(x,10);
x=x-porch;  
%x([1:30 924:1024],:)=0;
%x=fcnsmooth(x,19);
%H=filterV(x(:,1));  x(isnan(x))=0;  x=filter(H,x);
if darkCounts; x=-x; end %FOR DARK COUNTS

[~, aa, wa, ta, inta] = fcnpulsewidth(x,.5,[1 1000],[15 160],1000,'fraction'); %[10-60] width

%j = find(aa>5 & circshift(aa,-4)>5);
%fig; plot(x(:,j(3)+[0 4]))
%dt=fcnsigmarejection(ta(j)-ta(j+4),2.8,6); fig; fhistogram(dt,200)

z=zeros(n,4); a=z; w=z; t=z; int=z; m=z(:,1); I=cell(n,4);  cch=z;  dt=t(:,1);
for i=1:numel(event)
    r=event(i);
    c=ch(i)+1;
    cch(r,c)=ch(i);
    a(r,c)=aa(i)/1000*.333;
    w(r,c)=wa(i); %width
    t(r,c)=ta(i)*.360;

    int(r,c)=inta(i)/1000;
    m(r)=-motor(i)/400; %windings to mm (usually centers about m/400-300)
    if i<1000000
        I{r,c}=x(:,i)'/1000*.333; 
    end
    %if c==2 && w(r,c)>0
    %    dt(r) = 1024 - fcnSplinePulseMaxima(xcorr(x(:,i-1),x(:,i))); %XCORR DT
    %end
end
dt = t(:,2)-t(:,1);  m = m + 9.75; %mean3(m); %amplitudeCrossingPosition(a(:,1)./a(:,2),m);

i = a>.001 & w>10 & w<160 & a<.5 & abs(dt)<5 & t(:,1)~=0 & t(:,2)~=0 & m~=max(m) & m~=min(m);
i = all(i(:,1:2),2) & max(abs(int(:,2)./int(:,1)),abs(int(:,1)./int(:,2)))<20;
[~,j] = fcnsigmarejection(dt(i,:),3,3); i(i)=j; I=I(i,:); nj=sum(i); %REJECT OUTLIER DTs


% ha=fig(2,2,1);
% ai=(1:4:400)+0; bi=(2:4:400)+0; c1=fcndefaultcolors(1); c2=fcndefaultcolors(2);
% porch=reshape(porch,4,n)'; event=reshape(event,4,n)'; temp=reshape(temp,4,n)';
% sca; plot(A.x(ai,8:end)','Color',c1); plot(A.x(bi,8:end)','Color',c2); xyzlabel('T (bin)','mV','',[repmat(' ',1,60) str_(f)]); legend('BEFORE'); legend boxoff
% sca; plot(x(:,ai),'Color',c1); plot(x(:,bi),'Color',c2); xyzlabel('T (bin)','mV'); fcnlinewidth(1.5); legend('AFTER'); legend boxoff
% sca; for j=1:4; plot(event(:,j),porch(:,j),'.','Display',sprintf('CH%g %.1f mV',j,mean(porch(:,j)))); end; legend show; xyzlabel('event','porch (mV)'); 
% sca; plot(event(:,1),temp(:,1),'.'); xyzlabel('event','temp (^oC)')
% linkaxes(ha([1 2])); ha(1).XLim=[450 800]; ha(2).YLim=[-15 80]; fcnfontsize(14);

sprintf('%.2g candidate eff and %.2g mV mean amplitude\n',mean(i),mean3(a(i,1:2))*1E3)
if sum(i)<100 && darkCounts; return; end %not enough candidate waveforms
xa=m(i)*10;  T=xa*[1 0 0 0]; %T=[X T PE E]


%DARK COUNTS
%fig; for j=1; histogram(a(a(:,j)>0.000,j),linspace(0,.07,200)); end; fcntight; set(gca,'XLim',[0 .07]); set(gca,'YScale','log'); xyzlabel('Amplitude (V)','','',f)

if any(i)
    n=5;
    ha=fig(3,n,.8,.8,[1.5, 2, 2.2, 0.6, 1.6, 1.8]); for k=[1:n, round(nj/2)+(1:n), (nj-n):(nj-1)]; sca; for ci=1:2; xi=I{k,ci}; plot(xi); end; title(k); end
    fcntight('xyjoint'); fcnlinewidth(1.5); for j=1:(3*n); ha(j).CameraViewAngle=5.3; end
    ylabel(ha(1),'Start','fontsize',26); ylabel(ha(n+1),'Middle','fontsize',26); ylabel(ha(2*n+1),'End','fontsize',26); title(ha(1),sprintf('%s\n1',str_(f)))
    hb=ha([2:n, (2:n)+n, (2:n)+2*n]); yticklabels(hb,[]); 
    xticklabels(ha(1:2*n),[]); for j=1:(3*n); ha(j).CameraViewAngle=5.5; end
end

%TEMPERATURE DEPENDENCE PLOT
%fig; hist2(temp(:,1),porch(:,1),100,'profile'); xyzlabel('temp (C)','porch (mV)','',str_(f)); fcnlinewidth(2)

%NEURAL NET
%X=applyNN(I(i,:),T); %clear I


[um]=fcnunique(m); ha=fig(3,1,1,1.5);  ne=0;  xa=1:numel(i);
load timebase
for j=1:numel(um)
    k = i & m==um(j);
    dtj = t(k,2)-t(k,1);
    dtc = dtj-timebase(t(k,1));
    plot(xa(k),dtj,'.');
    %[std(dtj) std(dtc)]
end; xyzlabel('event','dt (ns)','',str_(f)); 
F=fit(xa(i)',dt(i),'poly1'); plot(xa(i),F(xa(i)),'-'); 
sca; h=hist2(xa(i),dt(i)-F(xa(i)),40,'profile'); xyzlabel('event','residuals (ns)')
sca; plot(h.XData,h.YPositiveDelta,'.-');  xyzlabel('event','1\sigma (ns)'); fcntight(gca,'y0'); fcnlinewidth(2);  fcnmarkersize(8); for j=1:3; ha(j).XDir='reverse'; end


clear h; tsigma=''; xsigma=''; slope=0;
if any(diff(m)); staticSource=false; else; staticSource=true; end
if staticSource
    ha=fig(2,4,1); nb=75;
    sca; b=fcnsigmarejection(t(i,1:2)); xa=linspace(minnonzero(b),max3(b),nb); h(1)=fhistogram(t(i,1),xa); h(2)=fhistogram(t(i,2),xa); xyzlabel('(ns)','','','time'); legend(sprintf('CH1 LEFT%s',h(1).DisplayName),sprintf('CH2 RIGHT%s',h(2).DisplayName),'Location','best');
    sca; b=fcnsigmarejection(a(i,1:2),4,2); xa=linspace(0,max3(b),nb); fhistogram(a(i,1),xa); fhistogram(a(i,2),xa); xyzlabel('(V)','','','amplitude')
    sca; xa=linspace(10,30,nb); fhistogram(w(i,1),xa); fhistogram(w(i,2),xa); xyzlabel('(bins)','','','width')
    sca; b=fcnsigmarejection(int(i,1:2),4,2); xa=linspace(0,max3(b),nb); fhistogram(int(i,1),xa); fhistogram(int(i,2),xa); xyzlabel('(V)','','','integral')
    sca; h=fhistogram((t(i,2)-t(i,1)),nb); xyzlabel('dt (ns)','','',''); h.FaceColor=[.7 .7 .7];
    sca; h=fhistogram(a(i,1)./a(i,2),nb); xyzlabel('ratio','','',''); h.FaceColor=[.7 .7 .7];
    sca; h=fhistogram(w(i,1)./w(i,2),nb); xyzlabel('ratio','','',''); h.FaceColor=[.7 .7 .7];
    sca; h=fhistogram(int(i,1)./(int(i,2)),nb); xyzlabel('ratio','','',''); h.FaceColor=[.7 .7 .7];
    for j=1:8; legend(ha(j),'show'); legend(ha(j),'boxoff'); ha(j).YTick=[]; ha(j).YColor=[1 1 1]; end
    text(ha(2),0,ha(2).YLim(2)*1.14,str_(f));  fcnlinewidth(1.5);  fcnfontsize(17); fcntight('x')
    %export_fig('-a1','-q90','-r250',[G.filename '.png'])
    
    fig(1,1,1.5); for j=1:4; X=x(:,ch==j-1);  fcnpulsetemplate(X(:,i),t(i,1)); end; fcntight; xyzlabel('T (sample)','Amplitude (mV)','',str_(f))
else
    ha=fig(2,4,1);  type='profile';  nb=1000; %round(diff(minmax3(xa))/10);  %nb=100;
    YL=[2 12; 0 .03; 10 30; 0 .8; -2.5 2.5;0 3; .8 1.2;  0 3];
    sca; for j=1:2; h(j)=hist2(xa,t(i,j),nb,type); end; xyzlabel('','(ns)','','time'); legend(sprintf('CH1 LEFT%s',h(1).DisplayName),sprintf('CH2 RIGHT%s',h(2).DisplayName),'Location','best');
    sca; for j=1:2; b=a(i,j); h=hist2(xa,b,nb,type); exp1fit(h,xa,b); end;  h=hist2(xa,a(i,1)+a(i,2),nb,type); h.DisplayName='sum'; xyzlabel('','(V)','','amplitude'); legend show; legend boxoff; legend('Location','northwest'); %delete(h)
    sca; hist2(xa,w(i,1),nb,type); hist2(xa,w(i,2),nb,type); xyzlabel('','(bins)','','width')    
    sca; for j=1:2; b=int(i,j); h=hist2(xa,b,nb,type); exp1fit(h,xa,b); end;  h=hist2(xa,int(i,1)+int(i,2),nb,type); h.DisplayName='sum'; xyzlabel('','(V)','','integral'); legend show; legend boxoff; legend('Location','northwest'); %delete(h)
    sca; b=dt(i,:); h=hist2(xa,b,nb,type); sx=''; sy=''; try h.Color=[.7 .7 .7]; sy=h.DisplayName; h=hist2(b,xa,nb,'profile'); sx=h.DisplayName; delete(h); end;  F=fit(xa,b,'poly1'); hl=plot(xa,F(xa),'Display',sprintf('%.1f mm/ns',1./abs(F.p1))); xyzlabel(sprintf('X (mm) %s',sx),sprintf('(ns) %s',sy),'','dt'); legend(hl,'Location','northwest'); xsigma=sx; tsigma=sy; slope=1./abs(F.p1);
    sca; b=a(i,1)./a(i,2); h=hist2(xa,b,nb,type); F=fit(xa,b,'exp1');  h.DisplayName=sprintf('\\lambda = %.0f mm',1/abs(F.b)); legend show; legend('Location','northwest'); sx=''; try h.Color=[.7 .7 .7]; h=hist2(b,xa,nb,type); sx=h.DisplayName; delete(h); end; xyzlabel(sprintf('X (mm) %s',sx),'','','ratio');
    sca; h=hist2(xa,w(i,1)./w(i,2),nb,type); xyzlabel('X (mm)','','','ratio'); try h.Color=[.7 .7 .7]; end
    sca; b=int(i,1)./int(i,2); h=hist2(xa,b,nb,type); F=fit(xa,b,'exp1');  h.DisplayName=sprintf('\\lambda = %.0f mm',1/abs(F.b)); legend show; legend('Location','northwest'); sx=''; try h.Color=[.7 .7 .7]; h=hist2(b,xa,nb,type); sx=h.DisplayName; delete(h); end; xyzlabel(sprintf('X (mm) %s',sx),'','','ratio');
    for j=1:8; ha(j).XLim=[-1 1]*100; end;
    for j=[5 6 7 8]; ha(j).YLim=YL(j,:); end;
    for j=[2 4 6 8]; ha(j).YLim(1)=0; end;
    text(ha(2),0,ha(2).YLim(2)*1.14,str_(f));  fcnlinewidth(1.5);  fcnfontsize(16);
    %export_fig('-a1','-q90','-r250',[A.filename '.png'])
    
    
    %[~,i]=sort(h.YData); F=griddedInterpolant(h.YData(i),h.XData(i),'spline','linear'); save F.mat F
    
    %h=get(gca,'Children'); legend(sprintf('CH1 LEFT%s',h(1).DisplayName),sprintf('CH2 RIGHT%s',h(2).DisplayName),'Location','best'); legend boxoff
    %h=gca; h.YLim=[-1 1]; legend('CH1','CH2','CH3','CH4');
    
    
%     fig; A=a(i,1); B=a(i,2); xA=xa+500; xB=-xa+500-65; hist2(xA/10,A/.015,nb,type); hist2(xB/10,B*1.076/.015,nb,type); xyzlabel('Distance (cm)','Amplitude (normalized at 30 cm)')
%     sx=10:10:300; sy=1/741.1*[ 826.0, 792.5, 741.1, 704.9, 681.8, 554.5, 546.3, 488.3, 483.5, 510.0, 446.5, 456.0, 432.6, 418.5, 414.3, 362.8, 387.0, 359.5, 364.1, 370.9, 322.2, 330.1, 307.8, 305.0, 304.5, 266.1, 283.8, 280.3, 280.4, 289.0]; plot(sx,sy,'.-'); fcnlinewidth(1.5)
%     legend('UH 0-80 cm','UH 20-100 cm','Saint Gobain 10-300 cm')
    
%    plotCurves(m(:,1),x,i);
end
''
% %PRINT OUTPUT
% try %#ok<TRYNC>
%     try load NTCruns.mat; catch; X={}; save NTCruns.mat X; end
%     Xa={A.filename,mode(board),abs(diff(minmax3(temp))),abs(diff(minmax3(m)))/10,mean(i),eval(xsigma(isnumericstr(xsigma))),eval(tsigma(isnumericstr(tsigma)))*1E3,slope,mean3(a(i,1:2))*1E3,mean3(w(i,1:2)),mean3(int(i,1:2))};
%     X(size(X,1)+1,:)=Xa;  save NTCruns.mat X;
% end
end


function []=runDRS4(A)
f=A.filename;  p=A.pathname;
if ~strcmp(f(end-3:end),'.mat') && ~exist([p f '.mat'],'file'); save([p f '.mat'],'A','-v6'); delete([p f]); end
darkCounts = any(regexpi(f,'dark'));

h=A.x(:,1:7); event=h(:,1); board=h(:,2); ch=h(:,3); tc=h(:,4);  motor=h(:,5);  temp=h(:,7);
x=A.x(:,8:end)';  A.x=A.x(1:400,:); n=numel(event)/4; %numel(unique(event));
CH34pedestal=(x(:,3:4:end)+x(:,4:4:end)*0)/1; for i=1:4; x(:,i:4:end)=x(:,i:4:end)-CH34pedestal; end;  clear CH34pedestal

x=removeSpikes(x,2);
x=removeSpikes(x,[4.5 1]);
[~,x]=triggerCell2timebase(event,board,ch,tc,x);
%porch=mean(x(31:200,:),'omitnan');  
porch=nanporch(x,3);
x=x-porch;  
x([1:30 924:1024],:)=0;
%x=fcnsmooth(x,19);
H=filterV(x(:,1));  x(isnan(x))=0;  x=filter(H,x);
if darkCounts; x=-x; end %FOR DARK COUNTS

[~, aa, wa, ta, inta] = fcnpulsewidth(x,.5,[1 1000],[15 160],500,'fraction'); %[10-60] width

%j = find(aa>5 & circshift(aa,-4)>5);
%fig; plot(x(:,j(3)+[0 4]))
%dt=fcnsigmarejection(ta(j)-ta(j+4),2.8,6); fig; fhistogram(dt,200)

z=zeros(n,4); a=z; w=z; t=z; int=z; m=z(:,1); I=cell(n,4);  cch=z;  dt=t(:,1);
for i=1:numel(event)
    r=event(i);
    c=ch(i)+1;
    cch(r,c)=ch(i);
    a(r,c)=aa(i)/1000;
    w(r,c)=wa(i); %width
    t(r,c)=ta(i)*.200;

    int(r,c)=inta(i)/1000;
    m(r)=-motor(i)/400; %windings to mm (usually centers about m/400-300)
    if i<900000; I{r,c}=x(:,i)'; end
    %if c==2 && w(r,c)>0
    %    dt(r) = 1024 - fcnSplinePulseMaxima(xcorr(x(:,i-1),x(:,i))); %XCORR DT
    %end
end
dt = t(:,2)-t(:,1);  m = m - dtCrossingPosition(dt,m);

i = a>.001 & w>=15 & w<160 & a<.148 & abs(dt)<8;
i = all(i(:,1:2),2) & max(abs(int(:,2)./int(:,1)),abs(int(:,1)./int(:,2)))<20;
[~,j] = fcnsigmarejection(dt(i,:),3,3); i(i)=j; I=I(i,:); nj=sum(j); %REJECT OUTLIER DTs


% ha=fig(2,2,1);
% ai=(1:4:400)+0; bi=(2:4:400)+0; c1=fcndefaultcolors(1); c2=fcndefaultcolors(2);
% porch=reshape(porch,4,n)'; event=reshape(event,4,n)'; temp=reshape(temp,4,n)';
% sca; plot(A.x(ai,8:end)','Color',c1); plot(A.x(bi,8:end)','Color',c2); xyzlabel('T (bin)','mV','',[repmat(' ',1,60) str_(f)]); legend('BEFORE'); legend boxoff
% sca; plot(x(:,ai),'Color',c1); plot(x(:,bi),'Color',c2); xyzlabel('T (bin)','mV'); fcnlinewidth(1.5); legend('AFTER'); legend boxoff
% sca; for j=1:4; plot(event(:,j),porch(:,j),'.','Display',sprintf('CH%g %.1f mV',j,mean(porch(:,j)))); end; legend show; xyzlabel('event','porch (mV)'); 
% sca; plot(event(:,1),temp(:,1),'.'); xyzlabel('event','temp (^oC)')
% linkaxes(ha([1 2])); ha(1).XLim=[450 800]; ha(2).YLim=[-15 80]; fcnfontsize(14);

sprintf('%.2g candidate eff and %.2g mV mean amplitude\n',mean(i),mean3(a(i,1:2))*1E3)
if sum(i)<100 && darkCounts; return; end %not enough candidate waveforms
xa=m(i);  T=xa*[1 0 0 0]; %T=[X T PE E]


%DARK COUNTS
%fig; for j=1; histogram(a(a(:,j)>0.000,j),linspace(0,.07,200)); end; fcntight; set(gca,'XLim',[0 .07]); set(gca,'YScale','log'); xyzlabel('Amplitude (V)','','',f)

if any(j)
    n=5;
    ha=fig(3,n,.8,.8,[1.5, 2, 2.2, 0.6, 1.6, 1.8]); for k=[1:n, round(nj/2)+(1:n), (nj-n):(nj-1)]; sca; for ci=1:2; xi=I{k,ci}; plot(xi); end; title(k); end
    fcntight('xyjoint'); fcnlinewidth(1.5); for j=1:(3*n); ha(j).CameraViewAngle=5.3; end
    ylabel(ha(1),'Start','fontsize',26); ylabel(ha(n+1),'Middle','fontsize',26); ylabel(ha(2*n+1),'End','fontsize',26); title(ha(1),sprintf('%s\n1',str_(f)))
    hb=ha([2:n, (2:n)+n, (2:n)+2*n]); yticklabels(hb,[]); 
    xticklabels(ha(1:2*n),[]); for j=1:(3*n); ha(j).CameraViewAngle=5.5; end
end
%TEMPERATURE DEPENDENCE PLOT
%fig; hist2(temp(:,1),porch(:,1),100,'profile'); xyzlabel('temp (C)','porch (mV)','',str_(f)); fcnlinewidth(2)

%NEURAL NET
%X=applyNN(I(i,:),T); %clear I


% [um]=fcnunique(m); ha=fig(3,1,1,1.5);  ne=0;  xa=1:numel(i);
% load timebase
% for j=1:numel(um)
%     k = i & m==um(j);
%     dtj = t(k,2)-t(k,1);
%     dtc = dtj-timebase(t(k,1));
%     plot(xa(k),dtj,'.');
%     %[std(dtj) std(dtc)]
% end; xyzlabel('event','dt (ns)','',str_(f)); 
% F=fit(xa(i)',dt(i),'poly1'); plot(xa(i),F(xa(i)),'-'); 
% sca; h=hist2(xa(i),dt(i)-F(xa(i)),40,'profile'); xyzlabel('event','residuals (ns)')
% sca; plot(h.XData,h.YPositiveDelta,'.-');  xyzlabel('event','1\sigma (ns)'); fcntight(gca,'y0'); fcnlinewidth(2);  fcnmarkersize(8); for j=1:3; ha(j).XDir='reverse'; end


clear h; tsigma=''; xsigma=''; slope=0;
if any(diff(m)); staticSource=false; else; staticSource=true; end
if staticSource
    ha=fig(2,4,1); nb=75;
    sca; b=fcnsigmarejection(t(i,1:2)); xa=linspace(minnonzero(b),max3(b),nb); h(1)=fhistogram(t(i,1),xa); h(2)=fhistogram(t(i,2),xa); xyzlabel('(ns)','','','time'); legend(sprintf('CH1 LEFT%s',h(1).DisplayName),sprintf('CH2 RIGHT%s',h(2).DisplayName),'Location','best');
    sca; b=fcnsigmarejection(a(i,1:2),4,2); xa=linspace(0,max3(b),nb); fhistogram(a(i,1),xa); fhistogram(a(i,2),xa); xyzlabel('(V)','','','amplitude')
    sca; xa=linspace(20,60,nb); fhistogram(w(i,1),xa); fhistogram(w(i,2),xa); xyzlabel('(bins)','','','width')
    sca; b=fcnsigmarejection(int(i,1:2),4,2); xa=linspace(0,max3(b),nb); fhistogram(int(i,1),xa); fhistogram(int(i,2),xa); xyzlabel('(V)','','','integral')
    sca; h=fhistogram((t(i,2)-t(i,1)),nb); xyzlabel('dt (ns)','','',''); h.FaceColor=[.7 .7 .7];
    sca; h=fhistogram(a(i,1)./a(i,2),nb); xyzlabel('ratio','','',''); h.FaceColor=[.7 .7 .7];
    sca; h=fhistogram(w(i,1)./w(i,2),nb); xyzlabel('ratio','','',''); h.FaceColor=[.7 .7 .7];
    sca; h=fhistogram(int(i,1)./(int(i,2)),nb); xyzlabel('ratio','','',''); h.FaceColor=[.7 .7 .7];
    for j=1:8; legend(ha(j),'show'); legend(ha(j),'boxoff'); ha(j).YTick=[]; ha(j).YColor=[1 1 1]; end
    text(ha(2),0,ha(2).YLim(2)*1.14,f);  fcnlinewidth(1.5);  fcnfontsize(17); fcntight('x')
    export_fig('-a1','-q90','-r250',[A.filename '.png'])
else
    ha=fig(2,4,1);  type='profile';  nb=round(diff(minmax3(xa))/10);  %nb=100;
    YL=[2 12; 0 .03; 15 40; 0 .8; -10 10;0 5; .6 1.6;  0 5];
    sca; for j=1:2; h(j)=hist2(xa,t(i,j),nb,type); end; xyzlabel('','(ns)','','time'); legend(sprintf('CH1 LEFT%s',h(1).DisplayName),sprintf('CH2 RIGHT%s',h(2).DisplayName),'Location','best');
    sca; for j=1:2; b=a(i,j); h=hist2(xa,b,nb,type); exp1fit(h,xa,b); end;  h=hist2(xa,a(i,1)+a(i,2),nb,type); h.DisplayName='sum'; xyzlabel('','(V)','','amplitude'); legend show; legend boxoff; legend('Location','northwest'); delete(h)
    sca; hist2(xa,w(i,1),nb,type); hist2(xa,w(i,2),nb,type); xyzlabel('','(bins)','','width')    
    sca; for j=1:2; b=int(i,j); h=hist2(xa,b,nb,type); exp1fit(h,xa,b); end;  h=hist2(xa,a(i,1)+a(i,2),nb,type); h.DisplayName='sum'; xyzlabel('','(V)','','integral'); legend show; legend boxoff; legend('Location','northwest'); delete(h)
    sca; b=dt(i,:); h=hist2(xa,b,nb,type); sx=''; sy=''; try h.Color=[.7 .7 .7]; sy=h.DisplayName; h=hist2(b,xa,nb,'profile'); sx=h.DisplayName; delete(h); end;  F=fit(xa,b,'poly1'); hl=plot(xa,F(xa),'Display',sprintf('%.1f mm/ns',1./abs(F.p1))); xyzlabel(sprintf('X (mm) %s',sx),sprintf('(ns) %s',sy),'','dt'); legend(hl,'Location','northwest'); xsigma=sx; tsigma=sy; slope=1./abs(F.p1);
    sca; b=a(i,1)./a(i,2); h=hist2(xa,b,nb,type); F=fit(xa,b,'exp1');  h.DisplayName=sprintf('\\lambda = %.0f mm',1/abs(F.b)); legend show; legend('Location','northwest'); sx=''; try h.Color=[.7 .7 .7]; h=hist2(b,xa,nb,type); sx=h.DisplayName; delete(h); end; xyzlabel(sprintf('X (mm) %s',sx),'','','ratio');
    sca; h=hist2(xa,w(i,1)./w(i,2),nb,type); xyzlabel('X (mm)','','','ratio'); try h.Color=[.7 .7 .7]; end
    sca; b=int(i,1)./int(i,2); h=hist2(xa,b,nb,type); F=fit(xa,b,'exp1');  h.DisplayName=sprintf('\\lambda = %.0f mm',1/abs(F.b)); legend show; legend('Location','northwest'); sx=''; try h.Color=[.7 .7 .7]; h=hist2(b,xa,nb,type); sx=h.DisplayName; delete(h); end; xyzlabel(sprintf('X (mm) %s',sx),'','','ratio');
    for j=1:8; ha(j).XLim=[-1 1]*500; end;
    for j=[5 6 7 8]; ha(j).YLim=YL(j,:); end;
    for j=[2 4 6 8]; ha(j).YLim(1)=0; end;
    text(ha(2),0,ha(2).YLim(2)*1.14,f);  fcnlinewidth(1.5);  fcnfontsize(16);
    %export_fig('-a1','-q90','-r250',[A.filename '.png'])
    
    %[~,i]=sort(h.YData); F=griddedInterpolant(h.YData(i),h.XData(i),'spline','linear'); save F.mat F
    
    %fig(1,1,1.5); for j=1:4; X=x(:,ch==j-1);  fcnpulsetemplate(X(:,i),t(i,1)); end; fcntight; xyzlabel('T (sample)','Amplitude (mV)','',f)
    %h=get(gca,'Children'); legend(sprintf('CH1 LEFT%s',h(1).DisplayName),sprintf('CH2 RIGHT%s',h(2).DisplayName),'Location','best'); legend boxoff
    %h=gca; h.YLim=[-1 1]; legend('CH1','CH2','CH3','CH4');
    
    
%     fig; A=a(i,1); B=a(i,2); xA=xa+500; xB=-xa+500-65; hist2(xA/10,A/.015,nb,type); hist2(xB/10,B*1.076/.015,nb,type); xyzlabel('Distance (cm)','Amplitude (normalized at 30 cm)')
%     sx=10:10:300; sy=1/741.1*[ 826.0, 792.5, 741.1, 704.9, 681.8, 554.5, 546.3, 488.3, 483.5, 510.0, 446.5, 456.0, 432.6, 418.5, 414.3, 362.8, 387.0, 359.5, 364.1, 370.9, 322.2, 330.1, 307.8, 305.0, 304.5, 266.1, 283.8, 280.3, 280.4, 289.0]; plot(sx,sy,'.-'); fcnlinewidth(1.5)
%     legend('UH 0-80 cm','UH 20-100 cm','Saint Gobain 10-300 cm')
    
%    plotCurves(m(:,1),x,i);
end

%PRINT OUTPUT
try %#ok<TRYNC>
    try load NTCruns.mat; catch; X={}; save NTCruns.mat X; end
    Xa={A.filename,mode(board),abs(diff(minmax3(temp))),abs(diff(minmax3(m)))/10,mean(i),eval(xsigma(isnumericstr(xsigma))),eval(tsigma(isnumericstr(tsigma)))*1E3,slope,mean3(a(i,1:2))*1E3,mean3(w(i,1:2)),mean3(int(i,1:2))};
    X(size(X,1)+1,:)=Xa;  save NTCruns.mat X;
end
end

function x=dtCrossingPosition(dt,m)
rows=min(1000,numel(dt));
dt=abs(dt);  dt(dt==0)=inf;  dt(:,2)=1:numel(dt);
a=sortrows(dt,1);  x=mean(m(a(1:rows,2)));
end

function x=amplitudeCrossingPosition(ar,m)
rows=min(1000,numel(ar));
ar=abs(ar-1);  ar(ar==0)=inf;  ar(:,2)=1:numel(ar);
a=sortrows(ar,1);  x=mean(m(a(1:rows,2)));
end

function []=exp1fit(h,xa,b)
F=fit(xa,b,'exp1');  Fb=1./abs(F.b)/10;
h.DisplayName=sprintf('%.0f mV, %.0f cm',mean(b)*1E3,Fb);
end

function []=exp2fit(h,xa,b)
F=fit(xa,b,'exp2');  Fb=1./abs(F.b)/10;  Fd=1./abs(F.d)/10;
if Fd>Fb 
    h.DisplayName=sprintf('%.0f mV, %.0f-%.0f cm',mean(b)*1E3,Fb,Fd);
else
    h.DisplayName=sprintf('%.0f mV, %.0f-%.0f cm',mean(b)*1E3,Fd,Fb);
end
end



function fitGMdist
load x.mat
a=x(:,1)+x(:,2);
options=statset; options.MaxIter=500;
GMModel=fitgmdist(a,1,'Options',options);
xa=linspace(0,4,300)';
y=pdf(GMModel,xa);

fig; h=histogram(a,xa); 
gain = y(1:end-1)./h.Values(:);  gain=1./mean(gain(h.Values(:)>100));
plot(xa,y*gain)
end

function plotCurves(m,x,ci)
xb=min(m):10:max(m); nb=numel(xb);
x=permute(reshape(x,1024,4,size(x,2)/4),[1 3 2]); %x2=[event, 1-1024, ch]
ha=fig(2,1,'19cm');
for i=1:nb
    j=m==xb(i) & ci;
    for hi=1:2
        c=fcnlightencolor(fcndefaultcolors(hi),i/nb*.9); plot(ha(hi),mean(x(:,j,hi),2),'Color',c,'DisplayName',sprintf('%g mm',xb(i)));
    end
end
linkaxes(ha); fcnlinewidth(1.5); fcntight('xyjoint')
sca(ha(1)); legend(ha(1).Children([1 end])); xyzlabel('T (bin)','mV','','LEFT')
sca(ha(2)); legend(ha(2).Children([1 end])); xyzlabel('T (bin)','mV','','RIGHT')
end


function X=applyNN(I,T)
nw=1024;  inputs=cell2mat(I(:,1:2))/1000; vi{1}=1:nw; vi{2}=vi{1}+nw;  n=size(I,1);
Ia=zeros(n,128*2);  v128=1:128;  v129=v128+128;  t=zeros(n,2);  a=t; w=t; int=w; %F=lowPassFilter(100E6,200E6,5E9);
for i=1:2
    %I = inputs(:,vi{i}); I = filter(F,I')';  nb=2^12;  I=round(I*nb)/nb;  inputs(:,vi{i})=I;  %12-bit digitized
    [~,a(:,i),w(:,i),t(:,i),int(:,i)]=fcnpulsewidth(inputs(:,vi{i})',.5,[0 1000],[0 1000],1E6,'fraction');
end
%ci=abs(diff(t,1,2))<50;  t=t(ci,:); a=a(ci,:); w=w(ci,:); int=int(ci,:); inputs=inputs(ci,:); T=T(ci,:); Ia=Ia(ci,:); n=sum(ci);

tbias=zeros(n,3);
for i=1:n
    [~,j]=max(a(i,:));
    ia=round(min(t(i,j)));  ia=floorandceil(ia-20,0,nw-128-1);  Ia(i,:)=inputs(i,ia+[v128 v128+nw]); tbias(i,2)=ia*.2;
end
Ia(:,v128)=Ia(:,v128)./max(Ia(:,v128),[],2); %NORMALIZE
Ia(:,v129)=Ia(:,v129)./max(Ia(:,v129),[],2); %NORMALIZE

Ib=inputs(:,1:nw*2);
ha=fig(2,1,'19x19cm'); plot(Ib(30-(1:20),:)','.-'); sca; plot(Ia(30-(1:20),:)','.-'); fcntight; xyzlabel(ha,'T (5 GSPS sample)','V'); ha(2).YLabel.String='Normalized';
Ic=Ia; Ic=[a w t int];  Ic=[Ic Ic(:,1:2:end)./(Ic(:,2:2:end)+eps), Ic(:,1:2:end)-Ic(:,2:2:end)];

applyFlag=0;
if applyFlag
    %name='NN.TS FiberPointSource - 3L SENSL J60035 5V 128ch BCF10 100k.mat'; %MC-TRAINED
    name='NN.DRS4 29.12.2016 - 75k.mat'; %REAL-WORLD TRAINED
    load(name);
    pos=net{1}(Ic')';  Id=[NNwaveformStats(Ib(:,1:nw),2) NNwaveformStats(Ib(:,257:end),2) pos(:,1)];
else
    name=sprintf('NN.DRS4 %i.%i.%i - %.0fk.mat',day(now),month(now),year(now),n/1E3);
    [net{1}, perf(1), sigmas(1)] = Copy_of_NNtrain(Ic,T(:,1)-tbias(:,1)*0,[20 20],1);
    pos=net{1}(Ic')';  Id=[NNwaveformStats(Ib(:,1:nw),2) NNwaveformStats(Ib(:,257:end),2) pos(:,1)];   [net{2}, perf(2), sigmas(3)]  =  NNtrain(Id,T(:,4),[20 20],1);
    save(name,'net','perf','sigmas');
end
X(1:2) = conventionalFiberFit(inputs, T);
X{3} = [pos pos*0 net{2}(Id')'] + tbias;

%PLOT
hf=findobj(0,'Name','DRS4plot');  if isempty(hf);  [ha,hf]=fig(2,3,'20x30cm');  hf.Name='DRS4plot';  else;  figure(hf);  end
nb=74; clear S;  sx={'X (mm)','T (ns)','E (MeV)'};  st={'Position','Time','Energy'};  sd={'FT 2mV','CFD 50%','NN'}; ha(2).Title.String=name;
ci=T(:,1)>-500 & T(:,1)<500;   xa=T(ci,1);
for i=1:3 %1:numel(X)
    e = X{i} - T(:,[1 2 4]);
    x=linspace(-500,500,nb);
    for j=1:3; sca; if(i)<3 && j>2; continue; end; if(i)==1; cla; end;  y=e(ci,j);  h=hist2(xa,y,nb,'profile');  h.DisplayName=sd{i};   [~,~,S{i,j}] = movingMean(xa,y,nb,0);  xyzlabel('X_{motor} (mm)',sx{j},'',st{j}); end
    for j=1:3; sca; if(i)<3 && j>2; continue; end; plot(x,S{i,j}.s(x),'.-','Display',sd{i});  xyzlabel('X_{motor} (mm)',sx{j}); end
end
fcntight('xjoint'); fcnmarkersize(ha(1:3),.1); ha(1).YLim=[-500 500]; legend(ha(1),'show')
end