function [] = DRS4stats()
clc; close all
%filename='scintillator.air.coupling.5gsps.xml';
%filename='standardTermination.xml';
%filename = 'scintillator.airCoupling.5gsps.50ohmTerm.ch2.xml';
[filename, pathname] = uigetfile('/Users/glennjocher/downloads/DRS4/*.xml','Select File:');
%filename='noise.4channels.xml'; pathname='/Users/glennjocher/Downloads/';

C=loadDRS4xml(filename, pathname);

% for i=1:100
%     dx(i,:) = getTimebaseConstants(C,i);
% end
% ''


C=C(1:2:end,:);
ne=size(C,1);
timem=(0:1023)'*.2; %monotonic 5 GSPS time vector
headers=cell2mat(C(:,1)')';  hours=etime(datevec(headers(end,2)),datevec(headers(1,2)))/3600;
T=cell2mat(C(:,2)');  GSPS=1024./(T(end,:)-T(1,:));
fprintf('%g events, %.2f GSPS, %.1f hours\n',ne,mean(GSPS),hours);


%x=0:5; y=poisspdf(x,2.88/5); fig; stem(x,y,'filled','o','linewidth',3); xyzlabel('SiPM Dark Counts per DRS4 5GSPS Timebase','Probability')

%RUN STATS
V=-cell2mat(C(:,3)'); %3-raw, 4-resampled
V=removeSpikes(V,2);
V=removeSpikes(V,[4.5 1]);
%porch=mean(V(100:200,:),'omitnan');  V=V-porch;
V=fcnsmooth(V,30);
porch = nanporch(V,3);
V=V-porch;
V([1:50 974:1024],:)=0;

fig(3,3); for i=1:9; sca; plot(V(20:end-50,i+20)); end


[j, a, w, t, integral, x] = fcnpulsewidth(single(V),8,[4 1100],[5 80],1100,'amplitude');% k=t==0;  a(k)=nan;  w(k)=nan; integral(k)=nan;
mu = [j a w t integral];


%PLOTTING
ha=fig(1,5);
tstr=filename;
str={'Candidates','Amplitude (mV)','Width (bins)','Time (bins)','Integral (bins)'};
j=true(size(mu,1),1); for i=1:5; [~,jmu]=fcnsigmarejection(mu(:,i),6,1);  j=j & jmu;  hv{i,1}=linspace(min(0,min(mu(j,i))),max(mu(j,i)),120); end
for i=1:5; sca(ha(i)); histogramN(mu(:,i),hv{i,1}); xlabel(str{i}); axis tight; end
title(ha(3),tstr)

fig; k=a>4; plot(timem,V(:,k),'.-'); xyzlabel('t (ns)','amplitude (mV)','',''); fcnmarkersize(1); axis tight

fig(1,3,1.5);
mu=nanmean(V'); s=nanstd(V');  Vn=V./a'; %#ok<UDIM>
sca; plot(timem,V(:,:),'.-'); xyzlabel('t (ns)','amplitude (mV)','','');
sca; plot(timem,Vn(:,:),'.-'); xyzlabel('t (ns)','amplitude (Normalized)','')
sca; errorarea(timem',mu,s,'b'); fcntight; xyzlabel('t (ns)','amplitude (mV)'); legend show
end


function histogramN(X,bins)
sigma=(max(bins)-min(bins))*.01;

mu=linspace(min(bins),max(bins),500);
Y = nansum(pdf('norm',mu',X',sigma),2);

plot(mu,Y);
end


function dxmu=getTimebaseConstants(C,ei)
a=cell2mat(C(:,1)'); ch=a(6,:);  TC=a(4,:);  nc=numel(unique(ch));   DT=cell(1,nc);  T=DT; dxmu=zeros(4,1);
%DTT=load('timebase.2675.mat');

%fig(1,1,'15cm');
for i=1:nc
    Ca=C(i==ch,:);  n=size(Ca,1); t=cell(1,n); d=t; du=t;  tc=TC(i==ch);
    for j=1:n
        t{j}=[Ca{j,2}; nan];
        d{j}=diff(t{j});
        du{j}=circshift(d{j},tc(j));
    end
    DT{i}=round(nanmean(cell2mat(du),2),5);
    %DT{i}=DTT.DT(:,i);
    
    dthat=circshift(DT{i},1024-tc(ei));  x=[t{ei}(1:end-1),  cumsum([0; d{ei}(1:end-1)]),  cumsum([0; dthat(1:end-1)])];
    T{i}=x(:,3);
    
    %CH OFFSET
    if i>1
        j=mod(700-tc(ei),1024)+1;
        T{i} = T{i} - (round(T{i}(j)-T{1}(j),3));
    end
    
    %plot(x(:,1)-T{i},'Display',sprintf('CH%g %.4f',i,mean(x(:,1)-T{i})));
    dxmu(i)=mean(x(:,1)-T{i});
end
that = triggerCell2timebase(ei*ones(4,1),2675*ones(4,1),(0:3)',tc(ei)*ones(4,1));
minmax3(cell2mat(T) - that')

dxmu
%fcntight('xy joint'); fcntight('jointeveny'); xyzlabel('T (bin)','error (ns)'); legend show
%DT=single(cell2mat(DT));  save timebase.2674.mat
end