function [that, X] = triggerCell2timebase(ei,board,ch,tc,X)
ei=ei(:);  board=board(:);  ch=ch(:);  tc=tc(:);  if any(ch==0); ch=ch+1; end

ne=numel(ei);
bi=zeros(ne,1);  bi(board==2674)=1;  bi(board==2675)=2; 

load timebase.2674.mat;  dt=DT';
load timebase.2675.mat;  dt=[dt; DT'];  dtc=cell(8,1); for i=1:8; dtc{i}=dt(i,:); end

T=dtc((bi-1)*4+ch,:);
dthat = cellfun(@circshift,T,num2cell(1024-tc+1),'UniformOutput',false);
that=cat(1,dthat{:});  that(:,1)=0;  that = cumsum(that,2);


%CH OFFSET
s=size(that);
j=mod(700-tc,1024)+1;
[~,b,c]=unique([ei bi],'rows'); ri0=b(c);
offsets = that(sub2ind(s,(1:ne)',j)) - that(sub2ind(s,ri0,j));
that = that - offsets;

%RESAMPLE DATA
if nargin==5
    timem=single((0:1023)'/5.12); %monotonic 5.12 GSPS time vector (1024 samples in 200 ns)
    for i=1:ne
       X(:,i)=interp1cfloat(that(i,:),X(:,i),timem);
    end
end

