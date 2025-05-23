% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [] = motionNTC()
ha=gca;
p=ha.CurrentPoint;  x1=p(1,:);  x2=p(2,:); %mouse pointing vector
if x1(1)<ha.XLim(1) || x1(1)>ha.XLim(2) || x1(2)<ha.YLim(1) || x1(2)>ha.YLim(2) || isempty(ha.Children); return; end
deleteh(findobj(gcf,'Tag','motion'));

for i=1:numel(ha.Children)  
    try 
        X{i} = [ha.Children(i).XData(:), ha.Children(i).YData(:), ha.Children(i).ZData(:)];  %#ok<*AGROW>
    catch
        X{i} = [ha.Children(i).XData(:), ha.Children(i).YData(:), zeros(size(ha.Children(i).YData(:)))];
    end

    [~, dx10] = fcnrange(x1, X{i});
    dx21=ones(size(dx10,1),1)*(x2-x1);
    r = fcnrange(cross(dx21,dx10))/norm(x2-x1);
    
    [mr(i),mi(i)]=min(r);
end
[~, i] = min(mr); mi=mi(i); X=X{i};

j=find(any(isnan(X),2))-mi;
if any(j)
    first = find(j<=0,1,'last');
    last = find(j>0,1,'first');
    j=(j(first):j(last))+mi;
else
    j=1:size(X,1);
end
plot3(X(j,1),X(j,2),X(j,3),'Color',ha.Children(i).Color','LineWidth',4,'Tag','motion');






% [~, dx10] = fcnrange(x1, x(v0,13:15));
% dx21=ones(size(dx10,1),1)*(x2-x1);
% r2 = fcnrange(cross(dx21,dx10))/norm(x2-x1);
% r = min(r1,r2);
% 
% [r,i]=min(r); i=v0(i);
% pid=x(:,7);
% upid = pid(i);
% if r>200
%    return 
% end
% 
% tid=x(:,8);
% utid = tid(i);
% v1=find(tid==utid);  nv1=numel(v1);  v11=v1(1);
% p = [x(v1,10:12); x(v1(end),13:15)];
% ps = [x(i,10:12); x(i,13:15)];
% t1 = x(v1,16);  t1i=x(i,16);
% 
% %parents
% name = particlename{v11};
% processname = processname{v11};
% parentID = x(v11,9);
% if parentID==0
%     parentname = 'Particle Gun';
% else
%     parentname = particlename{find(tid==parentID,1)};
% end
% energy = sum(x(v1,20));
% 
% rall=sum(fcnrange(x(v1,10:12)-x(v1,13:15)));
% r1=fcnrange(ps(1,:),ps(2,:));
% de1=x(i,20);
% ke1=x(i,18);
% 
% t1 = sort(t1,'ascend');
% nstep=find(t1i==t1,1);
% 
% 
% %children
% v1 = find(x(:,9)==utid);
% [tidc,i] = unique(tid(v1));  v1=v1(i);  nc=numel(tidc);
% upidc=unique(x(v1,7));  nuc=numel(upidc);
% str = sprintf(' %s %s, %.3gns, %.3gMeVke\n track %.0f (%.0f points, %.3gmm, %.3gMeV)\n step %.0f (%.3gmm, %.3gMeV)\n %s parent (track %.0f)', ...
%     strrep(name,'_','\_'),processname,t1i,ke1,utid,nv1,rall,energy,nstep,r1,de1,strrep(parentname,'_','\_'),parentID);
% 
% if nc>0
%     str = [str sprintf('\n %.0f children',nc)];
%     for j=1:nuc
%         p1 = upidc(j);
%         n = sum(x(v1,7)==p1);
%         str = [str sprintf(', %.0f %s',n,particlename{find(x(:,7)==p1,1)})];
%     end
% else
%    str = [str sprintf('\n no children')];
% end
% 
% h(1)=plot3(p(:,1),p(:,2),p(:,3),'.-','MarkerSize',20,'LineWidth',2,'Color',[1 .98 0]);
% h(2)=plot3(ps(:,1),ps(:,2),ps(:,3),'.-','MarkerSize',35,'LineWidth',3.5,'Color',[1 .8 0]);
% h(3)=text(x1(1),x1(2),x1(3), str, 'HorizontalAlignment','Left','VerticalAlignment','Bottom','fontsize',8,'BackgroundColor','none');
% 




