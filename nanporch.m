function porch = nanporch(V,threshold)
%Best practice is to remove spikes before subtracting porch!
%porch defined as baseline of every point before first crossing of threshold
if nargin==1
    porch=nanmedian(V,1);
else
    porch=zeros(size(V,1),1);
    for iter=1:2
        outliers = cummax(V>(threshold+porch),1);
        Vi=V;  Vi(outliers)=nan;  porch=nanmedian(Vi,1);
    end
end


% %STEP 2 (OPTIONAL): TAKE MEAN OF REMAINING INLIERS FOR EACH CHANNEL
% V=Vi;
% srl = 2; %sigma rejection level
% inliers = ~outliers;
% for iter=1:4
%     np=sum(inliers,1);
%     mu = sum(V,1,'omitnan')./np;
%     if iter<4
%         xms=(V-mu).^2;  sigma = sqrt((1./(np-1)).*sum(xms,1,'omitnan'));
%         inliers = inliers & xms<(srl*sigma).^2;
%         V(~inliers) = nan;
%     end
% end
% porch=mu;





