function [normdata,normcolumn] = normCx(depth,data,numberchans)
%
%
%
%
%
%% Parse the inputs
% Parameters TBD

%%
normcolumn = [0:1/numberchans:1-1/numberchans];

%
newdepth = [];
newdata = [];
for i = 1:size(depth,3)
    newdepth = cat(2,newdepth,depth(:,:,i));
    newdata = cat(2,newdata,data(:,:,i));
end

%
[B,I] = sort(newdepth);
newdata = newdata(:,I(~isnan(B)));
B = B(~isnan(B));
[C,ia,ic] = unique(B);

%
mndata = NaN(size(newdata,1),length(C));
for i = 1:max(ic)
    mndata(:,i) = nanmean(newdata(:,find(ic == i)),2);
end

%
normdata = NaN(size(newdata,1),length(normcolumn));
for ii = 1:size(mndata,1)
    normdata(ii,:) = interp1(C,mndata(ii,:),normcolumn,'nearest');
end

% Spatial smoothening
spat_sm = 20;
for t = 1:size(normdata,1)
    normdata(t,:) = smooth(normdata(t,:),spat_sm,'lowess');
end

end
