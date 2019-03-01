function [normdata] = rescaleCx( basePath,analysisname )
%UNTITLED Summary of this function goes here
%
%
%
% WMunoz - 02/28/2019
%% Scales according to H3 probe
chandist = [0:20:1275];
lnorm = [0 0.1 0.35 0.5 0.6 0.75 0.9 1];

%%
[baseFolder,baseName] = fileparts(basePath);
load(fullfile(basePath,[baseName,'.',analysisname,'.lfp.mat']));

%%
load(fullfile(basePath,[baseName,'.sessionInfo.mat']));
channels = sessionInfo.channels;
usechannels = sessionInfo.AnatGrps.Channels;

%lborders = [NaN 20 18 62 16 36 33 32];
lborders = sessionInfo.layerborders;

%%
truelayer1 = chandist(usechannels == lborders(end))*lnorm(2);
chandist = chandist-(chandist(usechannels == lborders(2))-truelayer1); 
truecolumn = chandist(usechannels == lborders(end));
normdepth = NaN(1,length(usechannels));

for i = 2:length(lnorm)-1
    lb1 = find(usechannels == lborders(i));
    lb2 = find(usechannels == lborders(i+1));
    lfactor1 = lnorm(i);
    lfactor2 = lnorm(i+1);
    normdepth(lb1:lb2) = ((chandist(lb1:lb2) - chandist(lb1)).* (lfactor2 - lfactor1))./(chandist(lb2) - chandist(lb1)) + lfactor1;
end

%Dealing w/ pia:L1
if chandist(1) < 0
    lb1 = find(chandist >= 0);
else
    lb1 = 1;
end
lb2 = find(usechannels == lborders(2));
lfactor1 = lnorm(1);
lfactor2 = lnorm(2);
temp = (([0 chandist(lb1(1):lb2)] - 0).* (lfactor2 - lfactor1))./(chandist(lb2) - 0) + lfactor1;
normdepth(lb1(1):lb2) = temp(2:end);

%% CONSIDER BELOW!!!
%%
badchannels = sessionInfo.badchannels;
badidx = ismember(usechannels,badchannels);
usechannels(badidx) = [];
chandist(badidx) = [];
normdepth(badidx) = [];

%%
smlopspec = laminarpspec.LOdata;
lof = laminarpspec.LOfreqs;
wm = find(usechannels == lborders(end));

figure;
imagesc(log10(lof),normdepth(2:wm),smlopspec(:,2:wm)')
axis xy
LogScale('x',10)
colormap('jet'); caxis([min(min(smlopspec)) max(max(smlopspec))]);
c = colorbar;
c.Label.String = 'power (dB)';
set(gca,'YDir','reverse');
set(gca,'Xtick',[log10([0.5 1 5 10 25 50 100 250])]);
set(gca,'Xticklabel',{'0.5','1','5','10','25','50','100','250'});
xlim([log10([0.5 250])]);
xtickangle(45);
xlabel('frequency (Hz)');
set(gca,'Ytick',lnorm);
set(gca,'Yticklabel',{'pia','L1/2','L3/4','L4/5a','L5a-5b','L5b-6a','L6a-6b','WM'});
set(gca,'YGrid','on','Layer','top','GridColor',[0 0 0]);
ylim([0 1]);

%%
nanchans = isnan(normdepth);
normdepth = normdepth(~nanchans);

%%
numberchans = 30;
normcolumn = [0:1/numberchans:1];

[bincounts,ind]=histc(normdepth,normcolumn);

normpspec = NaN(size(smlopspec,1),ind(end));
for i = 1:ind(end) 
    tempidx = find(ind == i);
    normpspec(:,i) = nanmean(smlopspec(:,tempidx),2);  
end
    
figure;
imagesc(log10(lof),normcolumn,normpspec')
axis xy
LogScale('x',10)
colormap('jet'); caxis([min(min(normpspec)) max(max(normpspec))]);
c = colorbar;
c.Label.String = 'power (dB)';
set(gca,'YDir','reverse');
set(gca,'Xtick',[log10([0.5 1 5 10 25 50 100 250])]);
set(gca,'Xticklabel',{'0.5','1','5','10','25','50','100','250'});
xlim([log10([0.5 250])]);
xtickangle(45);
xlabel('frequency (Hz)');
set(gca,'Ytick',lnorm);
set(gca,'Yticklabel',{'pia','L1/2','L3/4','L4/5a','L5a-5b','L5b-6a','L6a-6b','WM'});
set(gca,'YGrid','on','Layer','top','GridColor',[0 0 0]);
ylim([0 1]);
    

%%
zsmlopspec1 = zeros(size(normpspec1,1),size(normpspec1,2));
for n = 1:size(normpspec1,1)
    zsmlopspec1(n,:) = (normpspec1(n,:)-nanmean(normpspec1(n,:)))./nanstd(normpspec1(n,:));
end

zsmlopspec = zeros(size(normpspec,1),size(normpspec,2));
for n = 1:size(normpspec,1)
    zsmlopspec(n,:) = (normpspec(n,:)-nanmean(normpspec(n,:)))./nanstd(normpspec(n,:));
end

collectspec = cat(3,zsmlopspec,zsmlopspec1);

collectspec = cat(3,normpspec,normpspec1);
mnormpspec = nanmean(collectspec,3);

figure;
imagesc(log10(lof),normcolumn,mnormpspec')
axis xy
LogScale('x',10)
colormap('jet'); caxis([min(min(mnormpspec)) max(max(mnormpspec))]);
c = colorbar;
c.Label.String = 'power (dB)';
set(gca,'YDir','reverse');
set(gca,'Xtick',[log10([0.5 1 5 10 25 50 100 250])]);
set(gca,'Xticklabel',{'0.5','1','5','10','25','50','100','250'});
xlim([log10([0.5 250])]);
xtickangle(45);
xlabel('frequency (Hz)');
set(gca,'Ytick',lnorm);
set(gca,'Yticklabel',{'pia','L1/2','L3/4','L4/5a','L5a-5b','L5b-6a','L6a-6b','WM'});
set(gca,'YGrid','on','Layer','top','GridColor',[0 0 0]);
ylim([0 1]);

