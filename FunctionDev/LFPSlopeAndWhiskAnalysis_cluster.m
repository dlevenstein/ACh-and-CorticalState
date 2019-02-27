function [WhiskPSScorr,PSSstatsdepth,PSS,pupildist,...
    pupilPSSdist,pupilEMGdist,PSSEMGdist,whiskPETH,...
    timelockedPSS] = LFPSlopeAndWhiskAnalysis_cluster( basePath,figfolder )
%UNTITLED3 Summary of this function goes here
%
%
%
%%
baseName = bz_BasenameFromBasepath(basePath);
sessionInfo = bz_getSessionInfo(basePath,'noPrompts',true);

%% Load Pupil
pupildilation = bz_LoadBehavior(basePath,'pupildiameter');

smoothwin = 0.5; %s
pupildilation.data = smooth(pupildilation.data,smoothwin.*pupildilation.samplingRate,'moving');
nantimes = isnan(pupildilation.data);
pupildilation.data = pupildilation.data(~isnan(pupildilation.data));

if length(pupildilation.data) < 1
    warning('Not enough pupil data >)'); 
    return
end

smoothwin = 2; %s
pupildilation.dpdt = diff(smooth(pupildilation.data,smoothwin.*pupildilation.samplingRate,'moving')).*pupildilation.samplingRate;
pupildilation.dpdt = smooth(pupildilation.dpdt,smoothwin.*pupildilation.samplingRate,'moving');
pupildilation.timestamps = pupildilation.timestamps(~nantimes);

%% Load the PSS
%load([basePath,filesep,baseName,'.PowerSpectrumSlope.lfp.mat'])

%% Load Whisks
EMGwhisk = bz_LoadBehavior(basePath,'EMGwhisk');

%% Get the pupil phase at each point in time
lowfilter = [0.01 0.1];
pupil4filter = pupildilation;

lowpupildata = bz_Filter(pupil4filter,'passband',lowfilter,'filter' ,'fir1','order',3);
%highpupildata = bz_Filter(pupil4filter,'passband',highfilter,'filter' ,'fir1');

%Get the pupil phase/power of each whisk start
EMGwhisk.phase = interp1(lowpupildata.timestamps,lowpupildata.phase,...
    EMGwhisk.ints.Wh(:,1),'nearest');
EMGwhisk.power = interp1(lowpupildata.timestamps,log10(lowpupildata.amp),...
    EMGwhisk.ints.Wh(:,1),'nearest');
%Get the closest PSS timepoint to each whisk... (dt issue...)

%%
% figure
% plot(EMGwhisk.phase,EMGwhisk.power,'.')
% hold on
% plot(EMGwhisk.phase+2*pi,EMGwhisk.power,'.')

%% Wh-PSS relation by channel
WhiskPSScorr.corr = zeros(size(sessionInfo.AnatGrps.Channels));
WhiskPSScorr.pup = zeros(size(sessionInfo.AnatGrps.Channels));
WhiskPSScorr.dpdt = zeros(size(sessionInfo.AnatGrps.Channels));
WhiskPSScorr.phasecoupling = zeros(size(sessionInfo.AnatGrps.Channels));

nbins = 50;
PSSstatsdepth.bins = linspace(-2,0,nbins);
PSSstatsdepth.dist = zeros(length(sessionInfo.AnatGrps.Channels),nbins);

for cc = 1:length(sessionInfo.AnatGrps.Channels)
    cc
    channum = sessionInfo.AnatGrps.Channels(cc);
    %channum = 31;
    WhiskPSScorr.channum(cc) = channum;
    WhiskPSScorr.chanpos(cc) = cc;
    
    PSSstatsdepth.channum(cc) = channum;
    PSSstatsdepth.chanpos(cc) = cc;
    %%
    lfp = bz_GetLFP(channum,'basepath',basePath,'noPrompts',true);
    
    %%
    dt = 0.2;
    winsize = 1;
    [PSS] = bz_PowerSpectrumSlope(lfp,winsize,dt,'showfig',false);
    
    PSSstatsdepth.dist(cc,:) = hist(PSS.data,PSSstatsdepth.bins);
    PSSstatsdepth.dist(cc,:)./sum(PSSstatsdepth.dist(cc,:));
    
    %%
    PSS.EMG = interp1(EMGwhisk.timestamps,EMGwhisk.EMGenvelope,PSS.timestamps);
    PSS.pupilsize = interp1(pupildilation.timestamps,pupildilation.data,...
        PSS.timestamps,'nearest');
    PSS.dpdt = interp1(pupildilation.timestamps(1:end-1),pupildilation.dpdt,...
        PSS.timestamps,'nearest');
    PSS.pupilphase = interp1(lowpupildata.timestamps,lowpupildata.phase,...
        PSS.timestamps,'nearest');
    
    %%
    [WhiskPSScorr.EMG(cc),WhiskPSScorr.EMG_p(cc)] =...
        corr(log10(PSS.EMG),PSS.data,...
        'type','spearman','rows','complete');
    [WhiskPSScorr.pup(cc),WhiskPSScorr.pup_p(cc)] =...
        corr(log10(PSS.pupilsize),PSS.data,...
        'type','spearman','rows','complete');
    [WhiskPSScorr.dpdt(cc),WhiskPSScorr.dpdt_p(cc)] =...
        corr(PSS.dpdt,PSS.data,...
        'type','spearman','rows','complete');
    WhiskPSScorr.phasecoupling(cc) = ...
        abs(nanmean((PSS.data./nanmean(PSS.data)).*exp(1i.*PSS.pupilphase)));
    
    clear lfp
    
end

%%
[~, bestchans.EMG] = max(WhiskPSScorr.EMG);
bestchans.EMG = WhiskPSScorr.channum(bestchans.EMG);

[~, bestchans.pup] = max(WhiskPSScorr.phasecoupling);
bestchans.pup = WhiskPSScorr.channum(bestchans.pup);

%%
figure;

subplot(2,3,2:3);
plot(WhiskPSScorr.EMG,-WhiskPSScorr.chanpos,'b','linewidth',2)
hold on
plot(WhiskPSScorr.pup,-WhiskPSScorr.chanpos,'k','linewidth',2)
plot(WhiskPSScorr.dpdt,-WhiskPSScorr.chanpos,'k--','linewidth',1)
legend('EMG','Pupil Area','dpdt','location','eastoutside')
xlabel('PSS Correlation');
axis tight
%xlim([0 0.6])
box off

subplot(2,3,1);
imagesc(PSSstatsdepth.bins,PSSstatsdepth.chanpos,PSSstatsdepth.dist)
xlabel('PSS')
ylabel('Channel by Depth')

subplot(2,3,4);
plot(WhiskPSScorr.phasecoupling,-WhiskPSScorr.chanpos,'k','linewidth',2)
xlabel('Phase Coupling');
axis tight
%xlim([0 0.6])
box off

NiceSave('PSSCorrbyDepth',figfolder,baseName)

%% Finding best params for PSS
lfp = bz_GetLFP(bestchans.EMG,'basepath',basePath,'noPrompts',true);

dt = 0.05;
winsize = 1;
[PSS] = bz_PowerSpectrumSlope(lfp,winsize,dt,'showfig',true);

%%
PSS.EMG = interp1(EMGwhisk.timestamps,EMGwhisk.EMGenvelope,PSS.timestamps);
PSS.pupphase = interp1(lowpupildata.timestamps,lowpupildata.phase,...
    PSS.timestamps,'nearest');
PSS.pupmag = interp1(lowpupildata.timestamps,log10(lowpupildata.amp),...
    PSS.timestamps);

PSS.pupthresh = nanmedian(PSS.pupmag);% -0.75;
PSS.highpup = PSS.pupmag>PSS.pupthresh;

%%
pupildist.edges = {linspace(-pi,pi,20),linspace(-1.5,0,30)};
[pupildist.counts,pupildist.bins] = hist3([PSS.pupphase,PSS.pupmag],...
    'Edges',pupildist.edges);
pupildist.joint = pupildist.counts./sum(pupildist.counts(:));

pupildist.conditional = bsxfun(@rdivide,...
    pupildist.counts,sum(pupildist.counts,2));

cosx = linspace(-pi,3*pi,100);

figure;

subplot(2,2,1);
hist(PSS.pupmag)
hold on
plot(PSS.pupthresh.*[1 1],get(gca,'ylim'))

subplot(2,2,2);
imagesc(pupildist.bins{1},pupildist.bins{2},pupildist.joint')
hold on
imagesc(pupildist.bins{1}+2*pi,pupildist.bins{2},pupildist.joint')
plot(cosx,cos(cosx)./3-1.1,'w','linewidth',2)
plot([-pi 3*pi],PSS.pupthresh.*[1 1],'k--')
xlim([-pi 3*pi])
axis xy
xlabel('Pupil Phase');ylabel('Pupil Power')

subplot(2,2,3);
imagesc(pupildist.bins{1},pupildist.bins{2},pupildist.conditional')
hold on
imagesc(pupildist.bins{1}+2*pi,pupildist.bins{2},pupildist.conditional')
plot(cosx,cos(cosx)./3-1.1,'w','linewidth',2)
plot([-pi 3*pi],PSS.pupthresh.*[1 1],'k--')
xlim([-pi 3*pi])
axis xy
xlabel('Pupil Phase');ylabel('Pupil Power')

NiceSave('PSSPhasePow',figfolder,baseName)

%% Distribution of PSS given pupil phase
pupilPSSdist.edges = {linspace(-pi,pi,40),linspace(-2,0,50)};
[pupilPSSdist.counts,pupilPSSdist.bins] = hist3([PSS.pupphase,PSS.data],...
    'Edges',pupilPSSdist.edges);
pupilPSSdist.joint = pupilPSSdist.counts./sum(pupilPSSdist.counts(:));

pupilPSSdist.conditional = bsxfun(@rdivide,...
    pupilPSSdist.counts,sum(pupilPSSdist.counts,2));

[pupilPSSdist.counts_high] = hist3([PSS.pupphase(PSS.highpup),PSS.data(PSS.highpup)],...
    'Edges',pupilPSSdist.edges);
pupilPSSdist.conditional_high = bsxfun(@rdivide,...
    pupilPSSdist.counts_high,sum(pupilPSSdist.counts_high,2));

[pupilPSSdist.counts_low] = hist3([PSS.pupphase(~PSS.highpup),PSS.data(~PSS.highpup)],...
    'Edges',pupilPSSdist.edges);
pupilPSSdist.counts_low = bsxfun(@rdivide,...
    pupilPSSdist.counts_low,sum(pupilPSSdist.counts_low,2));

%% Figure: PSS by pupil phase
figure;

subplot(2,2,1);
imagesc(pupilPSSdist.bins{1},pupilPSSdist.bins{2},pupilPSSdist.joint')
hold on
imagesc(pupilPSSdist.bins{1}+2*pi,pupilPSSdist.bins{2},pupilPSSdist.joint')
plot(cosx,cos(cosx)./3-1.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
xlabel('Pupil Phase');ylabel('PSS')

subplot(2,2,3);
imagesc(pupilPSSdist.bins{1},pupilPSSdist.bins{2},pupilPSSdist.conditional')
hold on
imagesc(pupilPSSdist.bins{1}+2*pi,pupilPSSdist.bins{2},pupilPSSdist.conditional')
plot(cosx,cos(cosx)./3-1.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
xlabel('Pupil Phase');ylabel('PSS')

subplot(3,3,6);
imagesc(pupilPSSdist.bins{1},pupilPSSdist.bins{2},pupilPSSdist.counts_high')
hold on
imagesc(pupilPSSdist.bins{1}+2*pi,pupilPSSdist.bins{2},pupilPSSdist.counts_high')
plot(cosx,cos(cosx)./3-1.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
xlabel('Pupil Phase');ylabel('PSS')

subplot(3,3,9);
imagesc(pupilPSSdist.bins{1},pupilPSSdist.bins{2},pupilPSSdist.counts_low')
hold on
imagesc(pupilPSSdist.bins{1}+2*pi,pupilPSSdist.bins{2},pupilPSSdist.counts_low')
plot(cosx,cos(cosx)./3-1.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
xlabel('Pupil Phase');ylabel('PSS')

NiceSave('PSSbyPupilPhase',figfolder,baseName)

%% Distribution of EMG given pupil phase
pupilEMGdist.edges = {linspace(-pi,pi,40),linspace(-2,1,50)};
[pupilEMGdist.counts,pupilEMGdist.bins] = hist3([PSS.pupphase,log10(PSS.EMG)],...
    'Edges',pupilEMGdist.edges);
pupilEMGdist.joint = pupilEMGdist.counts./sum(pupilEMGdist.counts(:));

pupilEMGdist.conditional = bsxfun(@rdivide,...
    pupilEMGdist.counts,sum(pupilEMGdist.counts,2));

[pupilEMGdist.counts_high] = hist3([PSS.pupphase(PSS.highpup),log10(PSS.EMG(PSS.highpup))],...
    'Edges',pupilEMGdist.edges);
pupilEMGdist.conditional_high = bsxfun(@rdivide,...
    pupilEMGdist.counts_high,sum(pupilEMGdist.counts_high,2));

[pupilEMGdist.counts_low] = hist3([PSS.pupphase(~PSS.highpup),log10(PSS.EMG(~PSS.highpup))],...
    'Edges',pupilEMGdist.edges);
pupilEMGdist.counts_low = bsxfun(@rdivide,...
    pupilEMGdist.counts_low,sum(pupilEMGdist.counts_low,2));

%% Figure: EMG by pupil phase
cosx = linspace(-pi,3*pi,100);

figure;

subplot(2,2,1);
imagesc(pupilEMGdist.bins{1},pupilEMGdist.bins{2},pupilEMGdist.joint')
hold on
imagesc(pupilEMGdist.bins{1}+2*pi,pupilEMGdist.bins{2},pupilEMGdist.joint')
plot([-pi 3*pi],log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],'r--')
plot(cosx,cos(cosx)./3-1.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
xlabel('Pupil Phase');ylabel('EMG')
title('P(EMG,phase)')

subplot(2,2,3);
imagesc(pupilEMGdist.bins{1},pupilEMGdist.bins{2},pupilEMGdist.conditional')
hold on
imagesc(pupilEMGdist.bins{1}+2*pi,pupilEMGdist.bins{2},pupilEMGdist.conditional')
plot(cosx,cos(cosx)./3-1.6,'w','linewidth',2)
plot([-pi 3*pi],log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],'r--')
xlim([-pi 3*pi])
axis xy
xlabel('Pupil Phase');ylabel('EMG')
title('P(EMG|phase)')

subplot(3,3,6);
imagesc(pupilEMGdist.bins{1},pupilEMGdist.bins{2},pupilEMGdist.counts_high')
hold on
imagesc(pupilEMGdist.bins{1}+2*pi,pupilEMGdist.bins{2},pupilEMGdist.counts_high')
plot([-pi 3*pi],log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],'r--')
plot(cosx,cos(cosx)./3-1.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
xlabel('Pupil Phase');ylabel('EMG')
%title('High Pupil')

subplot(3,3,9);
imagesc(pupilEMGdist.bins{1},pupilEMGdist.bins{2},pupilEMGdist.counts_low')
hold on
imagesc(pupilEMGdist.bins{1}+2*pi,pupilEMGdist.bins{2},pupilEMGdist.counts_low')
plot([-pi 3*pi],log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],'r--')
plot(cosx,cos(cosx)./3-1.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
xlabel('Pupil Phase');ylabel('EMG')
%title('Low Pupil')

NiceSave('EMGbyPupilPhase',figfolder,baseName)

%% EMGPSS
PSSEMGdist.edges = {linspace(-2,1,50),linspace(-2,0,50)};
[PSSEMGdist.counts,PSSEMGdist.bins] = hist3([log10(PSS.EMG),PSS.data],...
    'Edges',PSSEMGdist.edges);
PSSEMGdist.joint = PSSEMGdist.counts./sum(PSSEMGdist.counts(:));

%%
figure;

subplot(2,2,1);
imagesc(PSSEMGdist.edges{1},PSSEMGdist.edges{2},PSSEMGdist.joint')
hold on
plot(log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],get(gca,'ylim'),'r--')
axis xy
xlabel('EMG');ylabel('PSS')

NiceSave('EMGPSS',figfolder,baseName)

%%
figure;

subplot(2,2,1);
plot(log10(PSS.EMG),PSS.data,'k.')
xlabel('EMG');ylabel('PSS');

subplot(2,2,2);
scatter(PSS.pupphase,PSS.data,3,log10(PSS.EMG))
hold on
scatter(PSS.pupphase+2*pi,PSS.data,3,log10(PSS.EMG))
colorbar
xlabel('Pupil Phase');ylabel('PSS')
% subplot(2,2,3)
%     plot(PSS.pupmag,PSS.data,'k.')
%     xlabel('Pupil Magnitude');ylabel('PSS')

subplot(2,2,4);
scatter(PSS.pupphase(~PSS.highpup),log10(PSS.EMG(~PSS.highpup)),2,PSS.data(~PSS.highpup))
hold on
scatter(PSS.pupphase(~PSS.highpup)+2*pi,log10(PSS.EMG(~PSS.highpup)),2,PSS.data(~PSS.highpup))
colorbar
ylim([-2 1]);xlim([-pi 3*pi])
caxis([-1.5 0])
xlabel('Pupil Phase');ylabel('EMG')

subplot(2,2,3);
scatter(PSS.pupphase(PSS.highpup),log10(PSS.EMG(PSS.highpup)),2,PSS.data(PSS.highpup))
hold on
scatter(PSS.pupphase(PSS.highpup)+2*pi,log10(PSS.EMG(PSS.highpup)),2,PSS.data(PSS.highpup))
colorbar
caxis([-1.5 0])
ylim([-2 1]);xlim([-pi 3*pi])
xlabel('Pupil Phase');ylabel('EMG')

NiceSave('PupilEMGPSS',figfolder,baseName)

%%
EMGwhisk.dur = diff(EMGwhisk.ints.Wh,1,2);
EMGwhisk.longwhisks = EMGwhisk.dur>1;

%% Example Whisks
% winsize = 40;
% exwhisk = randsample(EMGwhisk.ints.Wh(EMGwhisk.longwhisks,1),1);
% viewwin = exwhisk + winsize.*[-1 1];
% 
% figure
% subplot(3,1,1)
% imagesc(PSS.timestamps,log10(PSS.freqs),PSS.specgram)
% hold on
% plot(PSS.timestamps,PSS.data+2.2,'w','Linewidth',2)
% 
% axis xy
% xlim(viewwin)
% LogScale('y',10)
% 
% subplot(6,1,3)
% plot(pupildilation.timestamps,pupildilation.data,'k','LineWidth',2)
% hold on
% ylim([0 2.5])
% plot(EMGwhisk.timestamps,EMGwhisk.EMG./20,'b')
% box off
% xlim(viewwin)
% ylabel('Pupil')

% subplot(6,1,5)
% plot(EMGwhisk.timestamps,EMGwhisk.EMG,'k')
% hold on
% plot(EMGwhisk.timestamps,EMGwhisk.EMGenvelope,'b','linewidth',2)
% axis tight
% box off
% xlim(viewwin)

%%
%Find the high pupil/long whisk in the example window and the low
%pupil/short whisk in the example window
winsize = 3;
exwhisk = randsample(EMGwhisk.ints.Wh(EMGwhisk.longwhisks,1),1);
viewwin = exwhisk + winsize.*[-1 1];

figure;

subplot(3,1,1);
imagesc(PSS.timestamps,log10(PSS.freqs),PSS.specgram)
axis xy
xlim(viewwin)

subplot(6,1,3);
plot(lfp.timestamps,lfp.data)
xlim(viewwin)

subplot(6,1,4);
plot(PSS.timestamps,PSS.data)
xlim(viewwin)

subplot(6,1,5);
plot(EMGwhisk.timestamps,EMGwhisk.EMG,'k')
hold on
plot(EMGwhisk.timestamps,EMGwhisk.EMGenvelope,'b')
xlim(viewwin)

NiceSave('ExWhisk',figfolder,baseName)

%% Get PSS around Whisks
EMGwhisk.numwhisks = length(EMGwhisk.ints.Wh(:,1));
EMGwhisk.highpupil = EMGwhisk.power>PSS.pupthresh;
whiskPETH.window = [-5 5]; %s
whiskPETH.windex = 5*PSS.samplingRate; %s
whiskPETH.timestamps = whiskPETH.window(1):(1/PSS.samplingRate):whiskPETH.window(2);
timelockedPSS.data = zeros(length(whiskPETH.timestamps),EMGwhisk.numwhisks);

for ww = 1:EMGwhisk.numwhisks
    PSS.whidx(ww) = find(PSS.timestamps==interp1(PSS.timestamps,PSS.timestamps,EMGwhisk.ints.Wh(ww,1),'nearest'));
    
    timelockedPSS.data(:,ww) = PSS.data(PSS.whidx(ww)-whiskPETH.windex:PSS.whidx(ww)+whiskPETH.windex);
    timelockedPSS.timestamps(:,ww) = whiskPETH.timestamps;
    timelockedPSS.phases(:,ww) = ones(size(whiskPETH.timestamps)).*EMGwhisk.phase(ww);
    timelockedPSS.highpupil(:,ww) = true(size(whiskPETH.timestamps)).*EMGwhisk.highpupil(ww);
    %NEED TO TAKE OUT OTHER WHISKS!?
    %     prevwhisk = EMGwhisk.ints.Wh(ww-1,2) - EMGwhisk.ints.Wh(ww,1);
    %     nextwhisk = EMGwhisk.ints.Wh(ww+1,1) - EMGwhisk.ints.Wh(ww,2);
    %     timelockedPSS.otherwhisks(:,ww) = ...
    %         timelockedPSS.timestamps(:,ww)<prevwhisk |...
    %         timelockedPSS.timestamps(:,ww)>nextwhisk;
    
end
timelockedPSS.highpupil = logical(timelockedPSS.highpupil); %Why?

[phasePETH.high]=PairMatHist(timelockedPSS.data(timelockedPSS.highpupil),...%&~timelockedPSS.otherwhisks),...%
    [timelockedPSS.timestamps(timelockedPSS.highpupil),...%&~timelockedPSS.otherwhisks),...
    timelockedPSS.phases(timelockedPSS.highpupil)],...%&~timelockedPSS.otherwhisks)],...
    60,[-pi 4]);

[phasePETH.low]=PairMatHist(timelockedPSS.data(~timelockedPSS.highpupil),...%&~timelockedPSS.otherwhisks),...
    [timelockedPSS.timestamps(~timelockedPSS.highpupil),...%&~timelockedPSS.otherwhisks),...
    timelockedPSS.phases(~timelockedPSS.highpupil)],...%&~timelockedPSS.otherwhisks)],...
    60,[-pi 5]);

%%
figure;

subplot(2,2,1);
imagesc(phasePETH.high.bincenters,phasePETH.high.bincenters,phasePETH.high.mean')
hold on
imagesc(phasePETH.high.bincenters,phasePETH.high.bincenters+2*pi,phasePETH.high.mean')
plot(EMGwhisk.dur(EMGwhisk.highpupil),EMGwhisk.phase(EMGwhisk.highpupil),'r.','markersize',1)
plot(EMGwhisk.dur(EMGwhisk.highpupil),EMGwhisk.phase(EMGwhisk.highpupil)+2*pi,'r.','markersize',1)
plot(cos(cosx),cosx,'w')
plot([0 0],[-pi 3*pi],'r')
xlim([-1 4]);ylim([-pi 3*pi])
colorbar
axis xy
caxis([-1.4 -0.5])
xlabel('t (s, aligned to Wh Onset)');ylabel('Pupil Phase')
title('High Power Pupil')

subplot(2,2,2);
imagesc(phasePETH.low.bincenters,phasePETH.low.bincenters,phasePETH.low.mean')
hold on
imagesc(phasePETH.low.bincenters,phasePETH.low.bincenters+2*pi,phasePETH.low.mean')
plot(EMGwhisk.dur(~EMGwhisk.highpupil),EMGwhisk.phase(~EMGwhisk.highpupil),'r.','markersize',1)
plot(EMGwhisk.dur(~EMGwhisk.highpupil),EMGwhisk.phase(~EMGwhisk.highpupil)+2*pi,'r.','markersize',1)
plot(cos(cosx),cosx,'w')
plot([0 0],[-pi 3*pi],'r')
xlim([-1 4]);ylim([-pi 3*pi])
colorbar
axis xy
caxis([-1.4 -0.5])
xlabel('t (s, aligned to Wh Onset)');ylabel('Pupil Phase')
title('Low Power Pupil')

NiceSave('PETHbyPhase',figfolder,baseName)

%%
figure;

subplot(2,2,1);
scatter(timelockedPSS.timestamps(timelockedPSS.highpupil),...
    timelockedPSS.phases(timelockedPSS.highpupil),2,timelockedPSS.data(timelockedPSS.highpupil))
hold on
scatter(timelockedPSS.timestamps(timelockedPSS.highpupil),...
    timelockedPSS.phases(timelockedPSS.highpupil)+2*pi,2,timelockedPSS.data(timelockedPSS.highpupil))

plot(EMGwhisk.dur(EMGwhisk.highpupil),EMGwhisk.phase(EMGwhisk.highpupil),'r.','markersize',1)
plot(EMGwhisk.dur(EMGwhisk.highpupil),EMGwhisk.phase(EMGwhisk.highpupil)+2*pi,'r.','markersize',1)
plot(cos(cosx),cosx,'w')
xlim([-1 5]);ylim([-pi 3*pi])
colorbar
axis xy
caxis([-1.5 0])

subplot(2,2,2);
scatter(timelockedPSS.timestamps(~timelockedPSS.highpupil),...
    timelockedPSS.phases(~timelockedPSS.highpupil),2,timelockedPSS.data(~timelockedPSS.highpupil))
hold on
scatter(timelockedPSS.timestamps(~timelockedPSS.highpupil),...
    timelockedPSS.phases(~timelockedPSS.highpupil)+2*pi,2,timelockedPSS.data(~timelockedPSS.highpupil))

plot(EMGwhisk.dur(~EMGwhisk.highpupil),EMGwhisk.phase(~EMGwhisk.highpupil),'r.','markersize',1)
plot(EMGwhisk.dur(~EMGwhisk.highpupil),EMGwhisk.phase(~EMGwhisk.highpupil)+2*pi,'r.','markersize',1)
plot(cos(cosx),cosx,'w')
xlim([-1 5]);ylim([-pi 3*pi])
colorbar
caxis([-1.5 0])
axis xy

NiceSave('TimelockedPSS',figfolder,baseName)

%% Whisk Sorts
[~,whisksorts.phase] = sort(EMGwhisk.phase);
[~,whisksorts.dur] = sort(EMGwhisk.dur);

%%
figure;

subplot(2,2,1);
imagesc(whiskPETH.timestamps,[1 EMGwhisk.numwhisks],timelockedPSS.data')
hold on
plot([0 0],[1 EMGwhisk.numwhisks],'b')
plot(EMGwhisk.dur,1:EMGwhisk.numwhisks,'.r','markersize',1)
xlim([-2 5])
colorbar
%caxis([-1.5 -0.5])

subplot(2,2,2);
imagesc(whiskPETH.timestamps,[1 EMGwhisk.numwhisks],timelockedPSS.data(:,whisksorts.phase)')
hold on
plot(EMGwhisk.dur(whisksorts.phase),1:EMGwhisk.numwhisks,'.r','markersize',1)
plot([0 0],[0 1],'b')
xlim([-2 5])
colorbar

subplot(2,2,3);
imagesc(whiskPETH.timestamps,[1 EMGwhisk.numwhisks],timelockedPSS.data(:,whisksorts.dur)')
hold on
plot(EMGwhisk.dur(whisksorts.dur),1:EMGwhisk.numwhisks,'.r','markersize',1)
plot([0 0],[1 EMGwhisk.numwhisks],'b')
xlim([-2 5])
colorbar

NiceSave('PSSallWhisks',figfolder,baseName)

%% OLD STUFF BELOW

%% Take a look at the channels with dpdt and p
% [~, bestchans.pup] = max(WhiskPSScorr.pup);
% [~, bestchans.dpdt] = max(WhiskPSScorr.dpdt);
% repchans = [sessionInfo.AnatGrps.Channels(bestchans.pup) sessionInfo.AnatGrps.Channels(bestchans.dpdt)];
% %repchans = 42;
% lfp = bz_GetLFP(sessionInfo.AnatGrps.Channels(repchans(1)),'basepath',basePath,'noPrompts',true);
% 
% dt = 0.2;
% winsize = 2;
% [PSS,specgram] = bz_PowerSpectrumSlope(lfp,winsize,dt,'showfig',true,...
%     'saveMat',basePath);
% 
% PSS.pupilsize = interp1(pupildilation.timestamps,pupildilation.data,...
%     PSS.timestamps,'nearest');
% PSS.dpdt = interp1(pupildilation.timestamps(1:end-1),pupildilation.dpdt,...
%     PSS.timestamps,'nearest');
% 
% %Check residuals
% 
% %%
% numbins = 15;
% bins = linspace(-0.5,0.5,numbins+1);
% pupcyclePSS.bincenters = bins(1:end-1)+0.5.*diff(bins([1 2]));
% bins([1 end])=[-Inf Inf];
% [N,~,~,BINX,BINY] = histcounts2(log10(PSS.pupilsize),PSS.dpdt,...
%     bins,bins);
% pupcyclePSS.meanPSS = zeros(size(N));
% for xx = 1:numbins
%     for yy = 1:numbins
%         pupcyclePSS.meanPSS(xx,yy) = nanmean(PSS.data(BINX==xx & BINY==yy));
%     end
% end
% 
% nbinthresh = 10;  %Must have more than 10 time windows
% pupcyclePSS.meanPSS(N<nbinthresh) = nan;
% 
% %% PSS and UP/DOWN
% SlowWaves = bz_LoadEvents(basePath,'SlowWaves');
% 
% %%
% updown = {'DOWN','UP'};
% UDcolor = {'b','r'};
% for ss = 1:2
%     SlowWaves.dur.(updown{ss}) = diff(SlowWaves.ints.(updown{ss}),1,2);
%     SlowWaves.midpoint.(updown{ss}) = mean(SlowWaves.ints.(updown{ss}),2);
%     SlowWaves.PSS.(updown{ss}) = interp1(PSS.timestamps,PSS.data,SlowWaves.midpoint.(updown{ss}));
% end
% 
% %%
% numbins = 30;
% PSShist.bins = linspace(-2,0,numbins);
% 
% PSShist.hist = hist(PSS.data,PSShist.bins);
% 
% %%
% exwinsize = 300;
% exwin = bz_RandomWindowInIntervals(PSS.timestamps([1 end])',exwinsize);
% 
% winsize = 8;
% maxtime = pupildilation.timestamps(...
%     (pupildilation.timestamps>exwin(1) & pupildilation.timestamps<exwin(2)) & pupildilation.data==...
%     max(pupildilation.data(pupildilation.timestamps>exwin(1) & pupildilation.timestamps<exwin(2))));
% subsamplewin = maxtime(1)+[-0.5 0.5].*winsize;
% 
% mintime = pupildilation.timestamps(...
%     (pupildilation.timestamps>exwin(1) & pupildilation.timestamps<exwin(2)) & pupildilation.data==...
%     min(pupildilation.data(pupildilation.timestamps>exwin(1) & pupildilation.timestamps<exwin(2))));
% subsamplewin2 = mintime(1)+[-0.5 0.5].*winsize;
% 
% figure
% subplot(4,1,1)
% imagesc(specgram.timestamps,log2(specgram.freqs),specgram.amp)
% axis xy
% xlim(exwin)
% ylabel({'Specgram','f (Hz)'})
% LogScale('y',2)
% 
% subplot(6,1,3)
% plot(PSS.timestamps,PSS.data,'k')
% axis tight
% xlim(exwin)
% box off
% xlabel('t (s)')
% ylabel('PSS')
% subplot(6,1,4)
% plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2)
% hold on
% plot(pupildilation.timestamps(1:end-1),pupildilation.dpdt,'k')
% axis tight
% xlim(exwin)
% box off
% ylabel('Pupil')
% plot(exwin,[0 0],'k--')
% plot(subsamplewin,3.*ones(size(subsamplewin)),'r','linewidth',2)
% plot(subsamplewin2,3.*ones(size(subsamplewin2)),'r','linewidth',2)
% 
% subplot(6,2,11)
% bz_MultiLFPPlot(lfp,'timewin',subsamplewin)
% hold on
% % plot(slow.timestamps,10000.*slow.data,'g','linewidth',1)
% ylabel('High Pupil')
% xlim(subsamplewin)
% 
% subplot(6,2,12)
% bz_MultiLFPPlot(lfp,'timewin',subsamplewin2)
% hold on
% % plot(slow.timestamps,10000.*slow.data,'g','linewidth',1)
% ylabel('Low Pupil')
% xlim(subsamplewin2)
% 
% subplot(6,2,9)
% plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2)
% hold on
% plot(pupildilation.timestamps(1:end-1),pupildilation.dpdt,'k')
% plot(PSS.timestamps,PSS.data)
% axis tight
% xlim(subsamplewin)
% box off
% ylabel('Pupil')
% plot(subsamplewin,[0 0],'k--')
% 
% subplot(6,2,10)
% plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2)
% hold on
% plot(pupildilation.timestamps(1:end-1),pupildilation.dpdt,'k')
% plot(PSS.timestamps,PSS.data)
% axis tight
% xlim(subsamplewin2)
% box off
% ylabel('Pupil')
% plot(subsamplewin2,[0 0],'k--')
% 
% NiceSave('PSSexample',figfolder,baseName)
% 
% %%
% figure
% subplot(3,3,1)
% h = imagesc(pupcyclePSS.bincenters,pupcyclePSS.bincenters,pupcyclePSS.meanPSS');
% set(h,'AlphaData',~isnan(pupcyclePSS.meanPSS'));
% hold on
% plot(pupcyclePSS.bincenters([1 end]),[0 0],'k--')
% LogScale('x',10)
% axis xy
% colorbar
% xlabel('Pupil Area (med^-^1)')
% ylabel('dp/dt')
% LogScale('x',10)
% title('PSS')
% 
% subplot(6,3,12)
% plot(PSShist.bins,PSShist.hist,'k','linewidth',2)
% box off
% xlabel('PSS')
% axis tight
% xlim([-1.75 -0.25])
% 
% subplot(3,2,2)
% for ss = 1:2
%     plot(SlowWaves.PSS.(updown{ss}),log10(SlowWaves.dur.(updown{ss})),'.','color',UDcolor{ss},'markersize',3)
%     hold on
% end
% xlabel('PSS');ylabel('Dur (s)')
% axis tight
% box off
% LogScale('y',10)
% legend(updown{:},'location','eastoutside')
% 
% subplot(3,3,4)
% plot(log10(PSS.pupilsize),PSS.data,'k.','markersize',1)
% xlabel('Pupil Area (med^-^1)');ylabel('PSS')
% box off
% axis tight
% subplot(3,3,5)
% plot((PSS.dpdt),PSS.data,'k.','markersize',1)
% xlabel('dpdt (med^-^1s^-^1)');ylabel('PSS')
% box off
% axis tight
% 
% NiceSave('PSSandUPDOWN',figfolder,baseName)
% 
% %%
% samplewin = bz_RandomWindowInIntervals(lfp.timestamps([1 end])',300);
% winsize = 10;
% maxtime = pupildilation.timestamps(...
%     (pupildilation.timestamps>samplewin(1) & pupildilation.timestamps<samplewin(2)) & pupildilation.data==...
%     max(pupildilation.data(pupildilation.timestamps>samplewin(1) & pupildilation.timestamps<samplewin(2))));
% subsamplewin = maxtime(1)+[-0.5 0.5].*winsize;
% 
% mintime = pupildilation.timestamps(...
%     (pupildilation.timestamps>samplewin(1) & pupildilation.timestamps<samplewin(2)) & pupildilation.data==...
%     min(pupildilation.data(pupildilation.timestamps>samplewin(1) & pupildilation.timestamps<samplewin(2))));
% subsamplewin2 = mintime(1)+[-0.5 0.5].*winsize;
% 
% figure
% % subplot(3,3,1)
% %     plot(slopepupilcorr.dpdt,'k','linewidth',0.5)
% %     hold on
% %     plot(slopepupilcorr.p,'k','linewidth',2)
% %     plot(repchans(1),slopepupilcorr.dpdt(repchans(1)),'o')
% %     plot(repchans(2),slopepupilcorr.p(repchans(2)),'o')
% %     xlabel('Channel (sorted by depth)')
% %     axis tight
% %     box off
% 
% for cc = 1%:2
%     subplot(6,6,12+cc)
%     plot(log10(PSS.pupilsize),PSS.data(:,cc),'k.','markersize',3)
%     xlabel('Pupil Area (med^-^1)');ylabel('Spectrum Slope')
%     box off
%     axis tight
%     subplot(6,6,18+cc)
%     plot((PSS.dpdt),PSS.data(:,cc),'k.','markersize',1)
%     xlabel('dpdt (med^-^1s^-^1)');ylabel('Spectrum Slope')
%     box off
%     axis tight
% end
% 
% % subplot(6,1,4)
% %     bz_MultiLFPPlot(lfp,'timewin',samplewin)
% subplot(6,6,9:12)
% plot(PSS.timestamps,PSS.data)
% axis tight
% xlim(samplewin)
% box off
% xlabel('t (s)')
% ylabel('PSS')
% 
% subplot(6,6,3:6)
% plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2)
% hold on
% plot(pupildilation.timestamps(1:end-1),pupildilation.dpdt,'k')
% axis tight
% xlim(samplewin)
% box off
% ylabel('Pupil')
% plot(samplewin,[0 0],'k--')
% plot(subsamplewin,3.*ones(size(subsamplewin)),'r','linewidth',2)
% plot(subsamplewin2,3.*ones(size(subsamplewin2)),'r','linewidth',2)
% 
% 
% subplot(6,3,11:12)
% bz_MultiLFPPlot(lfp,'timewin',samplewin)
% hold on
% % plot(slow.timestamps,10000.*slow.data,'g','linewidth',1)
% ylabel('High Pupil')
% xlim(subsamplewin)
% 
% subplot(6,3,17:18)
% bz_MultiLFPPlot(lfp,'timewin',subsamplewin2)
% hold on
% % plot(slow.timestamps,10000.*slow.data,'g','linewidth',1)
% ylabel('Low Pupil')
% xlim(subsamplewin2)
% 
% subplot(6,3,8:9)
% plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2)
% hold on
% plot(pupildilation.timestamps(1:end-1),pupildilation.dpdt,'k')
% plot(PSS.timestamps,PSS.data)
% axis tight
% xlim(subsamplewin)
% box off
% ylabel('Pupil')
% plot(subsamplewin,[0 0],'k--')
% 
% subplot(6,3,14:15)
% plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2)
% hold on
% plot(pupildilation.timestamps(1:end-1),pupildilation.dpdt,'k')
% plot(PSS.timestamps,PSS.data)
% axis tight
% xlim(subsamplewin2)
% box off
% ylabel('Pupil')
% plot(subsamplewin2,[0 0],'k--')
% 
% NiceSave('PSSandPupil',figfolder,baseName)

end