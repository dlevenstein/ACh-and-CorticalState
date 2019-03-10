%%
basePath = pwd;
baseName = bz_BasenameFromBasepath(basePath);

sessionInfo = bz_getSessionInfo(basePath,'noPrompts',true);
badchannels = sessionInfo.badchannels;
usechannels = sessionInfo.AnatGrps.Channels;
usechannels(ismember(usechannels,badchannels))=[];
channels = sessionInfo.channels;

figfolder = fullfile(basePath,'AnalysisFigures');
savefile = fullfile(basePath,[baseName,'.PSSBehaviorAnalysis.mat']);

%% Loading behavior...
% Pupil diameter
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

% EMG
EMGwhisk = bz_LoadBehavior(basePath,'EMGwhisk');

%Specifying SPONT whisking
load(fullfile(basePath,[baseName,'.MergePoints.events.mat']),'MergePoints');
sidx = find(startsWith(MergePoints.foldernames,"Spont"));
sponttimes = [MergePoints.timestamps(sidx(1),1) MergePoints.timestamps(sidx(end),2)];

spontidx = find(EMGwhisk.ints.Wh(:,2) < sponttimes(2));
EMGwhisk.ints.Wh = EMGwhisk.ints.Wh(spontidx,:);

spontidx = find(EMGwhisk.ints.NWh(:,2) < sponttimes(2));
EMGwhisk.ints.NWh = EMGwhisk.ints.NWh(spontidx,:);

spontidx = find(EMGwhisk.timestamps < sponttimes(2));
EMGwhisk.timestamps = EMGwhisk.timestamps(spontidx);
EMGwhisk.EMGenvelope = EMGwhisk.EMGenvelope(spontidx);
EMGwhisk.EMG = EMGwhisk.EMG(spontidx);
EMGwhisk.EMGsm = EMGwhisk.EMGsm(spontidx);

EMGwhisk.pupiltime = interp1(EMGwhisk.timestamps,EMGwhisk.EMGsm,pupildilation.timestamps,'nearest');

% Pupil phase
lowfilter = [0.01 0.1];
pupil4filter = pupildilation;
lowpupildata = bz_Filter(pupil4filter,'passband',lowfilter,'filter' ,'fir1','order',3);

% Get the pupil phase/power of each whisk start
EMGwhisk.phase = interp1(lowpupildata.timestamps,lowpupildata.phase,...
    EMGwhisk.ints.Wh(:,1),'nearest');
EMGwhisk.power = interp1(lowpupildata.timestamps,log10(lowpupildata.amp),...
    EMGwhisk.ints.Wh(:,1),'nearest');

%% Load PSS and aligning behavior
load(fullfile(basePath,[baseName,'.PowerSpectrumSlope.lfp.mat']));

depthinfo = rescaleCx(basePath);
L56idx = find(depthinfo.ndepth >= 0.6 & depthinfo.ndepth <= 0.9);
L56idx = depthinfo.channels(L56idx);

tempPSSEMGcorr = PSpecSlope.Shortwin.PSScorr.EMG(L56idx+1);
bestchan = find(tempPSSEMGcorr == max(tempPSSEMGcorr));

% PSS.data = PSpecSlope.Shortwin.PSS(:,L56idx(bestchan)+1);
% PSS.timestamps = PSpecSlope.Shortwin.timestamps;
% PSS.samplingRate = 1/mean(diff(PSS.timestamps));

% Re-running PSS w/ better time resolution
lfp = bz_GetLFP(L56idx(bestchan),'basepath',basePath,'noPrompts',true);
lfp.data = lfp.data(spontidx);
lfp.timestamps = lfp.timestamps(spontidx);

% Deconstruction
movingwin = round([0.5 0.125].*lfp.samplingRate);
nwin = floor((length(lfp.data) - movingwin(1))/movingwin(2));
Frange = [2.5, 100]; % define frequency range for power-law fitting

sig = zeros(movingwin(1),nwin);
timestamp = zeros(nwin,1);
for i = 1 : nwin
    idx = [ceil((i-1)*movingwin(2))+1 : ceil((i-1)*movingwin(2))+movingwin(1)];
    sig(:,i) = lfp.data(idx);
    %figure out timestamp associated with window i
    timestamp(i) = mean(lfp.timestamps(idx));
end
clear lfp

Frac = amri_sig_fractal_gpu(sig,lfp.samplingRate,'detrend',1);
Frac.timestamps = timestamp;
Frac = amri_sig_plawfit(Frac,Frange);

PSS.data = Frac.Beta.*-1;
PSS.timestamps = Frac.timestamps;
PSS.samplingRate = 1/mean(diff(PSS.timestamps));

PSS.EMG = interp1(EMGwhisk.timestamps,EMGwhisk.EMGenvelope,...
    PSS.timestamps,'nearest');
PSS.pupil = interp1(pupildilation.timestamps,pupildilation.data,...
    PSS.timestamps,'nearest');
PSS.dpdt = interp1(pupildilation.timestamps(1:end-1),pupildilation.dpdt,...
    PSS.timestamps,'nearest');
PSS.pupphase = interp1(lowpupildata.timestamps,lowpupildata.phase,...
    PSS.timestamps,'nearest');
PSS.amp = interp1(lowpupildata.timestamps,log10(lowpupildata.amp),...
    PSS.timestamps,'nearest');

PSS.pupthresh = nanmedian(PSS.pupil);
PSS.highpup = PSS.pupil>PSS.pupthresh;

%% PSS-Pupil phase and dPdt codist
numbins = 15;
bins = linspace(-0.5,0.5,numbins+1);
pupcyclePSS.bincenters = bins(1:end-1)+0.5.*diff(bins([1 2]));
bins([1 end])=[-Inf Inf];
[N,~,~,BINX,BINY] = histcounts2(log10(PSS.pupil),PSS.dpdt,...
    bins,bins);
pupcyclePSS.meanPSS = zeros(size(N));

for xx = 1:numbins
    for yy = 1:numbins
        pupcyclePSS.meanPSS(xx,yy) = nanmean(PSS.data(BINX==xx & BINY==yy));
    end
end

nbinthresh = 10;  %Must have more than 10 time windows
pupcyclePSS.meanPSS(N<nbinthresh) = NaN;

% PSS distribution
numbins = 30;
PSShist.bins = linspace(-4,0,numbins);
PSShist.hist = hist(PSS.data,PSShist.bins);

%% PSS and UP/DOWN
SlowWaves = bz_LoadEvents(basePath,'SlowWaves');

updown = {'DOWN','UP'};
UDcolor = {'b','r'};
for ss = 1:2
    SlowWaves.dur.(updown{ss}) = diff(SlowWaves.ints.(updown{ss}),1,2);
    SlowWaves.midpoint.(updown{ss}) = mean(SlowWaves.ints.(updown{ss}),2);
    SlowWaves.PSS.(updown{ss}) = interp1(PSS.timestamps,PSS.data,SlowWaves.midpoint.(updown{ss}));
end

%%
figure;
for ss = 1:2
    plot(SlowWaves.PSS.(updown{ss}),log10(SlowWaves.dur.(updown{ss})),'.','color',UDcolor{ss},'markersize',10)
    hold on
end
xlabel('PSS'); ylabel('UP/DOWN duration (s)');
LogScale('y',10)
axis square
legend(updown{:},'location','eastoutside')

%NiceSave('PSS_SlowWaves',figfolder,baseName)

%% FIGURE:
figure;
subplot(2,2,1)
bar(PSShist.bins,PSShist.hist,'b')
xlabel('PSS')
axis tight
xlim([-3 0])

subplot(2,2,2)
h = imagesc(pupcyclePSS.bincenters,pupcyclePSS.bincenters,pupcyclePSS.meanPSS');
set(h,'AlphaData',~isnan(pupcyclePSS.meanPSS'));
hold on
plot(pupcyclePSS.bincenters([1 end]),[0 0],'k--')
LogScale('x',10)
axis xy
colorbar
xlabel('Pupil Area (med^-^1)')
ylabel('dp/dt')
LogScale('x',10)
title('PSS')

subplot(2,2,3)
plot(log10(PSS.pupil),PSS.data,'k.','markersize',1)
xlabel('Pupil Area (med^-^1)');ylabel('PSS')
axis tight

subplot(2,2,4)
plot((PSS.dpdt),PSS.data,'k.','markersize',1)
xlabel('dpdt (med^-^1s^-^1)');ylabel('PSS')
axis tight

%NiceSave('PSS_Phase_dPdt',figfolder,baseName)

%%
pupildist.edges = {linspace(-pi,pi,20),linspace(-1.5,0,30)};
[pupildist.counts,pupildist.bins] = hist3([PSS.pupphase,PSS.amp],...
    'Edges',pupildist.edges);
pupildist.joint = pupildist.counts./sum(pupildist.counts(:));

pupildist.conditional = bsxfun(@rdivide,...
    pupildist.counts,sum(pupildist.counts,2));

%% FIGURE:
cosx = linspace(-pi,3*pi,100);

figure;
subplot(2,2,1);
hist(PSS.pupil); hold on;
plot(PSS.pupthresh.*[1 1],get(gca,'ylim'));
xlabel('Pupil diameter'); ylabel('counts (au)');

subplot(2,2,3);
imagesc(pupildist.bins{1},pupildist.bins{2},pupildist.joint')
hold on
imagesc(pupildist.bins{1}+2*pi,pupildist.bins{2},pupildist.joint')
plot(cosx,cos(cosx)./3-1.1,'w','linewidth',2)
plot([-pi 3*pi],PSS.pupthresh.*[1 1],'k--')
xlim([-pi 3*pi])
axis xy
xlabel('Pupil Phase');ylabel('Pupil Power')
title('P(Pupil power,phase)')

subplot(2,2,4);
imagesc(pupildist.bins{1},pupildist.bins{2},pupildist.conditional')
hold on
imagesc(pupildist.bins{1}+2*pi,pupildist.bins{2},pupildist.conditional')
plot(cosx,cos(cosx)./3-1.1,'w','linewidth',2)
plot([-pi 3*pi],PSS.pupthresh.*[1 1],'k--')
xlim([-pi 3*pi])
axis xy
xlabel('Pupil Phase');ylabel('Pupil Power')
title('P(Pupil power|phase)')

%NiceSave('PSS_PhasePow',figfolder,baseName)

%% Pupil phase-PSS codistribution
pupilPSSdist.edges = {linspace(-pi,pi,40),linspace(-4,0,50)};
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

%% FIGURE:
figure;
subplot(2,2,1);
imagesc(pupilPSSdist.bins{1},pupilPSSdist.bins{2},pupilPSSdist.joint')
hold on
imagesc(pupilPSSdist.bins{1}+2*pi,pupilPSSdist.bins{2},pupilPSSdist.joint')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
xlabel('Pupil Phase');ylabel('PSS')
title('P(PSS,phase)')

subplot(2,2,3);
imagesc(pupilPSSdist.bins{1},pupilPSSdist.bins{2},pupilPSSdist.conditional')
hold on
imagesc(pupilPSSdist.bins{1}+2*pi,pupilPSSdist.bins{2},pupilPSSdist.conditional')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
xlabel('Pupil phase');ylabel('PSS')
title('P(PSS|phase)')

subplot(2,2,2);
imagesc(pupilPSSdist.bins{1},pupilPSSdist.bins{2},pupilPSSdist.counts_high')
hold on
imagesc(pupilPSSdist.bins{1}+2*pi,pupilPSSdist.bins{2},pupilPSSdist.counts_high')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
xlabel('Pupil phase');ylabel('PSS')
title('>median pupil');

subplot(2,2,4);
imagesc(pupilPSSdist.bins{1},pupilPSSdist.bins{2},pupilPSSdist.counts_low')
hold on
imagesc(pupilPSSdist.bins{1}+2*pi,pupilPSSdist.bins{2},pupilPSSdist.counts_low')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
xlabel('Pupil Phase');ylabel('PSS')
title('<median pupil');

%NiceSave('PSSbyPupilPhase',figfolder,baseName)

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

%% FIGURE:
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
xlabel('Pupil phase');ylabel('EMG')
title('P(EMG,phase)')

subplot(2,2,3);
imagesc(pupilEMGdist.bins{1},pupilEMGdist.bins{2},pupilEMGdist.conditional')
hold on
imagesc(pupilEMGdist.bins{1}+2*pi,pupilEMGdist.bins{2},pupilEMGdist.conditional')
plot(cosx,cos(cosx)./3-1.6,'w','linewidth',2)
plot([-pi 3*pi],log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],'r--')
xlim([-pi 3*pi])
axis xy
xlabel('Pupil phase');ylabel('EMG')
title('P(EMG|phase)')

subplot(2,2,2);
imagesc(pupilEMGdist.bins{1},pupilEMGdist.bins{2},pupilEMGdist.counts_high')
hold on
imagesc(pupilEMGdist.bins{1}+2*pi,pupilEMGdist.bins{2},pupilEMGdist.counts_high')
plot([-pi 3*pi],log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],'r--')
plot(cosx,cos(cosx)./3-1.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
xlabel('Pupil phase');ylabel('EMG')
title('>median pupil');

subplot(2,2,4);
imagesc(pupilEMGdist.bins{1},pupilEMGdist.bins{2},pupilEMGdist.counts_low')
hold on
imagesc(pupilEMGdist.bins{1}+2*pi,pupilEMGdist.bins{2},pupilEMGdist.counts_low')
plot([-pi 3*pi],log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],'r--')
plot(cosx,cos(cosx)./3-1.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
xlabel('Pupil phase');ylabel('EMG')
title('<median pupil');

%NiceSave('EMGbyPupilPhase',figfolder,baseName)

%% EMG-PSS codistribution:
PSSEMGdist.edges = {linspace(-2,1,50),linspace(-3,0,50)};
[PSSEMGdist.counts,PSSEMGdist.bins] = hist3([log10(PSS.EMG),PSS.data],...
    'Edges',PSSEMGdist.edges);
PSSEMGdist.joint = PSSEMGdist.counts./sum(PSSEMGdist.counts(:));

PSSEMGdist.conditional = bsxfun(@rdivide,...
    PSSEMGdist.counts,sum(PSSEMGdist.counts,2));

%% FIGURE:
figure;
subplot(1,2,1);
imagesc(PSSEMGdist.edges{1},PSSEMGdist.edges{2},PSSEMGdist.joint')
hold on
plot(log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],get(gca,'ylim'),'r--')
axis square
xlabel('EMG');ylabel('PSS')
title('P(EMG,PSS)')

subplot(1,2,2);
imagesc(PSSEMGdist.edges{1},PSSEMGdist.edges{2},PSSEMGdist.conditional')
hold on
plot(log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],get(gca,'ylim'),'r--')
axis square
xlabel('EMG');ylabel('PSS')
title('P(EMG|PSS)')

%NiceSave('EMGPSS',figfolder,baseName)

%% FIGURE:
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

subplot(2,2,3);
scatter(PSS.pupphase(PSS.highpup),log10(PSS.EMG(PSS.highpup)),2,PSS.data(PSS.highpup))
hold on
scatter(PSS.pupphase(PSS.highpup)+2*pi,log10(PSS.EMG(PSS.highpup)),2,PSS.data(PSS.highpup))
colorbar
caxis([-1.5 0])
ylim([-2 1]);xlim([-pi 3*pi])
xlabel('Pupil Phase');ylabel('EMG')

subplot(2,2,4);
scatter(PSS.pupphase(~PSS.highpup),log10(PSS.EMG(~PSS.highpup)),2,PSS.data(~PSS.highpup))
hold on
scatter(PSS.pupphase(~PSS.highpup)+2*pi,log10(PSS.EMG(~PSS.highpup)),2,PSS.data(~PSS.highpup))
colorbar
ylim([-2 1]);xlim([-pi 3*pi])
caxis([-1.5 0])
xlabel('Pupil Phase');ylabel('EMG')

%NiceSave('PupilEMGPSS',figfolder,baseName)

%% Get PSS around Whisks
% Only LONG whisks selected...
EMGwhisk.dur = diff(EMGwhisk.ints.Wh,1,2);
EMGwhisk.longwhisks = EMGwhisk.dur>1;

EMGwhisk.numwhisks = length(EMGwhisk.ints.Wh(:,1));
EMGwhisk.highpupil = EMGwhisk.power>nanmedian(EMGwhisk.power);
whiskPETH.window = [-5 5]; %s
whiskPETH.windex = 5*PSS.samplingRate; %s
whiskPETH.timestamps = whiskPETH.window(1):(1/PSS.samplingRate):whiskPETH.window(2);
timelockedPSS.data = zeros(length(whiskPETH.timestamps),EMGwhisk.numwhisks);

for ww = 1:EMGwhisk.numwhisks
    PSS.whidx(ww) = find(PSS.timestamps==interp1(PSS.timestamps,PSS.timestamps,EMGwhisk.ints.Wh(ww,1),'nearest'));
    
    if PSS.whidx(ww)-whiskPETH.windex > 0 && PSS.whidx(ww)+whiskPETH.windex < length(PSS.data)
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
    
end
timelockedPSS.highpupil = logical(timelockedPSS.highpupil); %Why?

[phasePETH.high]=PairMatHist(timelockedPSS.data(timelockedPSS.highpupil),...%&~timelockedPSS.otherwhisks),...%
    [timelockedPSS.timestamps(timelockedPSS.highpupil),...%&~timelockedPSS.otherwhisks),...
    timelockedPSS.phases(timelockedPSS.highpupil)],...%&~timelockedPSS.otherwhisks)],...
    25,[-pi 4]);

[phasePETH.low]=PairMatHist(timelockedPSS.data(~timelockedPSS.highpupil),...%&~timelockedPSS.otherwhisks),...
    [timelockedPSS.timestamps(~timelockedPSS.highpupil),...%&~timelockedPSS.otherwhisks),...
    timelockedPSS.phases(~timelockedPSS.highpupil)],...%&~timelockedPSS.otherwhisks)],...
    25,[-pi 5]);

% Sorting Wh by phase/duration
[~,whisksorts.phase] = sort(EMGwhisk.phase);
[~,whisksorts.dur] = sort(EMGwhisk.dur);

% Saving to struct

%% FIGURE:
figure;
subplot(1,2,1);
imagesc(phasePETH.high.bincenters,phasePETH.high.bincenters,phasePETH.high.mean')
hold on
imagesc(phasePETH.high.bincenters,phasePETH.high.bincenters+2*pi,phasePETH.high.mean')
plot(EMGwhisk.dur(EMGwhisk.highpupil),EMGwhisk.phase(EMGwhisk.highpupil),'r.','markersize',2)
plot(EMGwhisk.dur(EMGwhisk.highpupil),EMGwhisk.phase(EMGwhisk.highpupil)+2*pi,'r.','markersize',2)
plot(cos(cosx),cosx,'w','linewidth',2)
plot([0 0],[-pi 3*pi],'r')
colorbar
axis xy
caxis([-3 -1])
xlim([-1 4]);ylim([-pi 3*pi])
xlabel('t (s, aligned to Wh Onset)');ylabel('Pupil Phase')
title('>median pupil')

subplot(1,2,2);
imagesc(phasePETH.low.bincenters,phasePETH.low.bincenters,phasePETH.low.mean'); hold on
imagesc(phasePETH.low.bincenters,phasePETH.low.bincenters+2*pi,phasePETH.low.mean')
plot(EMGwhisk.dur(~EMGwhisk.highpupil),EMGwhisk.phase(~EMGwhisk.highpupil),'r.','markersize',1)
plot(EMGwhisk.dur(~EMGwhisk.highpupil),EMGwhisk.phase(~EMGwhisk.highpupil)+2*pi,'r.','markersize',1)
plot(cos(cosx),cosx,'w','linewidth',2)
plot([0 0],[-pi 3*pi],'r')
colorbar
axis xy
caxis([-3 -1])
xlim([-1 4]);ylim([-pi 3*pi])
xlabel('t (s, aligned to Wh Onset)');ylabel('Pupil Phase')
title('<median pupil')

%NiceSave('PETHbyPhase',figfolder,baseName)

%% FIGURE:
figure;
subplot(1,2,1);
imagesc(whiskPETH.timestamps,[1 EMGwhisk.numwhisks],timelockedPSS.data(:,whisksorts.phase)')
hold on
plot(EMGwhisk.dur(whisksorts.phase),1:EMGwhisk.numwhisks,'.r','markersize',1)
plot([0 0],[0 1],'b')
axis square
xlim([-2 5])
colorbar; caxis([-3 -1])
xlabel('t (s, aligned to Wh Onset)');ylabel('trial no.')
title('Epochs sorted by pupil phase');

subplot(1,2,2);
imagesc(whiskPETH.timestamps,[1 EMGwhisk.numwhisks],timelockedPSS.data(:,whisksorts.dur)')
hold on
plot(EMGwhisk.dur(whisksorts.dur),1:EMGwhisk.numwhisks,'.r','markersize',1)
plot([0 0],[1 EMGwhisk.numwhisks],'b')
axis square
xlim([-2 5])
colorbar; caxis([-3 -1])
xlabel('t (s, aligned to Wh Onset)');ylabel('trial no.')
title('Epochs sorted by Wh duration');

%NiceSave('PSSallWhisks',figfolder,baseName)