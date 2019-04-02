%%
basePath = pwd;
baseName = bz_BasenameFromBasepath(basePath);

sessionInfo = bz_getSessionInfo(basePath,'noPrompts',true);
badchannels = sessionInfo.badchannels;
usechannels = sessionInfo.AnatGrps.Channels;
usechannels(ismember(usechannels,badchannels))=[];
channels = sessionInfo.channels;

figfolder = fullfile(basePath,'AnalysisFigures');
savefile = fullfile(basePath,[baseName,'.PSSBehaviorAnalysis2.mat']);

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

%% PUPIL dilation ONsets/OFFsets
% Set the Pupil dilation thresholds by troughs 
tempPup = pupildilation.dpdt;
tempPup(tempPup<0) = nan;
%Pupz = NormToInt(tempPup,'modZ'); %Modified Z score - robust to outliers
Pupz = (tempPup - nanmean(tempPup))./nanstd(tempPup,0,1);

% find by "gradient descent"(ish) from initial guess (0.5)
tPupbins = linspace(-1.5,2,100);
tPuphist = hist(log10(Pupz),tPupbins);
Pupgrad = smooth(gradient(tPuphist),4);

% Find troughs (gradient crossing from - to +)
troughidx = find(diff(Pupgrad<0)==1);
troughs = 10.^tPupbins(troughidx);

pupthreshold = troughs(find(troughs>0.5,1,'first'));
pup_thresh = Pupz > pupthreshold;
pup_on = round(find(pup_thresh(2:end) > pup_thresh(1:end-1))+1); %Pup onsets (si)
pup_off = round(find(pup_thresh(2:end) < pup_thresh(1:end-1))+1); %Pup offsets (si)

% If data starts/ends in the middle of an epoch, drop first/last trigger
if pup_off(1)<pup_on(1)
    pup_off = pup_off(2:end);
end
if pup_off(end) < pup_on(end)
    pup_on = pup_on(1:end-1);
end

dPpeak = NaN(length(pup_on),1);
for i = 1:length(pup_on)
    tidx = find(pupildilation.dpdt(pup_on(i):pup_off(i))...
        == max(pupildilation.dpdt(pup_on(i):pup_off(i))),1,'first');
    dPpeak(i) = pupildilation.dpdt(pup_on(i)+tidx);
    pup_on(i) = pup_on(i)+tidx;
end

% figure; 
% plot(pupildilation.dpdt,'k'); hold on;
% plot(pup_on,pupildilation.dpdt(pup_on),'ro')

%temppup_on = pup_on;
for i = 1:length(pup_on)
    tidx = find(pupildilation.dpdt(1:pup_on(i)) < 0,1,'last');
    if ~isempty(tidx)
        pup_on(i) = tidx;
    else
    end
end

% figure; 
% plot(pupildilation.dpdt,'k'); hold on;
% plot(temppup_on,pupildilation.dpdt(temppup_on),'ro')
% plot(pup_on,pupildilation.dpdt(pup_on),'go')

for i = 1:length(pup_on)
    tidx = find(pupildilation.data(pup_on(i):end-1) < pupildilation.data(pup_on(i)),1,'first');
    if ~isempty(tidx)
        pup_off(i) = pup_on(i)+tidx;
    else
        pup_off(i) = length(pupildilation.data)-1;
    end
end

% Convert to seconds, merge and exclude transients
pup_on = pupildilation.timestamps(pup_on);
pup_off = pupildilation.timestamps(pup_off);
% [ Pupints ] = MergeSeparatedInts( [pup_on,pup_off],0.1 );
Pupints = [pup_on,pup_off];
% [pup_on,pup_off] = MinEpochLength(Pupints(:,1),Pupints(:,2),0.1,1);
% Pupints = [pup_on,pup_off];

% figure;
% plot(pupildilation.timestamps,pupildilation.data,'k'); hold on;
% plot(Pupints',nanmean(pupildilation.data).*ones(size(Pupints))','g','linewidth',2)

% Durations and peaks
Pupdur = pup_off-pup_on;
pup_peak = length(pup_on);
for i = 1:length(pup_on)
    tsidx = find(pupildilation.timestamps == pup_on(i));
    teidx = find(pupildilation.timestamps == pup_off(i));
    pup_peak(i) = max(pupildilation.data(tsidx:teidx));
end

% figure; 
% scatter(log10(Pupdur),pup_peak,'k.');
% LogScale('x',10);

% Get the pupil phase/power of each whisk start
Pupon.phase = interp1(lowpupildata.timestamps,lowpupildata.phase,...
    pup_on,'nearest');
Pupon.power = interp1(lowpupildata.timestamps,log10(lowpupildata.amp),...
    pup_on,'nearest');

%% Load PSS and aligning behavior
%load(fullfile(basePath,[baseName,'.PowerSpectrumSlope.mat']));
load(fullfile(basePath,[baseName,'.PowerSpectrumSlope.lfp.mat']));

depthinfo = rescaleCx(basePath);
L56idx = find(depthinfo.ndepth >= 0.6 & depthinfo.ndepth <= 0.9);
L56idx = depthinfo.channels(L56idx);

tempPSSEMGcorr = PSpecSlope.Shortwin.PSScorr.EMG(L56idx+1);
bestchan = find(tempPSSEMGcorr == max(tempPSSEMGcorr));

PSS.data = PSpecSlope.Shortwin.PSS(:,L56idx(bestchan)+1);
PSS.timestamps = PSpecSlope.Shortwin.timestamps;
PSS.samplingRate = 1/mean(diff(PSS.timestamps));

% Re-running PSS w/ better time resolution
% lfp = bz_GetLFP(L56idx(bestchan),'basepath',basePath,'noPrompts',true);
% lfp.data = lfp.data(spontidx);
% lfp.timestamps = lfp.timestamps(spontidx);
% 
% movingwin = round([0.5 0.125].*lfp.samplingRate);
% nwin = floor((length(lfp.data) - movingwin(1))/movingwin(2));
% Frange = [2.5, 100]; % define frequency range for power-law fitting
% 
% sig = zeros(movingwin(1),nwin);
% timestamp = zeros(nwin,1);
% for i = 1 : nwin
%     idx = [ceil((i-1)*movingwin(2))+1 : ceil((i-1)*movingwin(2))+movingwin(1)];
%     sig(:,i) = lfp.data(idx);
%     %figure out timestamp associated with window i
%     timestamp(i) = mean(lfp.timestamps(idx));
% end
% 
% Frac = amri_sig_fractal_gpu(sig,lfp.samplingRate,'detrend',1);
% Frac.timestamps = timestamp;
% Frac = amri_sig_plawfit(Frac,Frange);
% 
% PSS.data = Frac.Beta.*-1;
% PSS.timestamps = Frac.timestamps;
% PSS.samplingRate = 1/mean(diff(PSS.timestamps));
% clear lfp Frac

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

[ ~,~,~,~,~, binsig,threshsig] = PhaseAmpCouplingByAmp2(PSS.pupphase,PSS.amp,...
PSS.EMG,'showFig',true,'AmpBounds',[-1.25 -0.25],...
'shufflesig',true,'AmpZNorm',false,'numAmpbins',10);

%the two posible thresholds
%binsig.sigthresh   %more sensitive to spurious significant bins
%threshsig.sigthresh  %more conservataive

PSS.pupthresh = binsig.sigthresh;
%PSS.pupthresh = threshsig.sigthresh;
PSS.highpup = PSS.amp>PSS.pupthresh; 

%% PSS-Pupil phase and dPdt codist
% numbins = 15;
% bins = linspace(-0.5,0.5,numbins+1);
% pupcyclePSS.bincenters = bins(1:end-1)+0.5.*diff(bins([1 2]));
% bins([1 end])=[-Inf Inf];
% [N,~,~,BINX,BINY] = histcounts2(log10(PSS.pupil),PSS.dpdt,...
%     bins,bins);
% pupcyclePSS.meanPSS = zeros(size(N));
% 
% for xx = 1:numbins
%     for yy = 1:numbins
%         pupcyclePSS.meanPSS(xx,yy) = nanmean(PSS.data(BINX==xx & BINY==yy));
%     end
% end
% 
% nbinthresh = 10;  %Must have more than 10 time windows
% pupcyclePSS.meanPSS(N<nbinthresh) = NaN;

% PSS distribution
numbins = 30;
PSShist.bins = linspace(-4,0,numbins);
PSShist.hist = hist(PSS.data,PSShist.bins);
PSShist.hist = PSShist.hist./sum(PSShist.hist);

% PSS-Pupil dynamics histo3
PSSPupdynhist.bins={linspace(-3,-0.5,75) linspace(-0.5,0.5,75)};
PSSPupdynhist.Pupcounts = hist3([PSS.data,log10(PSS.pupil)],...
    PSSPupdynhist.bins);
PSSPupdynhist.Pupcounts = PSSPupdynhist.Pupcounts./sum(PSSPupdynhist.Pupcounts(:));
PSSPupdynhist.dPcounts = hist3([PSS.data,PSS.dpdt],...
    PSSPupdynhist.bins);
PSSPupdynhist.dPcounts = PSSPupdynhist.dPcounts./sum(PSSPupdynhist.dPcounts(:));

% PSS Pup cycle conditional histo
[pupcyclePSS.meanZ,pupcyclePSS.N,pupcyclePSS.Xbins,...
    pupcyclePSS.Ybins ] = ConditionalHist3( log10(PSS.pupil),PSS.dpdt,...
    PSS.data,'minXY',0,'Xbounds',[-0.5 0.5],'Ybounds',[-0.5 0.5],...
    'numXbins',100,'numYbins',100);

% Saving to struct
PSSBehavior.PSShist = PSShist;
PSSBehavior.PSSPupdynhist = PSSPupdynhist;
PSSBehavior.pupcyclePSS = pupcyclePSS;

%% PSS and UP/DOWN
SlowWaves = bz_LoadEvents(basePath,'SlowWaves');

updown = {'DOWN','UP'};
UDcolor = {'b','r'};
for ss = 1:2
    SlowWaves.dur.(updown{ss}) = diff(SlowWaves.ints.(updown{ss}),1,2);
    SlowWaves.midpoint.(updown{ss}) = mean(SlowWaves.ints.(updown{ss}),2);
    SlowWaves.PSS.(updown{ss}) = interp1(PSS.timestamps,PSS.data,SlowWaves.midpoint.(updown{ss}));
end

UPdur.bins = linspace(log10(0.01),log10(100),50);
UPdur.hist = hist(log10(SlowWaves.dur.UP),UPdur.bins);
UPdur.hist = UPdur.hist./sum(UPdur.hist);

DOWNdur.bins = linspace(log10(0.01),log10(100),50);
DOWNdur.hist = hist(log10(SlowWaves.dur.DOWN),DOWNdur.bins);
DOWNdur.hist = DOWNdur.hist./sum(DOWNdur.hist);

PSSUPDOWNhist.bins={linspace(-3,-0.5,50) linspace(log10(0.01),log10(100),50)};
PSSUPDOWNhist.UPcounts = hist3([SlowWaves.PSS.UP,log10(SlowWaves.dur.UP)],...
    PSSUPDOWNhist.bins);
PSSUPDOWNhist.UPcounts = PSSUPDOWNhist.UPcounts./sum(PSSUPDOWNhist.UPcounts(:));
PSSUPDOWNhist.DOWNcounts = hist3([SlowWaves.PSS.DOWN,log10(SlowWaves.dur.DOWN)],...
    PSSUPDOWNhist.bins);
PSSUPDOWNhist.DOWNcounts = PSSUPDOWNhist.DOWNcounts./sum(PSSUPDOWNhist.DOWNcounts(:));

% Saving to struct
PSSBehavior.UPdur = UPdur;
PSSBehavior.DOWNdur = DOWNdur;
PSSBehavior.PSSUPDOWNhist = PSSUPDOWNhist;

%% FIGURE:
% figure;
% for ss = 1:2
%     plot(SlowWaves.PSS.(updown{ss}),log10(SlowWaves.dur.(updown{ss})),'.','color',UDcolor{ss},'markersize',10)
%     hold on
% end
% xlabel('PSS'); ylabel('UP/DOWN duration (s)');
% LogScale('y',10)
% axis square
% legend(updown{:},'location','eastoutside')

figure;
subplot(2,2,1);
bar(UPdur.bins,UPdur.hist,'r'); hold on;
bar(DOWNdur.bins,DOWNdur.hist,'b'); hold on;
xlabel('UP/DOWN duration (s)');
LogScale('x',10);
ylabel('norm counts');

subplot(2,2,2);
imagesc(PSSUPDOWNhist.bins{1},PSSUPDOWNhist.bins{2},PSSUPDOWNhist.UPcounts'); hold on;
colormap(gca,'autumn')
ColorbarWithAxis([min(min(PSSUPDOWNhist.UPcounts)) max(max(PSSUPDOWNhist.UPcounts))/4],['counts (au)'])
xlabel('PSS (au)'); ylabel('UP duration (s)')
LogScale('y',10);
axis square
axis xy

subplot(2,2,4);
imagesc(PSSUPDOWNhist.bins{1},PSSUPDOWNhist.bins{2},PSSUPDOWNhist.DOWNcounts'); hold on;
colormap(gca,'winter')
ColorbarWithAxis([min(min(PSSUPDOWNhist.DOWNcounts)) max(max(PSSUPDOWNhist.DOWNcounts))/4],['counts (au)'])
xlabel('PSS (au)'); ylabel('DOWN duration (s)')
LogScale('y',10);
axis square
axis xy

NiceSave('PSS_SlowWaves',figfolder,baseName)

%% FIGURE:
figure;
subplot(2,2,1)
bar(PSShist.bins,PSShist.hist,'k')
xlabel('PSS'); ylabel('norm. count');
axis tight
xlim([-3 0])

subplot(2,2,3)    
h = imagesc(pupcyclePSS.Xbins,pupcyclePSS.Ybins,pupcyclePSS.meanZ');
colormap(gca,'jet')
set(h,'AlphaData',~isnan(pupcyclePSS.meanZ'));
hold on
plot(pupcyclePSS.Xbins([1 end]),[0 0],'k--')
ColorbarWithAxis([min(min(pupcyclePSS.meanZ)) max(max(pupcyclePSS.meanZ))],['PSS (au)'])
LogScale('x',10)
axis xy
colorbar
xlabel('Pupil Area (med^-^1)')
ylabel('dp/dt')

subplot(2,2,2)
imagesc(PSSPupdynhist.bins{1},PSSPupdynhist.bins{2},log10(PSSPupdynhist.Pupcounts)'); hold on;
colormap(gca,'jet')
%ColorbarWithAxis([min(min(log10(PSSPupdynhist.Pupcounts))) max(max(log10(PSSPupdynhist.Pupcounts)))/2],['counts (au)'])
LogScale('c',10)
xlabel('PSS (au)'); ylabel('Pupil Area (med^-^1)')
axis square
axis xy

subplot(2,2,4)
imagesc(PSSPupdynhist.bins{1},PSSPupdynhist.bins{2},log10(PSSPupdynhist.dPcounts)'); hold on;
colormap(gca,'jet')
%ColorbarWithAxis([min(min(log10(PSSPupdynhist.dPcounts))) max(max(log10(PSSPupdynhist.dPcounts)))/2],['counts (au)'])
LogScale('c',10)
xlabel('PSS (au)'); ylabel('dPdt (med^-^1s^-^1)')
ylim([-0.4 0.4])
axis square
axis xy

NiceSave('PSS_Phase_dPdt',figfolder,baseName)

%%
pupildist.edges = {linspace(-pi,pi,20),linspace(-1.5,0,30)};
[pupildist.counts,pupildist.bins] = hist3([PSS.pupphase,PSS.amp],...
    'Edges',pupildist.edges);
pupildist.joint = pupildist.counts./sum(pupildist.counts(:));

pupildist.conditional = bsxfun(@rdivide,...
    pupildist.counts,sum(pupildist.counts,2));

% Saving to struct
PSSBehavior.pupildist = pupildist;

%% FIGURE:
cosx = linspace(-pi,3*pi,100);

figure;
subplot(2,2,1);
hist(PSS.amp); hold on;
plot(PSS.pupthresh.*[1 1],get(gca,'ylim'));
xlim([-3 1]);
xlabel('Filt-Pupil amplitude'); ylabel('counts (au)');

subplot(2,2,3);
imagesc(pupildist.bins{1},pupildist.bins{2},pupildist.joint')
hold on
imagesc(pupildist.bins{1}+2*pi,pupildist.bins{2},pupildist.joint')
plot(cosx,cos(cosx)./3-1.1,'w','linewidth',2)
plot([-pi 3*pi],PSS.pupthresh.*[1 1],'k--')
xlim([-pi 3*pi])
axis xy
xlabel('Pupil phase');ylabel('Pupil amp')
title('P(Pupil amp,phase)')

subplot(2,2,4);
imagesc(pupildist.bins{1},pupildist.bins{2},pupildist.conditional')
hold on
imagesc(pupildist.bins{1}+2*pi,pupildist.bins{2},pupildist.conditional')
plot(cosx,cos(cosx)./3-1.1,'w','linewidth',2)
plot([-pi 3*pi],PSS.pupthresh.*[1 1],'k--')
xlim([-pi 3*pi])
axis xy
xlabel('Pupil phase');ylabel('Pupil amp')
title('P(Pupil amp|phase)')

NiceSave('PSS_PhasePow',figfolder,baseName)

%% Pupil phase-PSS codistribution
%Can replace with ConditionalHist here
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
pupilPSSdist.conditional_low = bsxfun(@rdivide,...
    pupilPSSdist.counts_low,sum(pupilPSSdist.counts_low,2));

% sorting nonwhisking/lohipup
tempidx = interp1(PSS.timestamps,PSS.timestamps,EMGwhisk.ints.Wh,'nearest');
PSS.whisktimes = zeros(length(PSS.timestamps),1);
for i = 1:size(EMGwhisk.ints.Wh,1)
    temps = find(PSS.timestamps == tempidx(i,1));
    tempe = find(PSS.timestamps == tempidx(i,2));
    
    if ~isempty(temps)
        if ~isempty(tempe)
        PSS.whisktimes(temps:tempe) = 1;
        else
        end
    else
    end
end

[pupilPSSdist.counts_highnow] = hist3([PSS.pupphase(logical(PSS.highpup.*~PSS.whisktimes)),...
    PSS.data(logical(PSS.highpup.*~PSS.whisktimes))],...
    'Edges',pupilPSSdist.edges);
pupilPSSdist.conditional_highnow = bsxfun(@rdivide,...
    pupilPSSdist.counts_highnow,sum(pupilPSSdist.counts_highnow,2));

[pupilPSSdist.counts_lownow] = hist3([PSS.pupphase(logical(~PSS.highpup.*~PSS.whisktimes)),...
    PSS.data(logical(~PSS.highpup.*~PSS.whisktimes))],...
    'Edges',pupilPSSdist.edges);
pupilPSSdist.conditional_lownow = bsxfun(@rdivide,...
    pupilPSSdist.counts_lownow,sum(pupilPSSdist.counts_lownow,2));

[pupilPSSdist.counts_highw] = hist3([PSS.pupphase(logical(PSS.highpup.*PSS.whisktimes)),...
    PSS.data(logical(PSS.highpup.*PSS.whisktimes))],...
    'Edges',pupilPSSdist.edges);
pupilPSSdist.conditional_highw = bsxfun(@rdivide,...
    pupilPSSdist.counts_highw,sum(pupilPSSdist.counts_highw,2));

[pupilPSSdist.counts_loww] = hist3([PSS.pupphase(logical(~PSS.highpup.*PSS.whisktimes)),...
    PSS.data(logical(~PSS.highpup.*PSS.whisktimes))],...
    'Edges',pupilPSSdist.edges);
pupilPSSdist.conditional_loww = bsxfun(@rdivide,...
    pupilPSSdist.counts_loww,sum(pupilPSSdist.counts_loww,2));

% Saving to struct
PSSBehavior.pupilPSSdist = pupilPSSdist;

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
imagesc(pupilPSSdist.bins{1},pupilPSSdist.bins{2},pupilPSSdist.conditional_high')
hold on
imagesc(pupilPSSdist.bins{1}+2*pi,pupilPSSdist.bins{2},pupilPSSdist.conditional_high')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
xlabel('Pupil phase');ylabel('PSS')
title('sig pupil amp');

subplot(2,2,4);
imagesc(pupilPSSdist.bins{1},pupilPSSdist.bins{2},pupilPSSdist.conditional_low')
hold on
imagesc(pupilPSSdist.bins{1}+2*pi,pupilPSSdist.bins{2},pupilPSSdist.conditional_low')
plot(cosx,0.2.*cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
xlabel('Pupil Phase');ylabel('PSS')
title('ns pupil amp');

NiceSave('PSSbyPupilPhase',figfolder,baseName)

%% FIGURE:
figure;
subplot(2,2,1);
imagesc(pupilPSSdist.bins{1},pupilPSSdist.bins{2},pupilPSSdist.conditional_highw')
hold on
imagesc(pupilPSSdist.bins{1}+2*pi,pupilPSSdist.bins{2},pupilPSSdist.conditional_highw')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
xlabel('Pupil phase');ylabel('PSS')
title('sig pupil amp, whisky times');

subplot(2,2,3);
imagesc(pupilPSSdist.bins{1},pupilPSSdist.bins{2},pupilPSSdist.conditional_loww')
hold on
imagesc(pupilPSSdist.bins{1}+2*pi,pupilPSSdist.bins{2},pupilPSSdist.conditional_loww')
plot(cosx,0.2.*cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
xlabel('Pupil Phase');ylabel('PSS')
title('ns pupil amp, whisky times');

subplot(2,2,2);
imagesc(pupilPSSdist.bins{1},pupilPSSdist.bins{2},pupilPSSdist.conditional_highnow')
hold on
imagesc(pupilPSSdist.bins{1}+2*pi,pupilPSSdist.bins{2},pupilPSSdist.conditional_highnow')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
xlabel('Pupil phase');ylabel('PSS')
title('sig pupil amp, no whisky times');

subplot(2,2,4);
imagesc(pupilPSSdist.bins{1},pupilPSSdist.bins{2},pupilPSSdist.conditional_lownow')
hold on
imagesc(pupilPSSdist.bins{1}+2*pi,pupilPSSdist.bins{2},pupilPSSdist.conditional_lownow')
plot(cosx,0.2.*cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
xlabel('Pupil Phase');ylabel('PSS')
title('ns pupil amp, no whisky times');

NiceSave('PSSbyPupilPhase_Whisksep',figfolder,baseName)

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
pupilEMGdist.conditional_low = bsxfun(@rdivide,...
    pupilEMGdist.counts_low,sum(pupilEMGdist.counts_low,2));

% Saving to struct
PSSBehavior.pupilEMGdist = pupilEMGdist;

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
imagesc(pupilEMGdist.bins{1},pupilEMGdist.bins{2},pupilEMGdist.conditional_high')
hold on
imagesc(pupilEMGdist.bins{1}+2*pi,pupilEMGdist.bins{2},pupilEMGdist.conditional_high')
plot([-pi 3*pi],log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],'r--')
plot(cosx,cos(cosx)./3-1.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
xlabel('Pupil phase');ylabel('EMG')
title('sig pupil amp');

subplot(2,2,4);
imagesc(pupilEMGdist.bins{1},pupilEMGdist.bins{2},pupilEMGdist.conditional_low')
hold on
imagesc(pupilEMGdist.bins{1}+2*pi,pupilEMGdist.bins{2},pupilEMGdist.conditional_low')
plot([-pi 3*pi],log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],'r--')
plot(cosx,0.2.*cos(cosx)./3-1.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
xlabel('Pupil phase');ylabel('EMG')
title('ns pupil amp');

NiceSave('EMGbyPupilPhase',figfolder,baseName)

%% EMG-PSS codistribution:
PSSEMGdist.edges = {linspace(-2,1,50),linspace(-3,0,50)};
[PSSEMGdist.counts,PSSEMGdist.bins] = hist3([log10(PSS.EMG),PSS.data],...
    'Edges',PSSEMGdist.edges);
PSSEMGdist.joint = PSSEMGdist.counts./sum(PSSEMGdist.counts(:));

PSSEMGdist.conditional = bsxfun(@rdivide,...
    PSSEMGdist.counts,sum(PSSEMGdist.counts,2));

% Saving to struct
PSSBehavior.PSSEMGdist = PSSEMGdist;

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

NiceSave('EMGPSS',figfolder,baseName)

%% Conditional histos PSS/EMG/Pupil dynamics
[EMGPupPSS.meanZ,EMGPupPSS.N,EMGPupPSS.Xbins,...
    EMGPupPSS.Ybins ] = ConditionalHist3( PSS.pupphase,PSS.data,log10(PSS.EMG),...
    'minXY',0,'Xbounds',[-pi pi],'Ybounds',[-3 1],...
    'numXbins',100,'numYbins',100);

% PSS by lo/hi pupil instances
[PSSbyhiPup.meanZ,PSSbyhiPup.N,PSSbyhiPup.Xbins,...
    PSSbyhiPup.Ybins ] = ConditionalHist3( PSS.pupphase(PSS.highpup),...
    log10(PSS.EMG(PSS.highpup)),...
    PSS.data(PSS.highpup),'minXY',0,'Xbounds',[-pi pi],'Ybounds',[-3 1],...
    'numXbins',100,'numYbins',100);

[PSSbyloPup.meanZ,PSSbyloPup.N,PSSbyloPup.Xbins,...
    PSSbyloPup.Ybins ] = ConditionalHist3( PSS.pupphase(~PSS.highpup),...
    log10(PSS.EMG(~PSS.highpup)),...
    PSS.data(~PSS.highpup),'minXY',0,'Xbounds',[-pi pi],'Ybounds',[-3 1],...
    'numXbins',100,'numYbins',100);

% Saving to struct
PSSBehavior.EMGPupPSS = EMGPupPSS;
PSSBehavior.PSSbyhiPup = PSSbyhiPup;
PSSBehavior.PSSbyloPup = PSSbyloPup;

%% FIGURE:
figure;
subplot(2,2,1);
a = imagesc([EMGPupPSS.Xbins EMGPupPSS.Xbins+2*pi],...
    EMGPupPSS.Ybins,[EMGPupPSS.meanZ; EMGPupPSS.meanZ]');
colormap(gca,'jet')
alpha(a,double(~isnan([EMGPupPSS.meanZ; EMGPupPSS.meanZ]')))
%ColorbarWithAxis([min(min(EMGPupPSS.meanZ)) max(max(EMGPupPSS.meanZ))],'EMG')
%caxis([min(min(EMGPupPSS.meanZ)) max(max(EMGPupPSS.meanZ))])
ylim([-3 0]);
xlabel('Pupil phase');ylabel('PSS')
axis xy
title('EMG by PSS/Pupil');

subplot(2,2,2);
a = imagesc([PSSbyhiPup.Xbins PSSbyhiPup.Xbins+2*pi],...
    PSSbyhiPup.Ybins,[PSSbyhiPup.meanZ; PSSbyhiPup.meanZ]');
colormap(gca,'jet')
alpha(a,double(~isnan([PSSbyhiPup.meanZ; PSSbyhiPup.meanZ]')))
% ColorbarWithAxis([min(min(PSSbyhiPup.meanZ)) max(max(PSSbyhiPup.meanZ))],'PSS (au)')
% caxis([min(min(PSSbyhiPup.meanZ)) max(max(PSSbyhiPup.meanZ))])
ylim([-2 1.5]);
xlabel('Pupil phase');ylabel('EMG')
axis xy
title('sig pupil amp');

subplot(2,2,4);
a = imagesc([PSSbyloPup.Xbins PSSbyloPup.Xbins+2*pi],...
    PSSbyloPup.Ybins,[PSSbyloPup.meanZ; PSSbyloPup.meanZ]');
colormap(gca,'jet')
alpha(a,double(~isnan([PSSbyloPup.meanZ; PSSbyloPup.meanZ]')))
% ColorbarWithAxis([min(min(PSSbyloPup.meanZ)) max(max(PSSbyloPup.meanZ))],'PSS (au)')
% caxis([min(min(PSSbyloPup.meanZ)) max(max(PSSbyloPup.meanZ))])
ylim([-2 1.5]);
xlabel('Pupil phase');ylabel('EMG')
axis xy
title('ns pupil amp');

%NiceSave('PupilEMGPSS',figfolder,baseName)

%% Get PSS around Whisks
% Only LONG whisks selected...
EMGwhisk.dur = diff(EMGwhisk.ints.Wh,1,2);
% EMGwhisk.longwhisks = EMGwhisk.dur>1;

% For clean baseline whisk epochs
% tempidx = [];
% for i = 1:size(EMGwhisk.ints.Wh,1)
%     idx = find(EMGwhisk.ints.NWh(:,1) > EMGwhisk.ints.Wh(i,1)-5 & EMGwhisk.ints.NWh(:,2) < EMGwhisk.ints.Wh(i,1));
%     if ~isempty(idx)
%         tempidx = cat(1,tempidx,i);
%     else
%     end
% end

EMGwhisk.numwhisks = length(EMGwhisk.ints.Wh(:,1));
%EMGwhisk.numwhisks = length(EMGwhisk.ints.Wh(tempidx,1));

EMGwhisk.highpupil = EMGwhisk.power>nanmedian(EMGwhisk.power);
whiskPETH.window = [-5 5]; %s
whiskPETH.windex = 5*PSS.samplingRate; %s
whiskPETH.timestamps = whiskPETH.window(1):(1/PSS.samplingRate):whiskPETH.window(2);
timelockedPSS.data = zeros(length(whiskPETH.timestamps),EMGwhisk.numwhisks);

for ww = 1:EMGwhisk.numwhisks
    tempidx = interp1(PSS.timestamps,PSS.timestamps,EMGwhisk.ints.Wh(ww,1),'nearest');
    if ~isnan(tempidx)
        PSS.whidx(ww) = find(PSS.timestamps==tempidx);
        %PSS.whidx(ww) = find(PSS.timestamps==interp1(PSS.timestamps,PSS.timestamps,EMGwhisk.ints.Wh(tempidx(ww),1),'nearest'));
        
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
    else
%         PSS.whidx(ww) = [];
%         timelockedPSS.data(:,ww) = [];
%         timelockedPSS.timestamps(:,ww) = NaN(length(whiskPETH.timestamps),1);
%         timelockedPSS.phases(:,ww) = NaN(size(whiskPETH.timestamps));
%         timelockedPSS.highpupil(:,ww) = NaN(size(whiskPETH.timestamps));
    end
end
timelockedPSS.highpupil = logical(timelockedPSS.highpupil); %Why?

[phasePETH.high]=PairMatHist(timelockedPSS.data(timelockedPSS.highpupil),...%&~timelockedPSS.otherwhisks),...%
    [timelockedPSS.timestamps(timelockedPSS.highpupil),...%&~timelockedPSS.otherwhisks),...
    timelockedPSS.phases(timelockedPSS.highpupil)],...%&~timelockedPSS.otherwhisks)],...
    25,[-pi 4]);
%phasePETH.high.mean = (phasePETH.high.mean-nanmean(phasePETH.high.mean(8:11,:),1))./nanmean(phasePETH.high.std(8:11,:),1);

[phasePETH.low]=PairMatHist(timelockedPSS.data(~timelockedPSS.highpupil),...%&~timelockedPSS.otherwhisks),...
    [timelockedPSS.timestamps(~timelockedPSS.highpupil),...%&~timelockedPSS.otherwhisks),...
    timelockedPSS.phases(~timelockedPSS.highpupil)],...%&~timelockedPSS.otherwhisks)],...
    25,[-pi 5]);
%phasePETH.low.mean = (phasePETH.low.mean-nanmean(phasePETH.low.mean(8:11,:),1))./nanmean(phasePETH.low.std(8:11,:),1);

% Sorting Wh by phase/duration
[whisksorts.phaseval,whisksorts.phase] = sort(EMGwhisk.phase);
[whisksorts.durval,whisksorts.dur] = sort(EMGwhisk.dur);

% 
tempidx = interp1(PSS.timestamps,PSS.timestamps,EMGwhisk.ints.Wh,'nearest');
whints_maxpss = NaN(size(EMGwhisk.ints.Wh,1),1);
for i = 1:size(EMGwhisk.ints.Wh,1)
    temps = find(PSS.timestamps == tempidx(i,1));
    tempe = find(PSS.timestamps == tempidx(i,2));
    
    if ~isempty(temps)
        if ~isempty(tempe)
        whints_maxpss(i) = max(PSS.data(temps:tempe));
        else
        end
    else
    end
end

tempidx = interp1(EMGwhisk.timestamps,EMGwhisk.timestamps,EMGwhisk.ints.Wh,'nearest');
whints_amp = NaN(size(EMGwhisk.ints.Wh,1),1);
for i = 1:size(EMGwhisk.ints.Wh,1)
    temps = find(EMGwhisk.timestamps == tempidx(i,1));
    tempe = find(EMGwhisk.timestamps == tempidx(i,2));
    
    if ~isempty(temps)
        if ~isempty(tempe)
        whints_amp(i) = max(EMGwhisk.EMGenvelope(temps:tempe));
        else
        end
    else
    end
end

% Saving to struct
PSSBehavior.whiskPETH = whiskPETH;
PSSBehavior.timelockedPSS = timelockedPSS;
PSSBehavior.phasePETH = phasePETH;
PSSBehavior.whiskidx_phase = EMGwhisk.phase;
PSSBehavior.whiskidx_dur = EMGwhisk.dur;
PSSBehavior.whiskidx_maxpss = whints_maxpss;
PSSBehavior.whiskidx_amp = whints_amp;

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
% caxis([max(max(phasePETH.high.mean))*-1 max(max(phasePETH.high.mean))])
xlim([-1 4]);ylim([-pi 3*pi])
xlabel('t (s, aligned to Wh Onset)');ylabel('Pupil Phase')
title('sig pupil amp')

subplot(1,2,2);
imagesc(phasePETH.low.bincenters,phasePETH.low.bincenters,phasePETH.low.mean'); hold on
imagesc(phasePETH.low.bincenters,phasePETH.low.bincenters+2*pi,phasePETH.low.mean')
plot(EMGwhisk.dur(~EMGwhisk.highpupil),EMGwhisk.phase(~EMGwhisk.highpupil),'r.','markersize',1)
plot(EMGwhisk.dur(~EMGwhisk.highpupil),EMGwhisk.phase(~EMGwhisk.highpupil)+2*pi,'r.','markersize',1)
plot(0.2.*cos(cosx),cosx,'w','linewidth',2)
plot([0 0],[-pi 3*pi],'r')
colorbar
axis xy
% caxis([max(max(phasePETH.high.mean))*-1 max(max(phasePETH.high.mean))])
xlim([-1 4]);ylim([-pi 3*pi])
xlabel('t (s, aligned to Wh Onset)');ylabel('Pupil Phase')
title('ns pupil amp')

NiceSave('zPSS_PETHbyPhase',figfolder,baseName)

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

% sorted by whisking amplitude
% sorted by pupil amplitude

NiceSave('PSS_sortedWhisks',figfolder,baseName)

%% Get PSS around Pupil dilation, sorted by Whisking onset
Pupon.whon = NaN(length(pup_on),2);
for i = 1:length(pup_on)
    %Pupon.whon(i) = interp1(EMGwhisk.ints.Wh(:,1),EMGwhisk.ints.Wh(:,1),pup_on(i),'nearest');
    tempidx = find(EMGwhisk.ints.Wh(:,1)<pup_on(i),1,'last');
    if ~isnan(tempidx)
        Pupon.whon(i,:) = EMGwhisk.ints.Wh(tempidx,:);
    else
    end
end
Pupon.pupwhdiff = Pupon.whon(:,1) - pup_on;
%Pupon.pupwhend = Pupon.pupwhdiff + Pupon.whon(:,2)- Pupon.whon(:,1);

Pupon.numpups = length(pup_on);
whiskPETH.window = [-10 10]; %s
whiskPETH.windex = 10*PSS.samplingRate; %s
whiskPETH.timestamps = whiskPETH.window(1):(1/PSS.samplingRate):whiskPETH.window(2);
puplockedPSS.data = zeros(length(whiskPETH.timestamps),Pupon.numpups);

for ww = 1:Pupon.numpups
    tempidx = interp1(PSS.timestamps,PSS.timestamps,pup_on(ww),'nearest');
    if ~isnan(tempidx)
        PSS.pupidx(ww) = find(PSS.timestamps==tempidx);
        
        if PSS.pupidx(ww)-whiskPETH.windex > 0 && PSS.pupidx(ww)+whiskPETH.windex < length(PSS.data)
            puplockedPSS.data(:,ww) = PSS.data(PSS.pupidx(ww)-whiskPETH.windex:PSS.pupidx(ww)+whiskPETH.windex);
            puplockedPSS.timestamps(:,ww) = whiskPETH.timestamps;
            puplockedPSS.pupwhdiff(:,ww) = ones(size(whiskPETH.timestamps)).*Pupon.pupwhdiff(ww);
        end
    else
    end
end

% Sorting trials by params
[pupsorts.Pupwhdiffval,pupsorts.Pupwhdiff] = sort(Pupon.pupwhdiff);
[pupsorts.dPval,pupsorts.dP] = sort(dPpeak);
[pupsorts.durval,pupsorts.dur] = sort(Pupdur);
[pupsorts.peakval,pupsorts.peak] = sort(pup_peak);

% Saving to struct
PSSBehavior.puplockedPSS = puplockedPSS;
PSSBehavior.pupidx_pupwhdiff = Pupon.pupwhdiff;
PSSBehavior.pupidx_dur = Pupdur;
PSSBehavior.pupidx_dPpeak = dPpeak;
PSSBehavior.pupidx_amp = pup_peak;
PSSBehavior.pupidx_phase = Pupon.phase; 
PSSBehavior.pupidx_pow = Pupon.power;

save(savefile,'PSSBehavior');

%% FIGURE:
figure; 
subplot(2,2,1);
imagesc(puplockedPSS.timestamps(:,1),[1:Pupon.numpups],puplockedPSS.data(:,pupsorts.Pupwhdiff)'); hold on;
plot(Pupon.pupwhdiff(pupsorts.Pupwhdiff),[1:Pupon.numpups],'r.','markersize',5)
plot([0 0],[1 Pupon.numpups],'b')
axis square
xlim([-10 10])
%colorbar; caxis([-3 -1])
xlabel('t (s, aligned to Pup Onset)'); ylabel('trial no.')
title('Epochs sorted by Wh-Pup onset diff');

subplot(2,2,2);
imagesc(puplockedPSS.timestamps(:,1),[1:Pupon.numpups],puplockedPSS.data(:,pupsorts.dP)'); hold on;
%plot(Pupon.pupwhdiff(pupsorts.dP),[1:Pupon.numpups],'r.','markersize',5)
plot([0 0],[1 Pupon.numpups],'b')
axis square
xlim([-10 10])
%colorbar; caxis([-3 -1])
xlabel('t (s, aligned to Pup Onset)'); ylabel('trial no.')
title('Epochs sorted by dP/dt');

subplot(2,2,3);
imagesc(puplockedPSS.timestamps(:,1),[1:Pupon.numpups],puplockedPSS.data(:,pupsorts.peak)'); hold on;
%plot(Pupon.pupwhdiff(pupsorts.dP),[1:Pupon.numpups],'r.','markersize',5)
plot([0 0],[1 Pupon.numpups],'b')
axis square
xlim([-10 10])
%colorbar; caxis([-3 -1])
xlabel('t (s, aligned to Pup Onset)'); ylabel('trial no.')
title('Epochs sorted by peak pupil area');

subplot(2,2,4);
imagesc(puplockedPSS.timestamps(:,1),[1:Pupon.numpups],puplockedPSS.data(:,pupsorts.dur)'); hold on;
%plot(pup_off(pupsorts.dur)-pup_on(pupsorts.dur),[1:Pupon.numpups],'r.','markersize',5)
plot([0 0],[1 Pupon.numpups],'b')
axis square
xlim([-10 10])
%colorbar; caxis([-3 -1])
xlabel('t (s, aligned to Pup Onset)'); ylabel('trial no.')
title('Epochs sorted by epoch duration');

NiceSave('PSS_sortedPupil',figfolder,baseName)