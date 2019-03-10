%%
basePath = pwd;
baseName = bz_BasenameFromBasepath(basePath);
sessionInfo = bz_getSessionInfo(basePath,'noPrompts',true);

badchannels = sessionInfo.badchannels;
usechannels = sessionInfo.AnatGrps.Channels;
usechannels(ismember(usechannels,badchannels))=[];
channels = sessionInfo.channels;

figfolder = fullfile(basePath,'AnalysisFigures');
savefile = fullfile(basePath,[baseName,'.BehaviorAnalysis.mat']);

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

%% PUPIL STATS
% Computing x-covariance
acgwin = 200; %s
[pupACG.ACG,pupACG.tlag] = xcov(pupildilation.data,...
    round(acgwin.*pupildilation.samplingRate),'coeff');
pupACG.tlag = pupACG.tlag./pupildilation.samplingRate;

% Pupil diameter histogram
puphist.bins = linspace(0,5,20);
puphist.counts = hist(pupildilation.data,puphist.bins);
puphist.counts =puphist.counts./sum(puphist.counts);

% dPupil/dt histogram
pupdthist.bins={linspace(-0.3,0.3,25),linspace(-0.3,0.3,25)};
pupdthist.counts = hist3([log10(pupildilation.data(1:end-1)),pupildilation.dpdt],pupdthist.bins);
pupdthist.counts = pupdthist.counts./sum(pupdthist.counts(:));

% Pupil PSD
frange = [0.001 10];
pupilspec = bz_WaveSpec(pupildilation.data,...
   'frange',frange,'nfreqs',100,'ncyc',3,...
   'samplingRate',pupildilation.samplingRate);
freqs = pupilspec.freqs;
spec = pupilspec.data;
spec = (abs(spec));
pupPSD.freqs = freqs;
pupPSD.psd = mean(log10(spec),1);
clear pupilspec spec

% Filtered Pupil
lowfilter = [0.01 0.1];
%highfilter = [0.3 0.8];

pupil4filter = pupildilation;
lowpupildata = bz_Filter(pupil4filter,'passband',lowfilter,'filter' ,'fir1','order',3);
%highpupildata = bz_Filter(pupil4filter,'passband',highfilter,'filter' ,'fir1');

% Saving to struct
PupEMG.pupACG = pupACG;
PupEMG.puphist = puphist;
PupEMG.pupdthist = pupdthist;
PupEMG.pupPSD = pupPSD;

%% FIGURE 1:
figure;
subplot(4,2,1:2);
plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2); hold on;
scatter(lowpupildata.timestamps,lowpupildata.data+0.5.*nanmean(pupildilation.data),4,lowpupildata.phase)
%plot(highpupildata.timestamps,highpupildata.data+nanmean(pupildilation.data),'r')
colormap(gca,'jet'); 
ColorbarWithAxis([min(lowpupildata.phase) max(lowpupildata.phase)],['pupil phase'])
h1 = plot(get(gca,'xlim'),[0 0],'k-');
plot(pupildilation.timestamps(1:end-1),pupildilation.dpdt,'r','linewidth',2)
xlim([100 350]); ylim([-0.5 2.5]);
xlabel('time (s)'); ylabel('Pupil');
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'diameter','phase','dPdt'},'location','northeast');

subplot(4,2,3);
bar(puphist.bins,puphist.counts,'facecolor','k')
axis tight;
xlim([0 3]);
xlabel('diameter (norm.)'); ylabel('counts (au)');

subplot(4,2,5);
plot(pupACG.tlag,pupACG.ACG,'k','linewidth',2); hold on;
axis tight
plot(get(gca,'xlim'),[0 0],'r')
xlabel('lag (s)')
ylabel('autocovariance')
legend({'diameter'},'location','northeast');

subplot(4,2,7);
plot(log10(pupPSD.freqs),pupPSD.psd,'k','linewidth',2)
hold on
plot(log10(lowfilter),[-1 -1],'r')
LogScale('x',10)
xlabel('f (Hz)'); ylabel('power (dB)')
axis tight
ylim([-2 max(pupPSD.psd)]);

subplot(4,2,4);
imagesc(pupdthist.bins{1},pupdthist.bins{2},pupdthist.counts'); hold on;
colormap(gca,[1 1 1; colormap('jet')])
plot(get(gca,'xlim'),[0 0],'r-')
ColorbarWithAxis([min(min(pupdthist.counts)) max(max(pupdthist.counts))],['counts (au)'])
xlabel('diameter (norm)'); ylabel('dP/dt')
LogScale('x',10);
axis square
ylim([min(pupildilation.dpdt) max(pupildilation.dpdt)]);

subplot(4,2,6)
colormap(gca,'jet')
scatter(log10(pupildilation.data(1:end-1)),pupildilation.dpdt,0.2,lowpupildata.phase(1:end-1))
LogScale('x',10)
ColorbarWithAxis([-pi pi],'phase')
xlabel('diameter');ylabel('dP/dt')
axis square; box on;
ylim([-0.7 0.7]);xlim([-0.5 0.5])

subplot(4,2,8);
%colormap(gca,distcolor)
scatter(log10(pupildilation.data(1:end-1)),pupildilation.dpdt,0.2,lowpupildata.amp(1:end-1))
%plot(pupildilation.data(1:end-1),pupildilation.dpdt,'k.','markersize',2)
LogScale('x',10)
ColorbarWithAxis([0 1.2],'power')
xlabel('diameter');ylabel('dP/dt')
axis square; box on;
ylim([-0.7 0.7]);xlim([-0.5 0.5])

NiceSave('PupilStats',figfolder,baseName)

%% EMG STATS
% EMG envelope histogram
EMGhist.bins = linspace(0,20,20);
EMGhist.counts = hist(EMGwhisk.EMGsm,EMGhist.bins);
EMGhist.counts =EMGhist.counts./sum(EMGhist.counts);

EMGhist.logbins = linspace(-1.5,1.5,20);
EMGhist.logcounts = hist(log10(EMGwhisk.EMGsm),EMGhist.logbins);
EMGhist.logcounts =EMGhist.logcounts./sum(EMGhist.logcounts);

EMGhist.threshold = EMGwhisk.detectorparms.Whthreshold;

% Whisk Durations
EMGwhisk.Whdurs = diff(EMGwhisk.ints.Wh,1,2);
EMGwhisk.InterWhdurs = EMGwhisk.ints.Wh(2:end,1)-EMGwhisk.ints.Wh(1:end-1,2);

Whdurhist.bins = linspace(-1,2,20);
Whdurhist.Whdurs = hist(log10(EMGwhisk.Whdurs),Whdurhist.bins);
Whdurhist.Whdurs = Whdurhist.Whdurs./diff(EMGwhisk.timestamps([1 end])); %Units: whisks/s
Whdurhist.InterWhdurs = hist(log10(EMGwhisk.InterWhdurs),Whdurhist.bins);
Whdurhist.InterWhdurs = Whdurhist.InterWhdurs./diff(EMGwhisk.timestamps([1 end])); %Units: whisks/s

% EMG PSD
frange = [0.001 10];
EMGspec = bz_WaveSpec(EMGwhisk.EMGsm,...
   'frange',frange,'nfreqs',100,'ncyc',3,...
   'samplingRate',EMGwhisk.samplingRate);
freqs = EMGspec.freqs;
spec = EMGspec.data;
spec = (abs(spec));
EMGPSD.freqs = freqs;
EMGPSD.psd = mean(log10(spec),1);
clear EMGspec spec

% Saving to struct
PupEMG.EMGhist = EMGhist;
PupEMG.Whdurhist = Whdurhist;
PupEMG.EMGPSD = EMGPSD;

%% FIGURE 2:
figure;
subplot(3,1,1);
plot(EMGwhisk.timestamps,EMGwhisk.EMG,'color',[0.5 0.5 0.5],'linewidth',0.5); hold on;
plot(pupildilation.timestamps,EMGwhisk.pupiltime,'b','linewidth',2);
plot(EMGwhisk.ints.Wh',...
    ones(size(EMGwhisk.ints.Wh))',...
    'g-','linewidth',1)
ylabel('EMG'); xlabel('time (s)');
axis tight
xlim([100 250]); ylim([-10 40]);

subplot(3,2,3);
bar(EMGhist.logbins,EMGhist.logcounts,'facecolor','b'); hold on
plot(log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],get(gca,'ylim'),'r--','linewidth',2)
LogScale('x',10);
axis tight
xlabel('EMG (norm.)');
legend({'distribution','Wh threshold'},'location','northeast');

subplot(3,2,4);
plot(Whdurhist.bins,Whdurhist.Whdurs,'b','linewidth',2)
hold on
plot(Whdurhist.bins,Whdurhist.InterWhdurs,'k','linewidth',2)
LogScale('x',10)
legend('Wh','NWh')
xlabel('duration (s)');

subplot(3,2,5);
plot(log10(EMGPSD.freqs),EMGPSD.psd,'k','linewidth',2)
hold on
LogScale('x',10)
xlabel('f (Hz)'); ylabel('power (dB)')
axis tight
ylim([min(EMGPSD.psd) max(EMGPSD.psd)]);

NiceSave('EMGStats',figfolder,baseName)

%% EMG/PUPIL correlation
% Codistribution
pupilEMGdist.bins = {linspace(-0.5,0.5,70),linspace(-1.5,1,70)};
[pupilEMGdist.counts,pupilEMGdist.bins] = hist3([log10(pupildilation.data),log10(EMGwhisk.pupiltime)],pupilEMGdist.bins);
pupilEMGdist.counts = pupilEMGdist.counts./sum(pupilEMGdist.counts(:));

% XCovariance
[pupilEMGcorr.xcorr,pupilEMGcorr.corrlags] = xcov(pupildilation.data,EMGwhisk.pupiltime,'unbiased');
pupilEMGcorr.corrlags = pupilEMGcorr.corrlags.*(1./pupildilation.samplingRate);

% Pupil/EMG at whisking onset
[pwCCG.pupil.WhOn,t_lag,~,alltrans.pupil.WhOn,skippedWh] = EventVsContinousCCG(pupildilation.data,pupildilation.timestamps,EMGwhisk.ints.Wh(:,1),20);
[pwCCG.pupilxy.WhOn(:,1),t_lag,~,~] = EventVsContinousCCG(pupildilation.pupilxy(:,1),pupildilation.timestamps,EMGwhisk.ints.Wh(:,1),20);
[pwCCG.pupilxy.WhOn(:,2),t_lag,~,~] = EventVsContinousCCG(pupildilation.pupilxy(:,2),pupildilation.timestamps,EMGwhisk.ints.Wh(:,1),20);
[pwCCG.EMG.WhOn,t_lag,~,alltrans.EMG.WhOn] = EventVsContinousCCG(EMGwhisk.pupiltime,pupildilation.timestamps,EMGwhisk.ints.Wh(:,1),20);
pwCCG.t_lag = t_lag;

% Histogram: Whisking by Pupil Dynamics
[pupildynamicsEMG,pupildynamicsEMG.bins]=PairMatHist(log10(EMGwhisk.pupiltime(1:end-1)),...
    [log10(pupildilation.data(1:end-1)),pupildilation.dpdt],15,[-0.5 0.5]);

[~,pupildilation.iswhisk] = RestrictInts(pupildilation.timestamps,EMGwhisk.ints.Wh);
[pWhisk]=PairMatHist(single(pupildilation.iswhisk(1:end-1)),[log10(pupildilation.data(1:end-1)),pupildilation.dpdt],pupildynamicsEMG.binedges);
pupildynamicsEMG.pWhisk = pWhisk.mean;

[pupilphaseEMG,pupilphaseEMG.bins]=PairMatHist(log10(EMGwhisk.pupiltime),...
    [log10(lowpupildata.amp),lowpupildata.phase],20,[-pi pi]);
[pWhisk]=PairMatHist(single(pupildilation.iswhisk),[log10(lowpupildata.amp),lowpupildata.phase],pupilphaseEMG.binedges);
pupilphaseEMG.pWhisk = pWhisk.mean;

% Wh onsets/offsets in pupil space
whints_pupil = interp1(pupildilation.timestamps,pupildilation.data,EMGwhisk.ints.Wh);
whints_pupildt = interp1(pupildilation.timestamps(1:end-1),pupildilation.dpdt,EMGwhisk.ints.Wh);
whints_pupilphase = interp1(lowpupildata.timestamps,lowpupildata.phase,EMGwhisk.ints.Wh);
whints_pupilamp = interp1(lowpupildata.timestamps,lowpupildata.amp,EMGwhisk.ints.Wh);

[pWhiskStart]=PairMatHist(EMGwhisk.Whdurs(:,1),[log10(whints_pupil(:,1)),whints_pupildt(:,1)],pupildynamicsEMG.binedges);
pupildynamicsEMG.meanWhdur = pWhiskStart.mean;
pupildynamicsEMG.numWhstarts = pWhiskStart.num;
pupildynamicsEMG.pWhstarts = pupildynamicsEMG.numWhstarts./sum(pupildynamicsEMG.numWhstarts(:));
pupildynamicsEMG.occupancy = pupildynamicsEMG.num./pupildilation.samplingRate;
pupildynamicsEMG.Whstartrate = pupildynamicsEMG.numWhstarts./pupildynamicsEMG.occupancy;

[pWhiskStart]=PairMatHist(EMGwhisk.Whdurs(:,1),[log10(whints_pupilamp(:,1)),whints_pupilphase(:,1)],pupilphaseEMG.binedges);
pupilphaseEMG.meanWhdur = pWhiskStart.mean;
pupilphaseEMG.numWhstarts = pWhiskStart.num;
pupilphaseEMG.pWhstarts = pupilphaseEMG.numWhstarts./sum(pupilphaseEMG.numWhstarts(:));
pupilphaseEMG.occupancy = pupilphaseEMG.num./pupildilation.samplingRate;
pupilphaseEMG.Whstartrate = pupilphaseEMG.numWhstarts./pupilphaseEMG.occupancy;

occupancythresh = 2; %s (don't count bins with less than 1s total time)
whthresh = 2;
pupildynamicsEMG.mean(pupildynamicsEMG.occupancy<=occupancythresh)=nan;
pupilphaseEMG.mean(pupilphaseEMG.occupancy<=occupancythresh)=nan;
pupildynamicsEMG.Whstartrate(pupildynamicsEMG.occupancy<=occupancythresh)=nan;
pupilphaseEMG.Whstartrate(pupilphaseEMG.occupancy<=occupancythresh)=nan;

pupildynamicsEMG.meanWhdur(pupildynamicsEMG.numWhstarts<=whthresh)=nan;
pupilphaseEMG.meanWhdur(pupilphaseEMG.numWhstarts<=whthresh)=nan;

% Saving to struct
PupEMG.pupilEMGdist = pupilEMGdist;
PupEMG.pupilEMGcorr = pupilEMGcorr;
PupEMG.pwCCG = pwCCG;
PupEMG.pWhiskStart = pWhiskStart;
PupEMG.pWhisk = pWhisk;
PupEMG.pupildynamicsEMG = pupildynamicsEMG; 
PupEMG.pupilphaseEMG = pupilphaseEMG;
PupEMG.whints_pupil = whints_pupil;
PupEMG.whints_pupildt = whints_pupildt;
PupEMG.whints_pupilphase = whints_pupilphase;
PupEMG.whints_pupilamp = whints_pupilamp;

%% FIGURE 3:
distcolor = [1 1 1; makeColorMap([0.7 0.7 0.7],[0 0.5 0],[0.7 0.6,0])];
emgcolor = [1 1 1;makeColorMap([0.5 0.5 0.5],[0 0 0.8])];

figure;
subplot(3,2,1:2)
plot(pupildilation.timestamps,...
    pupildilation.data./max(pupildilation.data),'k','linewidth',2); hold on;
h1 = plot(EMGwhisk.timestamps,EMGwhisk.EMG./max(EMGwhisk.EMG),...
    'color',[0.5 0.5 0.5],'linewidth',0.5);
hold on
plot(pupildilation.timestamps,EMGwhisk.pupiltime./max(EMGwhisk.pupiltime),'b')
plot(EMGwhisk.ints.Wh',...
    zeros(size(EMGwhisk.ints.Wh))',...
    'g-','linewidth',1)
axis tight
xlim([100 250]); ylim([0 1.25])
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'Pupil diameter','EMG','Wh'},'location','northeast');

subplot(3,2,3);
colormap(gca,[1 1 1;colormap]);
imagesc(pupilEMGdist.bins{1},pupilEMGdist.bins{2},...
    (pupilEMGdist.counts)'./max(pupilEMGdist.counts(:))); hold on
plot(get(gca,'xlim'),log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],'g--')
LogScale('y',10); axis xy
ColorbarWithAxis([0 0.6],['counts (au)'])
xlabel('Pupil diameter (norm.)'); ylabel('EMG');

subplot(3,2,5);
scatter(log10(pupildilation.data(1:end-1)),pupildilation.dpdt,0.2,log10(EMGwhisk.pupiltime(1:end-1)))
ylim([-0.5 0.5]);xlim([-0.5 0.5])
LogScale('x',10)
ColorbarWithAxis([-1.5 1.8],'EMG Envelope')
xlabel('Pupil Diameter');ylabel('dP/dt')

subplot(3,2,4);
scatter(lowpupildata.phase,log10(lowpupildata.amp),0.2,log10(EMGwhisk.pupiltime))
hold on
scatter(lowpupildata.phase+2*pi,log10(lowpupildata.amp),0.2,log10(EMGwhisk.pupiltime))
axis tight
caxis([-1.5 1.8])
xlabel('Phase');ylabel('Amplitude')

subplot(3,2,6);
colormap(gca,emgcolor(2:end,:))
scatter(lowpupildata.phase,log10(lowpupildata.amp),0.2,pupildilation.iswhisk)
hold on
scatter(lowpupildata.phase+2*pi,log10(lowpupildata.amp),0.2,pupildilation.iswhisk)
axis tight
xlabel('Phase');ylabel('Amplitude')

NiceSave('EMG_Pupil_Amp_Corr',figfolder,baseName)

%% FIGURE 4:
emgcolor = [1 1 1;makeColorMap([0.5 0.5 0.5],[0 0 0.8])];

figure;
subplot(4,2,1);
imagesc(pupildynamicsEMG.bins,pupildynamicsEMG.bins,(pupildynamicsEMG.mean)');
colormap(gca,emgcolor)
axis xy
xlim([-0.5 0.5]);ylim([-0.3 0.3])
colorbar
ylabel('dP/dt')
title('EMG')
hold on
caxis([-0.75 0.5])
LogScale('c',10)
plot(log10(whints_pupil(:,1)),whints_pupildt(:,1),'k.','markersize',0.05)
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color','r')

subplot(4,2,3);
imagesc(pupildynamicsEMG.bins,pupildynamicsEMG.bins,pupildynamicsEMG.pWhisk');
colormap(gca,emgcolor)
axis xy
xlim([-0.5 0.5]);ylim([-0.3 0.3])
caxis([-0.01 1])
colorbar
ylabel('dP/dt')
title('p(Whisking)')
hold on
plot(log10(whints_pupil(:,1)),whints_pupildt(:,1),'k.','markersize',0.05)
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color','r')

subplot(4,2,5);
imagesc(pupildynamicsEMG.bins,pupildynamicsEMG.bins,(pupildynamicsEMG.Whstartrate)');
colormap(gca,emgcolor)
axis xy
xlim([-0.5 0.5]);ylim([-0.3 0.3])
colorbar
ylabel('dP/dt')
title('p(Wh onset)')
hold on
caxis([-0.1 1])
plot(log10(whints_pupil(:,1)),whints_pupildt(:,1),'k.','markersize',0.05)
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color','r')

subplot(4,2,7);
imagesc(pupildynamicsEMG.bins,pupildynamicsEMG.bins,log10(pupildynamicsEMG.meanWhdur)');
axis xy
xlim([-0.5 0.5]);ylim([-0.3 0.3])
colorbar
xlabel('Pupil Diameter');ylabel('dP/dt')
title('Whisk Duration')
hold on
caxis([-1 0.5])
LogScale('c',10)
plot(log10(whints_pupil(:,1)),whints_pupildt(:,1),'k.','markersize',0.05)
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color','r')

subplot(4,2,2);
imagesc(pupilphaseEMG.bins,pupilphaseEMG.bins,(pupilphaseEMG.mean));
hold on
imagesc(pupilphaseEMG.bins+2*pi,pupilphaseEMG.bins,(pupilphaseEMG.mean));
colormap(gca,emgcolor)
axis xy
ylim([-2 0]);xlim([-pi 3*pi])
colorbar
ylabel('Power')
title('EMG')
hold on
caxis([-0.75 0.5])
LogScale('c',10)
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color',0.5.*[1 1 1])

subplot(4,2,4);
imagesc(pupilphaseEMG.bins,pupilphaseEMG.bins,(pupilphaseEMG.pWhisk));
hold on
imagesc(pupilphaseEMG.bins+2*pi,pupilphaseEMG.bins,(pupilphaseEMG.pWhisk));
colormap(gca,emgcolor)
axis xy
ylim([-2 0]);xlim([-pi 3*pi])
colorbar
ylabel('Power')
title('p(Whisking)')
hold on
caxis([-0.01 1])
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color',0.5.*[1 1 1])

subplot(4,2,6);
imagesc(pupilphaseEMG.bins,pupilphaseEMG.bins,(pupilphaseEMG.Whstartrate));
hold on
imagesc(pupilphaseEMG.bins+2*pi,pupilphaseEMG.bins,(pupilphaseEMG.Whstartrate));
colormap(gca,emgcolor)
axis xy
ylim([-2 0]);xlim([-pi 3*pi])
colorbar
ylabel('Power')
title('p(Wh onset)')
hold on
caxis([-0.1 1])
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color',0.5.*[1 1 1])

subplot(4,2,8);
imagesc(pupilphaseEMG.bins,pupilphaseEMG.bins,log10(pupilphaseEMG.meanWhdur));
hold on
imagesc(pupilphaseEMG.bins+2*pi,pupilphaseEMG.bins,log10(pupilphaseEMG.meanWhdur));
axis xy
ylim([-2 0]);xlim([-pi 3*pi])
colorbar
xlabel('Pupil phase');ylabel('Power')
title('Whisk duration')
hold on
caxis([-1 0.5])
LogScale('c',10)
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color',0.5.*[1 1 1])

NiceSave('EMG_Pupil_Space',figfolder,baseName);

%% EMG-Pupil Phase coupling
frange = [0.001 1];
[wavespec] = bz_WaveSpec(pupildilation.data,...
    'frange',frange,'nfreqs',100,'ncyc',4,'samplingRate',pupildilation.samplingRate);
wavespec.timestamps=pupildilation.timestamps;

whpow = EMGwhisk.pupiltime./mean(EMGwhisk.pupiltime);
powbins = 10;
coupling = zeros(wavespec.nfreqs,1);
for ff = 1:wavespec.nfreqs
    coupling(ff) = abs(mean(whpow.*exp(1i.*angle(wavespec.data(:,ff)))));
end
WPcoupling.coupling = coupling;
WPcoupling.freqs = wavespec.freqs;

% Saving to struct
PupEMG.WPcoupling = WPcoupling;
PupEMG.pupildilation = pupildilation;
PupEMG.EMGwhisk = EMGwhisk;

save(savefile,'PupEMG');

%% FIGURE 5:
lagwin = [-5 5];

figure;
subplot(3,2,1);
plot(pupilEMGcorr.corrlags,pupilEMGcorr.xcorr','k','linewidth',2); hold on
axis tight
xlim([-40 40])
plot([0 0],get(gca,'ylim'),'g-','linewidth',2)
xlabel('t lag (s)');ylabel('EMG')
title('Pupil-EMG CCG')

subplot(3,2,3);
plot(t_lag,pwCCG.pupil.WhOn,'k','linewidth',2); hold on; 
axis tight
plot([0 0],get(gca,'ylim'),'g-','linewidth',2)
xlim(lagwin)
ylabel('Pupil')
set(gca,'xticklabel',[])

subplot(3,2,5);
plot(t_lag,pwCCG.EMG.WhOn,'b','linewidth',2); hold on;
axis tight
plot([0 0],get(gca,'ylim'),'g-','linewidth',2)
xlim(lagwin)
xlabel('t (s - aligned to WhOn)')
ylabel('EMG')

subplot(3,2,2);
plot(log10(wavespec.freqs),coupling,'r','LineWidth',2)
LogScale('x',10);
xlabel('f (Hz)');
axis tight
title('Pupil-EMG phase coupling')

subplot(3,2,4);
scatter(log10(whints_pupil(:,1)),whints_pupildt(:,1),3,log10(EMGwhisk.Whdurs))
xlabel('Pupil diameter (norm.)');ylabel('dP/dt');
colormap(gca,'jet'); 
ColorbarWithAxis([min(min(log10(EMGwhisk.Whdurs)))...
    max(max(log10(EMGwhisk.Whdurs)))],['Wh duration (s)'])
LogScale('x',10)
axis tight

subplot(3,2,6);
scatter(whints_pupilphase(:,1),log10(whints_pupilamp(:,1)),3,log10(EMGwhisk.Whdurs))
hold on
scatter(whints_pupilphase(:,1)+2*pi,log10(whints_pupilamp(:,1)),1,log10(EMGwhisk.Whdurs))
colormap(gca,'jet'); 
ColorbarWithAxis([min(min(log10(EMGwhisk.Whdurs)))...
    max(max(log10(EMGwhisk.Whdurs)))],['Wh duration (s)'])
xlabel('Pupil phase');ylabel('Pupil amplitude')
LogScale('y',10)
axis tight

NiceSave('Pupil_EMG_Temp_Phase_corr',figfolder,baseName)