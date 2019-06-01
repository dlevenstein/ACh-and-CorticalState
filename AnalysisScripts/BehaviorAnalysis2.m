function [ ] = BehaviorAnalysis2(basePath,figfolder)

%Initiate Paths
reporoot = '/home/dlevenstein/ProjectRepos/ACh-and-CorticalState/';
%reporoot = '/Users/dlevenstein/Project Repos/ACh-and-CorticalState/';
basePath = '/mnt/proraidDL/Database/WMData/AChPupil/171209_WT_EM1M3/';
%basePath = '/mnt/proraidDL/Database/WMData/AChPupil/180706_WT_EM1M3/';
%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];


%%
baseName = bz_BasenameFromBasepath(basePath);
sessionInfo = bz_getSessionInfo(basePath,'noPrompts',true);

badchannels = sessionInfo.badchannels;
usechannels = sessionInfo.AnatGrps.Channels;
usechannels(ismember(usechannels,badchannels))=[];
channels = sessionInfo.channels;

%figfolder = fullfile(basePath,'AnalysisFigures');
%savefile = fullfile(basePath,[baseName,'.BehaviorAnalysis.mat']);

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

% Filtered Pupil
lowfilter = [0.01 0.1];
%highfilter = [0.3 0.8];

pupil4filter = pupildilation;
pupilcycle = bz_Filter(pupil4filter,'passband',lowfilter,'filter' ,'fir1','order',3);
%highpupildata = bz_Filter(pupil4filter,'passband',highfilter,'filter' ,'fir1');
pupilcycle.pupthresh = -0.8;
pupilcycle.highpup = log10(pupilcycle.amp)>pupilcycle.pupthresh; 

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

%% Get rid of recording start/stop artifact

maxtimejump = 1; %s
pupilcycle.amp = NanPadJumps( pupilcycle.timestamps,pupilcycle.amp,maxtimejump );
pupilcycle.phase = NanPadJumps( pupilcycle.timestamps,pupilcycle.phase,maxtimejump );
pupildilation.data = NanPadJumps( pupildilation.timestamps,pupildilation.data,maxtimejump );
EMGwhisk.EMGsm = NanPadJumps( EMGwhisk.timestamps,EMGwhisk.EMGsm,maxtimejump );





%% Pupil phase/amp and duration at whisks

EMGwhisk.whisks.pupphase = interp1(pupilcycle.timestamps,pupilcycle.phase,EMGwhisk.ints.Wh(:,1),'nearest');
EMGwhisk.whisks.pupamp = interp1(pupilcycle.timestamps,pupilcycle.amp,EMGwhisk.ints.Wh(:,1),'nearest');
EMGwhisk.whisks.pup = interp1(pupildilation.timestamps,pupildilation.data,EMGwhisk.ints.Wh(:,1),'nearest');
EMGwhisk.whisks.dur = diff(EMGwhisk.ints.Wh,[],2);
EMGwhisk.whisks.hipup = log10(EMGwhisk.whisks.pupamp)>pupilcycle.pupthresh;
EMGwhisk.whisks.lopup = ~EMGwhisk.whisks.hipup;

EMGwhisk.pupphase = interp1(pupilcycle.timestamps,pupilcycle.phase,EMGwhisk.timestamps,'nearest');
EMGwhisk.pupamp = interp1(pupilcycle.timestamps,pupilcycle.amp,EMGwhisk.timestamps,'nearest');
EMGwhisk.pup = interp1(pupilcycle.timestamps,pupilcycle.data,EMGwhisk.timestamps,'nearest');
EMGwhisk.hipup =  log10(EMGwhisk.pupamp)>pupilcycle.pupthresh;
EMGwhisk.lopup = ~EMGwhisk.hipup;

HILO = {'lopup','hipup'};
%%
EMGrange = [-2 1.2];
for pp= 1:2
    EMGdist.(HILO{pp})  = ConditionalHist(...
        EMGwhisk.pupphase(EMGwhisk.(HILO{pp})),...
        log10(EMGwhisk.EMGsm(EMGwhisk.(HILO{pp}))),...
        'Ybounds',EMGrange,'numYbins',100,'numXbins',20,'Xbounds',[-pi pi]);

end

EMGdist.pup  = ConditionalHist( EMGwhisk.pup,...
        log10(EMGwhisk.EMGsm),...
        'numXbins',20,'Xbounds',[-0.4 0.4],'Ybounds',EMGrange,'numYbins',100);

%%
cosx = linspace(-pi,pi,100);
cospamp = [0.08 0.8];
figure
subplot(2,2,1)
        for pp = 1:2
        imagesc( EMGdist.(HILO{pp}).Xbins+2*pi*(pp-1),...
            EMGdist.(HILO{pp}).Ybins,...
            EMGdist.(HILO{pp}).pYX')
        hold on; axis xy; box off
        plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp)-2,'k')
        end   
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        xlim([-pi 3*pi])
        xlabel('Pupil Phase');ylabel('EMG')
        crameri bilbao
        
subplot(2,2,2)

        imagesc( EMGdist.pup.Xbins,...
            EMGdist.pup.Ybins,...
            EMGdist.pup.pYX')
        hold on; axis xy; box off
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        xlabel('Pupil Size');ylabel('EMG')
   crameri bilbao







%% PUPIL STATS
% Pupil diameter histogram
puphist.bins = linspace(0,5,20);
puphist.counts = hist(pupildilation.data,puphist.bins);
puphist.counts = puphist.counts./sum(puphist.counts);

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

% Pupil dynamics conditional histogram
[phasedynamics.meanZ,phasedynamics.N,phasedynamics.Xbins,...
    phasedynamics.Ybins ] = ConditionalHist3( log10(pupildilation.data(1:end-1)),...
    pupildilation.dpdt,lowpupildata.phase(1:end-1),...
    'minXY',0,'Xbounds',[-0.5 0.5],'Ybounds',[-0.5 0.5],...
    'numXbins',100,'numYbins',100,'circstats',true);

[ampdynamics.meanZ,ampdynamics.N,ampdynamics.Xbins,...
    ampdynamics.Ybins ] = ConditionalHist3( log10(pupildilation.data(1:end-1)),...
    pupildilation.dpdt,lowpupildata.amp(1:end-1),...
    'minXY',0,'Xbounds',[-0.5 0.5],'Ybounds',[-0.5 0.5],...
    'numXbins',100,'numYbins',100);

[areadynamics.meanZ,areadynamics.N,areadynamics.Xbins,...
    areadynamics.Ybins] = ConditionalHist3( lowpupildata.phase,...
    log10(lowpupildata.amp),log10(pupildilation.data),...
    'minXY',0,'Xbounds',[-pi pi],'Ybounds',[-3 1],...
    'numXbins',100,'numYbins',100);

[dpdtdynamics.meanZ,dpdtdynamics.N,dpdtdynamics.Xbins,...
    dpdtdynamics.Ybins] = ConditionalHist3( lowpupildata.phase(1:end-1),...
    log10(lowpupildata.amp(1:end-1)),pupildilation.dpdt,...
    'minXY',0,'Xbounds',[-pi pi],'Ybounds',[-3 1],...
    'numXbins',100,'numYbins',100);

% Saving to struct
PupEMG.puphist = puphist;
PupEMG.pupdthist = pupdthist;
PupEMG.pupPSD = pupPSD;
PupEMG.phasedynamics = phasedynamics;
PupEMG.ampdynamics = ampdynamics;
PupEMG.areadynamics = areadynamics;
PupEMG.dpdtdynamics = dpdtdynamics;

%% FIGURE 1:
figure;
subplot(4,2,1:2);
plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2); hold on;
scatter(lowpupildata.timestamps,lowpupildata.data+0.5.*nanmean(pupildilation.data),4,lowpupildata.phase)
%plot(highpupildata.timestamps,highpupildata.data+nanmean(pupildilation.data),'r')
colormap(gca,hsv)
ColorbarWithAxis([min(lowpupildata.phase) max(lowpupildata.phase)],['pupil phase'])
h1 = plot(get(gca,'xlim'),[0 0],'k-');
plot(pupildilation.timestamps(1:end-1),pupildilation.dpdt,'r','linewidth',2)
xlim([100 350]); ylim([-0.5 2.5]);
xlabel('time (s)'); ylabel('Pupil');
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'diameter','phase','dPdt'},'location','northeast');

subplot(4,3,4);
bar(puphist.bins,puphist.counts,'facecolor','k')
axis tight;
xlim([0 3]);
xlabel('diameter (norm.)'); ylabel('counts (au)');

subplot(4,3,5);
plot(log10(pupPSD.freqs),pupPSD.psd,'k','linewidth',2)
hold on
plot(log10(lowfilter),[-1 -1],'r')
LogScale('x',10)
xlabel('f (Hz)'); ylabel('power (dB)')
axis tight
ylim([-2 max(pupPSD.psd)]);

subplot(4,3,6);
imagesc(pupdthist.bins{1},pupdthist.bins{2},pupdthist.counts'); hold on;
colormap(gca,'jet')
%alpha(a,~isnan(pupdthist.counts)')
plot(get(gca,'xlim'),[0 0],'r-')
ColorbarWithAxis([min(min(pupdthist.counts)) max(max(pupdthist.counts))],['counts (au)'])
xlabel('diameter (norm)'); ylabel('dP/dt')
LogScale('x',10);
axis square
ylim([-0.2 0.2]);

subplot(4,3,7:8);
a = imagesc([areadynamics.Xbins areadynamics.Xbins+2*pi],...
    areadynamics.Ybins,[areadynamics.meanZ; areadynamics.meanZ]');
colormap(gca,'jet')
alpha(a,double(~isnan([areadynamics.meanZ; areadynamics.meanZ]')))
ColorbarWithAxis([min(min(areadynamics.meanZ)) max(max(areadynamics.meanZ))],'pupil area (au)')
caxis([min(min(areadynamics.meanZ)) max(max(areadynamics.meanZ))])
ylim([-3 1]);
xlabel('Phase');ylabel('Amplitude')
axis xy

subplot(4,3,9)
a = imagesc(phasedynamics.Xbins,phasedynamics.Ybins,phasedynamics.meanZ');
colormap(gca,'hsv')
alpha(a,double(~isnan(phasedynamics.meanZ')))
ColorbarWithAxis([-pi pi],'phase')
axis xy
xlabel('diameter');ylabel('dP/dt')
axis square

subplot(4,3,10:11);
a = imagesc([dpdtdynamics.Xbins dpdtdynamics.Xbins+2*pi],...
    dpdtdynamics.Ybins,[dpdtdynamics.meanZ; dpdtdynamics.meanZ]');
colormap(gca,'jet')
alpha(a,double(~isnan([dpdtdynamics.meanZ; dpdtdynamics.meanZ]')))
ColorbarWithAxis([min(min(dpdtdynamics.meanZ)) max(max(dpdtdynamics.meanZ))],'dPdt')
caxis([min(min(dpdtdynamics.meanZ)) max(max(dpdtdynamics.meanZ))])
ylim([-3 1]);
xlabel('Phase');ylabel('Amplitude')
axis xy

subplot(4,3,12);
a = imagesc(ampdynamics.Xbins,ampdynamics.Ybins,ampdynamics.meanZ');
colormap(gca,'jet')
alpha(a,double(~isnan(ampdynamics.meanZ')))
ColorbarWithAxis([0 1.2],'power')
%caxis([0 1.2])
axis xy
xlabel('diameter');ylabel('dP/dt')
axis square

NiceSave('Pupil_Stats_PupSpace',figfolder,baseName)

%% EMG STATS
% EMG envelope histogram
EMGhist.bins = linspace(0,20,20);
EMGhist.counts = hist(EMGwhisk.EMGsm,EMGhist.bins);
EMGhist.counts = EMGhist.counts./sum(EMGhist.counts);
EMGhist.logbins = linspace(-1.5,1.5,20);
EMGhist.logcounts = hist(log10(EMGwhisk.EMGsm),EMGhist.logbins);
EMGhist.logcounts = EMGhist.logcounts./sum(EMGhist.logcounts);
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
plot(EMGwhisk.timestamps,EMGwhisk.EMG./max(EMGwhisk.EMG),'color',[0.5 0.5 0.5],'linewidth',1); hold on;
plot(pupildilation.timestamps,EMGwhisk.pupiltime./max(EMGwhisk.pupiltime),'b','linewidth',1);
plot(EMGwhisk.ints.Wh',...
    zeros(size(EMGwhisk.ints.Wh))',...
    'g-','linewidth',1)
ylabel('EMG'); xlabel('time (s)');
axis tight
xlim([100 150]); ylim([0 1]);

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

NiceSave('EMG_Stats',figfolder,baseName)

%% EMG in PupSpace
% Codistribution
pupilEMGdist.bins = {linspace(-0.5,0.5,70),linspace(-1.5,1,70)};
[pupilEMGdist.counts,pupilEMGdist.bins] = hist3([log10(pupildilation.data),log10(EMGwhisk.pupiltime)],pupilEMGdist.bins);
pupilEMGdist.counts = pupilEMGdist.counts./sum(pupilEMGdist.counts(:));

% Conditional histograms
[pupildynamicsEMG.meanZ,pupildynamicsEMG.N,pupildynamicsEMG.Xbins,...
    pupildynamicsEMG.Ybins] = ConditionalHist3( log10(pupildilation.data(1:end-1)),...
    pupildilation.dpdt,log10(EMGwhisk.pupiltime(1:end-1)),...
    'minXY',0,'Xbounds',[-0.5 0.5],'Ybounds',[-0.5 0.5],...
    'numXbins',100,'numYbins',100);
pupildynamicsEMG.occupancy = pupildynamicsEMG.N./pupildilation.samplingRate;

[pWhisk.meanZ,pWhisk.N,pWhisk.Xbins,pWhisk.Ybins] = ConditionalHist3( log10(pupildilation.data(1:end-1)),...
    pupildilation.dpdt,single(pupildilation.iswhisk(1:end-1)),...
    'minXY',0,'Xbounds',[-0.5 0.5],'Ybounds',[-0.5 0.5],...
    'numXbins',100,'numYbins',100);
pupildynamicsEMG.pWhisk = pWhisk.meanZ;

[pupilphaseEMG.meanZ,pupilphaseEMG.N,pupilphaseEMG.Xbins,...
    pupilphaseEMG.Ybins] = ConditionalHist3( lowpupildata.phase,...
    log10(lowpupildata.amp),log10(EMGwhisk.pupiltime),...
    'minXY',0,'Xbounds',[-pi pi],'Ybounds',[-3 1],...
    'numXbins',100,'numYbins',100);
pupilphaseEMG.occupancy = pupilphaseEMG.N./pupildilation.samplingRate;

[pWdynamics.meanZ,pWdynamics.N,pWdynamics.Xbins,...
    pWdynamics.Ybins] = ConditionalHist3( lowpupildata.phase,...
    log10(lowpupildata.amp),single(pupildilation.iswhisk),...
    'minXY',0,'Xbounds',[-pi pi],'Ybounds',[-3 1],...
    'numXbins',100,'numYbins',100);
pupilphaseEMG.pWhisk = pWdynamics.meanZ;

% Wh onsets/offsets in pupil space
whints_pupil = interp1(pupildilation.timestamps,pupildilation.data,EMGwhisk.ints.Wh);
whints_pupildt = interp1(pupildilation.timestamps(1:end-1),pupildilation.dpdt,EMGwhisk.ints.Wh);
whints_pupilphase = interp1(lowpupildata.timestamps,lowpupildata.phase,EMGwhisk.ints.Wh);
whints_pupilamp = interp1(lowpupildata.timestamps,lowpupildata.amp,EMGwhisk.ints.Wh);

tempidx = interp1(lowpupildata.timestamps,lowpupildata.timestamps,EMGwhisk.ints.Wh,'nearest');
whints_maxamp = NaN(size(EMGwhisk.ints.Wh,1),1);
for i = 1:size(EMGwhisk.ints.Wh,1)
    temps = find(lowpupildata.timestamps == tempidx(i,1));
    tempe = find(lowpupildata.timestamps == tempidx(i,2));
    
    if ~isempty(temps) 
        if ~isempty(tempe)
        whints_maxamp(i) = max(lowpupildata.amp(temps:tempe));
        else
        end
    else
    end
end

X = whints_pupilphase(:,1);
Y = log10(whints_pupilamp(:,1));
Y = Y(~isnan(X));
Z = EMGwhisk.Whdurs(:,1);
Z = Z(~isnan(X));
X = X(~isnan(X));

[pupilphaseEMG.meanWhdur,pupilphaseEMG.numWhstarts,~,~] = ConditionalHist3( X,Y,Z,...
    'minXY',0,'Xbounds',[-pi pi],'Ybounds',[-3 1],...
    'numXbins',100,'numYbins',100);

pupilphaseEMG.meanWhdur(pupilphaseEMG.numWhstarts < 0.02) = NaN;
pupilphaseEMG.numWhstarts(pupilphaseEMG.numWhstarts < 0.02) = NaN;
pupilphaseEMG.pWhstarts = pupilphaseEMG.numWhstarts./pupilphaseEMG.N;
pupilphaseEMG.Whstartrate = pupilphaseEMG.numWhstarts./pupilphaseEMG.occupancy;

% Saving to struct
PupEMG.pupilEMGdist = pupilEMGdist;
PupEMG.pupildynamicsEMG = pupildynamicsEMG;
PupEMG.pupilphaseEMG = pupilphaseEMG;
PupEMG.whints_pupil = whints_pupil;
PupEMG.whints_pupildt = whints_pupildt;
PupEMG.whints_pupilphase = whints_pupilphase;
PupEMG.whints_pupilamp = whints_pupilamp;
PupEMG.whints_maxamp = whints_maxamp;
PupEMG.whints_durs = EMGwhisk.Whdurs(:,1);

%% FIGURE 3:
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

subplot(2,2,3);
colormap(gca,[1 1 1;colormap('jet')]);
imagesc(pupilEMGdist.bins{1},pupilEMGdist.bins{2},...
    (pupilEMGdist.counts)'./max(pupilEMGdist.counts(:))); hold on
plot(get(gca,'xlim'),log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],'w--')
LogScale('y',10); axis xy
ColorbarWithAxis([0 0.6],['counts (au)'])
xlabel('Pupil diameter (norm.)'); ylabel('EMG');

subplot(3,2,4);
imagesc(pupildynamicsEMG.Xbins,pupildynamicsEMG.Ybins,pupildynamicsEMG.meanZ');
%colormap(gca,emgcolor)
colormap(gca,[1 1 1; colormap('jet')])
axis xy
xlim([-0.5 0.5]);ylim([-0.5 0.5])
colorbar
ylabel('dP/dt')
title('EMG')
hold on
%caxis([min(min(pupildynamicsEMG.meanZ)) max(max(pupildynamicsEMG.meanZ))])
%LogScale('c',10)
plot(log10(whints_pupil(:,1)),whints_pupildt(:,1),'k.','markersize',2)
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color','k')

subplot(3,2,6);
imagesc(pupildynamicsEMG.Xbins,pupildynamicsEMG.Ybins,pupildynamicsEMG.pWhisk');
%colormap(gca,emgcolor)
colormap(gca,[1 1 1; colormap('jet')])
axis xy
xlim([-0.5 0.5]);ylim([-0.5 0.5])
caxis([-0.01 1])
colorbar
ylabel('dP/dt')
title('p(Whisking)')
hold on
plot(log10(whints_pupil(:,1)),whints_pupildt(:,1),'k.','markersize',2)
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color','k')

NiceSave('EMG_PupSpace1',figfolder,baseName);

%% FIGURE 4:
%distcolor = [1 1 1; makeColorMap([0.7 0.7 0.7],[0 0.5 0],[0.7 0.6,0])];
emgcolor = [1 1 1;makeColorMap([0.5 0.5 0.5],[0 0 0.8])];

figure;

subplot(2,2,1);
a = imagesc([pupilphaseEMG.Xbins pupilphaseEMG.Xbins+2*pi],...
    pupilphaseEMG.Ybins,[pupilphaseEMG.meanZ; pupilphaseEMG.meanZ]');
colormap(gca,'jet')
alpha(a,double(~isnan([pupilphaseEMG.meanZ; pupilphaseEMG.meanZ]')))
colorbar
%ColorbarWithAxis([min(min(pupilphaseEMG.meanZ)) max(max(pupilphaseEMG.meanZ))],'EMG')
%caxis([min(min(pupilphaseEMG.meanZ)) max(max(pupilphaseEMG.meanZ))])
ylim([-2 0.25]);
xlabel('Phase');ylabel('Amplitude')
title('EMG')
axis xy

subplot(2,2,3);
a = imagesc([pupilphaseEMG.Xbins pupilphaseEMG.Xbins+2*pi],...
    pupilphaseEMG.Ybins,[pupilphaseEMG.pWhisk; pupilphaseEMG.pWhisk]');
colormap(gca,'jet')
alpha(a,double(~isnan([pupilphaseEMG.pWhisk; pupilphaseEMG.pWhisk]')))
%ColorbarWithAxis([min(min(pupilphaseEMG.pWhisk)) max(max(pupilphaseEMG.pWhisk))],'pWhisk')
%caxis([min(min(pupilphaseEMG.pWhisk)) max(max(pupilphaseEMG.pWhisk))])
ylim([-2 0.25]);
xlabel('Phase');ylabel('Amplitude')
title('p(Whisking)')
colorbar
axis xy

subplot(2,2,2);
a = imagesc([pupilphaseEMG.Xbins pupilphaseEMG.Xbins+2*pi],...
    pupilphaseEMG.Ybins,log10([pupilphaseEMG.meanWhdur; pupilphaseEMG.meanWhdur])');
colormap(gca,'jet')
alpha(a,double(~isnan([pupilphaseEMG.meanWhdur; pupilphaseEMG.meanWhdur]')))
%ColorbarWithAxis(log10([min(min(pupilphaseEMG.meanWhdur)) max(max(pupilphaseEMG.meanWhdur))]),'EMG Envelope')
ylim([-2 0.25]);  xlim([-pi 3*pi])
caxis([-1 0.5])
LogScale('c',10)
xlabel('Phase');ylabel('Amplitude')
colorbar
axis xy
title('Whisk duration (s)')

subplot(2,2,4);
a = imagesc([pupilphaseEMG.Xbins pupilphaseEMG.Xbins+2*pi],...
    pupilphaseEMG.Ybins,log10([pupilphaseEMG.Whstartrate; pupilphaseEMG.Whstartrate])');
colormap(gca,'jet')
alpha(a,double(~isnan([pupilphaseEMG.Whstartrate; pupilphaseEMG.Whstartrate]')))
%ColorbarWithAxis(log10([min(min(pupilphaseEMG.Whstartrate)) max(max(pupilphaseEMG.Whstartrate))]),'EMG Envelope')
ylim([-2 0.25]); xlim([-pi 3*pi])
LogScale('c',10)
xlabel('Phase');ylabel('Amplitude')
colorbar
axis xy
title('Whisk start rate')

NiceSave('EMG_PupSpace2',figfolder,baseName);

%% ACG/CCG
% Computing x-covariance
acgwin = 200; %s
[pupACG.ACG,pupACG.tlag] = xcov(pupildilation.data,...
    round(acgwin.*pupildilation.samplingRate),'coeff');
pupACG.tlag = pupACG.tlag./pupildilation.samplingRate;

% XCovariance
[pupilEMGcorr.xcorr,pupilEMGcorr.corrlags] = xcov(pupildilation.data,EMGwhisk.pupiltime,'unbiased');
pupilEMGcorr.corrlags = pupilEMGcorr.corrlags.*(1./pupildilation.samplingRate);

% Pupil/EMG at whisking onset
[pwCCG.pupil.WhOn,t_lag,~,alltrans.pupil.WhOn,skippedWh] = EventVsContinousCCG(pupildilation.data,...
    pupildilation.timestamps,EMGwhisk.ints.Wh(:,1),20);
[pwCCG.pupilxy.WhOn(:,1),t_lag,~,~] = EventVsContinousCCG(pupildilation.pupilxy(:,1),...
    pupildilation.timestamps,EMGwhisk.ints.Wh(:,1),20);
[pwCCG.pupilxy.WhOn(:,2),t_lag,~,~] = EventVsContinousCCG(pupildilation.pupilxy(:,2),...
    pupildilation.timestamps,EMGwhisk.ints.Wh(:,1),20);
[pwCCG.EMG.WhOn,t_lag,~,alltrans.EMG.WhOn] = EventVsContinousCCG(EMGwhisk.pupiltime,...
    pupildilation.timestamps,EMGwhisk.ints.Wh(:,1),20);
pwCCG.t_lag = t_lag;

% Pupil/EMG at pupil dilation onset
[wpCCG.pupil.PupOn,t_lag,~,alltrans.pupil.PupOn,skippedWh] = EventVsContinousCCG(pupildilation.data,...
    pupildilation.timestamps,pup_on,20);
[wpCCG.pupilxy.PupOn(:,1),t_lag,~,~] = EventVsContinousCCG(pupildilation.pupilxy(:,1),...
    pupildilation.timestamps,pup_on,20);
[wpCCG.pupilxy.PupOn(:,2),t_lag,~,~] = EventVsContinousCCG(pupildilation.pupilxy(:,2),...
    pupildilation.timestamps,pup_on,20);
[wpCCG.EMG.PupOn,t_lag,~,alltrans.EMG.PupOn] = EventVsContinousCCG(EMGwhisk.pupiltime,...
    pupildilation.timestamps,pup_on,20);
wpCCG.t_lag = t_lag;

% Saving to struct
PupEMG.pupACG = pupACG;
PupEMG.pupilEMGcorr = pupilEMGcorr;
PupEMG.pwCCG = pwCCG;
PupEMG.wpCCG = wpCCG;

%% FIGURE 5:
lagwin = [-5 5];

figure;
subplot(3,2,1);
plot(pupACG.tlag,pupACG.ACG,'k','linewidth',2); hold on;
axis tight
plot(get(gca,'xlim'),[0 0],'r')
xlabel('lag (s)')
ylabel('autocovariance')
title('Pupil diameter')

subplot(3,2,2);
plot(pupilEMGcorr.corrlags,pupilEMGcorr.xcorr','k','linewidth',2); hold on
axis tight
xlim([-40 40])
plot([0 0],get(gca,'ylim'),'b-','linewidth',2)
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

subplot(3,2,4);
plot(t_lag,wpCCG.pupil.PupOn,'k','linewidth',2); hold on;
axis tight
plot([0 0],get(gca,'ylim'),'r-','linewidth',2)
xlim(lagwin)
ylabel('Pupil')
set(gca,'xticklabel',[])

subplot(3,2,6);
plot(t_lag,wpCCG.EMG.PupOn,'b','linewidth',2); hold on;
axis tight
plot([0 0],get(gca,'ylim'),'r-','linewidth',2)
xlim(lagwin)
xlabel('t (s - aligned to PupOn)')
ylabel('EMG')

NiceSave('Pupil_EMG_ACG_CCG',figfolder,baseName)

%% EMG-Pupil Phase coupling
% Phase coupling
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

% Phase-Amp coupling
PhaseAmpCoup = PhaseAmpCouplingByAmp(lowpupildata.phase,log10(lowpupildata.amp),...
    log10(EMGwhisk.pupiltime),10);

% Saving to struct
PupEMG.WPcoupling = WPcoupling;
PupEMG.PhaseAmpCoup = PhaseAmpCoup;

save(savefile,'PupEMG');

%% FIGURE 6:
rwbcolormap = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);
plotx = linspace(-pi,3*pi,100);

figure
subplot(3,2,1)
imagesc(PhaseAmpCoup.phasebins,PhaseAmpCoup.ampbins,PhaseAmpCoup.phaseamphist); hold on;
imagesc(PhaseAmpCoup.phasebins+2*pi,PhaseAmpCoup.ampbins,PhaseAmpCoup.phaseamphist)
plot(PhaseAmpCoup.sig2prefangle,PhaseAmpCoup.ampbins,'.k')
plot(PhaseAmpCoup.sig2prefangle+2*pi,PhaseAmpCoup.ampbins,'.k')
plot(plotx,cos(plotx),'k')
colormap(gca,rwbcolormap)
axis xy
axis tight
ColorbarWithAxis([-0.5 0.5],['Phase-EMG dist'])
caxis([-0.5 0.5])
xlim([-pi 3*pi]); ylim(PhaseAmpCoup.ampbins([1 end]))
xlabel('Pupil phase'); ylabel('Pupil (Z)')

subplot(3,2,2);
plot(log10(wavespec.freqs),coupling,'r','LineWidth',2)
LogScale('x',10);
xlabel('f (Hz)');
axis tight
title('Pupil-EMG phase coupling')

subplot(3,2,3)
plot(PhaseAmpCoup.ampbins,PhaseAmpCoup.sig2powerskew,'k','LineWidth',1)
xlabel('Pupil (Z)');
ylabel('Phase-EMG MI (mrl)')
axis tight

subplot(3,2,4)
histogram(PhaseAmpCoup.sig1amp,PhaseAmpCoup.ampbins)
xlabel('Pupil (Z)');
ylabel('Occupancy')
axis tight

subplot(3,2,6)
histogram(PhaseAmpCoup.sig2amp,PhaseAmpCoup.ampbins)
xlabel('EMG (Z)');
ylabel('Occupancy')
axis tight

subplot(3,2,5)
imagesc(PhaseAmpCoup.ampbins,PhaseAmpCoup.ampbins,PhaseAmpCoup.gasphist)
colormap(gca,'jet');
axis xy
ylabel('EMG (Z)'); xlabel('Pupil (Z)')

NiceSave('Pupil_EMG_PhaseAmpCoup',figfolder,baseName)