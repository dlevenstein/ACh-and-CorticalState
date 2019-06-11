function [ EMGdur,EMGdist,pupilphaseEMG,pupildpEMG] = BehaviorAnalysis2(basePath,figfolder)

%Initiate Paths
%reporoot = '/home/dlevenstein/ProjectRepos/ACh-and-CorticalState/';
%reporoot = '/Users/dlevenstein/Project Repos/ACh-and-CorticalState/';
%basePath = '/mnt/proraidDL/Database/WMData/AChPupil/171209_WT_EM1M3/';
%basePath = '/mnt/proraidDL/Database/WMData/AChPupil/180706_WT_EM1M3/';
%basePath = pwd;
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];


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
pupildilation.dpdt = [pupildilation.dpdt; nan]; %To align to timestamps
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

EMGwhisk.iswhisk = InIntervals(EMGwhisk.timestamps,EMGwhisk.ints.Wh);

%% Get rid of recording start/stop artifact

maxtimejump = 1; %s
pupilcycle.amp = NanPadJumps( pupilcycle.timestamps,pupilcycle.amp,maxtimejump );
pupilcycle.phase = NanPadJumps( pupilcycle.timestamps,pupilcycle.phase,maxtimejump );
pupildilation.data = NanPadJumps( pupildilation.timestamps,pupildilation.data,maxtimejump );
pupildilation.dpdt = NanPadJumps( pupildilation.timestamps,pupildilation.dpdt,maxtimejump );
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
EMGwhisk.pup = interp1(pupildilation.timestamps,pupildilation.data,EMGwhisk.timestamps,'nearest');
EMGwhisk.dpdt = interp1(pupildilation.timestamps,pupildilation.dpdt,EMGwhisk.timestamps,'nearest');
EMGwhisk.hipup =  log10(EMGwhisk.pupamp)>pupilcycle.pupthresh;
EMGwhisk.lopup = log10(EMGwhisk.pupamp)<=pupilcycle.pupthresh;

HILO = {'lopup','hipup'};

%% Mean EMG in pupil space

[pupilphaseEMG.meanZ,pupilphaseEMG.N,pupilphaseEMG.Xbins,pupilphaseEMG.Ybins] = ...
    ConditionalHist3( EMGwhisk.pupphase(EMGwhisk.EMGsm~=0),...
    log10(EMGwhisk.pupamp(EMGwhisk.EMGsm~=0)),log10(EMGwhisk.EMGsm(EMGwhisk.EMGsm~=0)),...
    'minXY',500,'Xbounds',[-pi pi],'Ybounds',[-2.25 0.5],...
    'numXbins',30,'numYbins',30);

[pupildpEMG.meanZ,pupildpEMG.N,pupildpEMG.Xbins,pupildpEMG.Ybins] = ...
    ConditionalHist3( log10(EMGwhisk.pup(EMGwhisk.EMGsm~=0)),...
    (EMGwhisk.dpdt(EMGwhisk.EMGsm~=0)),log10(EMGwhisk.EMGsm(EMGwhisk.EMGsm~=0)),...
    'minXY',250,'Xbounds',[-0.5 0.5],'Ybounds',[-0.5 0.5],...
    'numXbins',40,'numYbins',40);




%%
EMGrange = [-1.75 1.1];
for pp= 1:2
    EMGdist.(HILO{pp})  = ConditionalHist(...
        EMGwhisk.pupphase(EMGwhisk.(HILO{pp})),...
        log10(EMGwhisk.EMGsm(EMGwhisk.(HILO{pp}))),...
        'Ybounds',EMGrange,'numYbins',80,'numXbins',20,'Xbounds',[-pi pi]);
    
    EMGdur.(HILO{pp})  = ConditionalHist(...
        EMGwhisk.whisks.pupphase(EMGwhisk.whisks.(HILO{pp})),...
        log10(EMGwhisk.whisks.dur(EMGwhisk.whisks.(HILO{pp}))),...
        'Ybounds',[-1 1.5],'numYbins',30,'numXbins',20,'Xbounds',[-pi pi],...
        'minX',10);
    
    EMGdur.(HILO{pp}).pWhisk = EMGdur.(HILO{pp}).Xhist./(EMGdist.(HILO{pp}).Xhist./EMGwhisk.samplingRate);
    

end

EMGdist.pup  = ConditionalHist( log10(EMGwhisk.pup),...
        log10(EMGwhisk.EMGsm),...
        'numXbins',20,'Xbounds',[-0.25 0.25],'Ybounds',EMGrange,'numYbins',80);
    
EMGdur.pup  = ConditionalHist( log10(EMGwhisk.whisks.pup),...
        log10(EMGwhisk.whisks.dur),...
        'numXbins',20,'Xbounds',[-0.25 0.25],'Ybounds',[-1 1.5],'numYbins',30,...
        'minX',10);
EMGdur.pup.pWhisk = EMGdur.pup.Xhist./(EMGdist.pup.Xhist./EMGwhisk.samplingRate);


%%
cosx = linspace(-pi,pi,100);
cospamp = [0.08 0.8];
figure


subplot(3,2,1)
a = imagesc(pupilphaseEMG.Xbins,pupilphaseEMG.Ybins,pupilphaseEMG.meanZ');
hold on
alpha(a,double(~isnan(pupilphaseEMG.meanZ')))

%imagesc(pupilphaseEMG.Xbins+2*pi,pupilphaseEMG.Ybins,pupilphaseEMG.meanZ')
crameri lapaz
plot([-pi 3*pi],pupilcycle.pupthresh.*[1 1],'w--')
plot(cosx,(cos(cosx)+1).*cospamp(pp)-2,'k')
axis xy
box off
%xlim([-pi 3*pi])
ColorbarWithAxis([-0.7 0.7],'Mean EMG')
LogScale('c',10)
LogScale('y',10)
xlabel('Pupil Phase');ylabel('Pupil Amplitude')

subplot(3,2,2)
a = imagesc(pupildpEMG.Xbins,pupildpEMG.Ybins,pupildpEMG.meanZ');
hold on
alpha(a,double(~isnan(pupildpEMG.meanZ')))
plot(pupildpEMG.Xbins([1 end]),[0 0],'k--')
crameri lapaz
ColorbarWithAxis([-0.7 0.7],'Mean EMG')
LogScale('c',10)
axis xy
box off
LogScale('c',10)
LogScale('x',10)



subplot(3,2,3)
        for pp = 1:2
        imagesc( EMGdist.(HILO{pp}).Xbins+2*pi*(pp-1),...
            EMGdist.(HILO{pp}).Ybins,...
            EMGdist.(HILO{pp}).pYX')
        hold on; axis xy; box off
        plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp)-2,'k')
        end   
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        colorbar
        xlim([-pi 3*pi])
        xlabel('Pupil Phase');ylabel('EMG')
        crameri bilbao
        
subplot(3,2,4)

        imagesc( EMGdist.pup.Xbins,...
            EMGdist.pup.Ybins,...
            EMGdist.pup.pYX')
        hold on; axis xy; box off
        colorbar
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        xlabel('Pupil Size');ylabel('EMG')
   crameri bilbao



subplot(3,2,5)
        for pp = 1:2
        imagesc( EMGdur.(HILO{pp}).Xbins+2*pi*(pp-1),...
            EMGdur.(HILO{pp}).Ybins,...
            EMGdur.(HILO{pp}).pYX')
        hold on; axis xy; box off
        plot(EMGdur.(HILO{pp}).Xbins+2*pi*(pp-1),EMGdur.(HILO{pp}).pWhisk,'k')
        plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp)-1.2,'k')
        end   
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        colorbar
        xlim([-pi 3*pi])
        LogScale('y',10)
        xlabel('Pupil Phase');ylabel('Dur')
        crameri bilbao
        
subplot(3,2,6)

        imagesc( EMGdur.pup.Xbins,...
            EMGdur.pup.Ybins,...
            EMGdur.pup.pYX')
        hold on; axis xy; box off
        plot(EMGdur.pup.Xbins,EMGdur.pup.pWhisk,'k')
        colorbar
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        xlabel('Pupil Size');ylabel('Dur')
        LogScale('y',10)
   crameri bilbao
   
   
        NiceSave('EMGPupBehavior',figfolder,baseName)

%% Example Figure
windows(1,:) = [100 250];
windows(2,:) = bz_RandomWindowInIntervals(pupildilation.timestamps([1 end]),150);
windows(3,:) = bz_RandomWindowInIntervals(pupildilation.timestamps([1 end]),150);

figure;
for ww = 1:3
subplot(3,1,ww);
%plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2); hold on;
scatter(pupildilation.timestamps,pupildilation.data,3,pupilcycle.phase,'filled')
hold on
%scatter(pupilcycle.timestamps,pupilcycle.amp,3,pupilcycle.phase,'filled')
plot(EMGwhisk.timestamps,EMGwhisk.EMG./max(EMGwhisk.EMG),...
    'color',[0.5 0.5 0.5],'linewidth',0.5);
plot(EMGwhisk.timestamps,5*EMGwhisk.EMGsm./(max(EMGwhisk.EMG)),...
    'color','k','linewidth',0.5);
plot(pupilcycle.timestamps(pupilcycle.highpup),3*ones(sum(pupilcycle.highpup),1),'r.')
plot(pupilcycle.timestamps(~pupilcycle.highpup),3*ones(sum(~pupilcycle.highpup),1),'k.')
plot(EMGwhisk.timestamps(EMGwhisk.iswhisk),2.8*ones(sum(EMGwhisk.iswhisk),1),'b.')
%plot(highpupildata.timestamps,highpupildata.data+nanmean(pupildilation.data),'r')
colormap(gca,hsv)
ColorbarWithAxis([min(pupilcycle.phase) max(pupilcycle.phase)],['pupil phase'])
%h1 = plot(get(gca,'xlim'),[0 0],'k-');
xlim(windows(ww,:)); ylim([-0.05 3]);
bz_ScaleBar('s')
ylabel('Pupil');
%set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%legend({'diameter','phase','dPdt'},'location','northeast');

end

NiceSave('BehaviorExamples',figfolder,baseName)