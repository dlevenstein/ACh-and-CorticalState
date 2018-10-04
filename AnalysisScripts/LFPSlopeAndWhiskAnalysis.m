function [] = LFPSlopeAndWhiskAnalysis( basePath,figfolder )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%basePath = '/mnt/proraidDL/Database/WMProbeData/180213_WT_M1M3_LFP_Layers_Pupil_EMG_Pole/180213_WT_M1M3_LFP_Layers_Pupil_EMG_180213_113045';
%basePath = '/home/dlevenstein/ProjectRepos/ACh-and-CorticalState/Dataset/180605_WT_M1M3_LFP_Layers_Pupil_EMG_180605_121846';

basePath = pwd;
figfolder = '/home/dlevenstein/ProjectRepos/ACh-and-CorticalState/AnalysisScripts/AnalysisFigs/LFPSlopeAndWhiskAnalysis';
baseName = bz_BasenameFromBasepath(basePath);
%%
sessionInfo = bz_getSessionInfo(basePath,'noPrompts',true);

%%
pupildilation = bz_LoadBehavior(basePath,'pupildiameter');

%pupildilation.dpdt = 
smoothwin =2;%s
pupildilation.dpdt = diff(smooth(pupildilation.data,smoothwin.*pupildilation.samplingRate,'moving')).*pupildilation.samplingRate;
pupildilation.dpdt = smooth(pupildilation.dpdt,smoothwin.*pupildilation.samplingRate,'moving');

nantimes = isnan(pupildilation.data);
pupildilation.interpdata = interp1(pupildilation.timestamps(~nantimes),...
    pupildilation.data(~nantimes),pupildilation.timestamps);

%% Load the PSS
load([basePath,filesep,baseName,'.PowerSpectrumSlope.lfp.mat'])

%% Load Whisks
EMGwhisk = bz_LoadStates(basePath,'EMGwhisk');

%% Get the pupil phase at each point in time
lowfilter = [0.01 0.1];

pupil4filter = pupildilation;
pupil4filter.data = pupildilation.interpdata(~isnan(pupildilation.interpdata));
%pupil4filter.t
pupil4filter.timestamps = pupil4filter.timestamps(~isnan(pupildilation.interpdata));
lowpupildata = bz_Filter(pupil4filter,'passband',lowfilter,'filter' ,'fir1','order',3);
%highpupildata = bz_Filter(pupil4filter,'passband',highfilter,'filter' ,'fir1');



%Get the pupil phase/power of each whisk start
EMGwhisk.phase = interp1(lowpupildata.timestamps,lowpupildata.phase,...
    EMGwhisk.ints.Wh(:,1),'nearest');
EMGwhisk.power = interp1(lowpupildata.timestamps,log10(lowpupildata.amp),...
    EMGwhisk.ints.Wh(:,1),'nearest');
%Get the closest PSS timepoint to each whisk... (dt issue...)




%%
figure
plot(EMGwhisk.phase,EMGwhisk.power,'.')
hold on
plot(EMGwhisk.phase+2*pi,EMGwhisk.power,'.')



%% Wh-PSS relation by channel
WhiskPSScorr.corr = zeros(size(sessionInfo.AnatGrps.Channels));
WhiskPSScorr.pup = zeros(size(sessionInfo.AnatGrps.Channels));
WhiskPSScorr.dpdt = zeros(size(sessionInfo.AnatGrps.Channels));
for cc = 1:length(sessionInfo.AnatGrps.Channels)
    cc
    channum = sessionInfo.AnatGrps.Channels(cc);
    %channum = 31;
    WhiskPSScorr.channum(cc) = channum;
    WhiskPSScorr.chanpos(cc) = cc;
    %%
    lfp = bz_GetLFP(channum,'basepath',basePath,'noPrompts',true);

    %%
    dt = 0.2;
    winsize = 1;
    [PSS] = bz_PowerSpectrumSlope(lfp,winsize,dt,'showfig',false);

    %%
    PSS.EMG = interp1(EMGwhisk.t,EMGwhisk.EMGenvelope,PSS.timestamps);
    PSS.pupilsize = interp1(pupildilation.timestamps,pupildilation.data,...
        PSS.timestamps,'nearest');
    PSS.dpdt = interp1(pupildilation.timestamps(1:end-1),pupildilation.dpdt,...
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
   
    clear lfp
end

%%
[~, bestchans.EMG] = max(WhiskPSScorr.EMG);
bestchans.EMG = WhiskPSScorr.channum(bestchans.EMG);

%%
figure
subplot(2,3,1:2)
plot(WhiskPSScorr.EMG,-WhiskPSScorr.chanpos,'b','linewidth',2)
hold on
plot(WhiskPSScorr.pup,-WhiskPSScorr.chanpos,'k','linewidth',2)
plot(WhiskPSScorr.dpdt,-WhiskPSScorr.chanpos,'k--','linewidth',1)
legend('EMG','Pupil Area','dpdt','location','eastoutside')
xlabel('PSS Correlation');ylabel('Channel by Depth')
axis tight
box off

NiceSave('PSSCorrbyDepth',figfolder,baseName)

%% Finding best params for PSS

lfp = bz_GetLFP(bestchans.EMG,'basepath',basePath,'noPrompts',true);

dt = 0.05;
winsize = 1;
[PSS] = bz_PowerSpectrumSlope(lfp,winsize,dt,'showfig',true);


%% 
PSS.EMG = interp1(EMGwhisk.t,EMGwhisk.EMGenvelope,PSS.timestamps);
PSS.pupphase = interp1(lowpupildata.timestamps,lowpupildata.phase,...
    PSS.timestamps,'nearest');
PSS.pupmag = interp1(lowpupildata.timestamps,log10(lowpupildata.amp),...
    PSS.timestamps);
%%
figure
subplot(2,2,1)
    plot(log10(PSS.EMG),PSS.data,'k.')
    xlabel('EMG');ylabel('PSS')
subplot(2,2,2)
    scatter(PSS.pupphase,PSS.data,3,log10(PSS.EMG))
    hold on
    scatter(PSS.pupphase+2*pi,PSS.data,3,log10(PSS.EMG))
    colorbar
    xlabel('Pupil Phase');ylabel('PSS')
% subplot(2,2,3)
%     plot(PSS.pupmag,PSS.data,'k.')
%     xlabel('Pupil Magnitude');ylabel('PSS')
subplot(2,2,4)
    scatter(PSS.pupphase,log10(PSS.EMG),3,PSS.data)
    hold on
    scatter(PSS.pupphase+2*pi,log10(PSS.EMG),3,PSS.data)
    colorbar
    xlabel('Pupil Phase');ylabel('EMG')
NiceSave('PupilEMGPSS',figfolder,baseName)
    
%%
EMGwhisk.dur = diff(EMGwhisk.ints.Wh,1,2);
EMGwhisk.longwhisks = EMGwhisk.dur>1;
%% Example Whisk
winsize = 3;
exwhisk = randsample(EMGwhisk.ints.Wh(EMGwhisk.longwhisks,1),1);
viewwin = exwhisk + winsize.*[-1 1];
figure
subplot(3,1,1)
imagesc(PSS.timestamps,log10(PSS.freqs),PSS.specgram)
axis xy
xlim(viewwin)

subplot(6,1,3)
plot(lfp.timestamps,lfp.data)
xlim(viewwin)

subplot(6,1,4)
plot(PSS.timestamps,PSS.data)
xlim(viewwin)

subplot(6,1,5)
plot(EMGwhisk.t,EMGwhisk.EMG,'k')
hold on
plot(EMGwhisk.t,EMGwhisk.EMGenvelope,'b')
xlim(viewwin)
NiceSave('ExWhisk',figfolder,baseName)

%% Get PSS around Whisks
EMGwhisk.numwhisks = length(EMGwhisk.ints.Wh(:,1));

whiskPETH.window = [-5 5]; %s
whiskPETH.windex = 5*PSS.samplingRate; %s
whiskPETH.timestamps = whiskPETH.window(1):(1/PSS.samplingRate):whiskPETH.window(2);
timelockedPSS.data = zeros(length(whiskPETH.timestamps),EMGwhisk.numwhisks);
for ww = 1:EMGwhisk.numwhisks
    PSS.whidx(ww) = find(PSS.timestamps==interp1(PSS.timestamps,PSS.timestamps,EMGwhisk.ints.Wh(ww,1),'nearest'));
    
    timelockedPSS.data(:,ww) = PSS.data(PSS.whidx(ww)-whiskPETH.windex:PSS.whidx(ww)+whiskPETH.windex);
end

%% Whisk Sorts
[~,whisksorts.phase] = sort(EMGwhisk.phase);
[~,whisksorts.dur] = sort(EMGwhisk.dur);


%%
figure
subplot(2,2,1)
    imagesc(whiskPETH.timestamps,[1 EMGwhisk.numwhisks],timelockedPSS.data')
    hold on 
    plot([0 0],[1 EMGwhisk.numwhisks],'b')
    plot(EMGwhisk.dur,1:EMGwhisk.numwhisks,'.r','markersize',1)
    xlim([-2 5])
    colorbar
    %caxis([-1.5 -0.5])
subplot(2,2,2)
    imagesc(whiskPETH.timestamps,[1 EMGwhisk.numwhisks],timelockedPSS.data(:,whisksorts.phase)')
    hold on
    plot(EMGwhisk.dur(whisksorts.phase),1:EMGwhisk.numwhisks,'.r','markersize',1)
    plot([0 0],[0 1],'b')
    xlim([-2 5])
    colorbar
subplot(2,2,3)
    imagesc(whiskPETH.timestamps,[1 EMGwhisk.numwhisks],timelockedPSS.data(:,whisksorts.dur)')
    hold on
    plot(EMGwhisk.dur(whisksorts.dur),1:EMGwhisk.numwhisks,'.r','markersize',1)
    plot([0 0],[1 EMGwhisk.numwhisks],'b')
    xlim([-2 5])
    colorbar

NiceSave('PSSallWHisks',figfolder,baseName)
   
%% Take a look at the channels with dpdt and p
[~, bestchans.pup] = max(WhiskPSScorr.pup);
[~, bestchans.dpdt] = max(WhiskPSScorr.dpdt);
repchans = [sessionInfo.AnatGrps.Channels(bestchans.pup) sessionInfo.AnatGrps.Channels(bestchans.dpdt)];
%repchans = 42;
lfp = bz_GetLFP(sessionInfo.AnatGrps.Channels(repchans(1)),'basepath',basePath,'noPrompts',true);


dt = 0.2;
winsize = 2;
[PSS,specgram] = bz_PowerSpectrumSlope(lfp,winsize,dt,'showfig',true,...
    'saveMat',basePath);

PSS.pupilsize = interp1(pupildilation.timestamps,pupildilation.data,...
    PSS.timestamps,'nearest');
PSS.dpdt = interp1(pupildilation.timestamps(1:end-1),pupildilation.dpdt,...
    PSS.timestamps,'nearest');

%Check residuals
%%
numbins = 15;
bins = linspace(-0.5,0.5,numbins+1);
pupcyclePSS.bincenters = bins(1:end-1)+0.5.*diff(bins([1 2]));
bins([1 end])=[-Inf Inf];
[N,~,~,BINX,BINY] = histcounts2(log10(PSS.pupilsize),PSS.dpdt,...
    bins,bins);
pupcyclePSS.meanPSS = zeros(size(N));
for xx = 1:numbins
    for yy = 1:numbins
        pupcyclePSS.meanPSS(xx,yy) = nanmean(PSS.data(BINX==xx & BINY==yy));
    end
end

nbinthresh = 10;  %Must have more than 10 time windows
pupcyclePSS.meanPSS(N<nbinthresh) = nan;

%% PSS and UP/DOWN

SlowWaves = bz_LoadEvents(basePath,'SlowWaves');


%%
updown = {'DOWN','UP'};
UDcolor = {'b','r'};
for ss = 1:2
    SlowWaves.dur.(updown{ss}) = diff(SlowWaves.ints.(updown{ss}),1,2);
    SlowWaves.midpoint.(updown{ss}) = mean(SlowWaves.ints.(updown{ss}),2);
    SlowWaves.PSS.(updown{ss}) = interp1(PSS.timestamps,PSS.data,SlowWaves.midpoint.(updown{ss}));
end

%%
numbins = 30;
PSShist.bins = linspace(-2,0,numbins);

PSShist.hist = hist(PSS.data,PSShist.bins);

%%
exwinsize = 300;
exwin = bz_RandomWindowInIntervals(PSS.timestamps([1 end])',exwinsize);

winsize = 8;
maxtime = pupildilation.timestamps(...
    (pupildilation.timestamps>exwin(1) & pupildilation.timestamps<exwin(2)) & pupildilation.data==...
    max(pupildilation.data(pupildilation.timestamps>exwin(1) & pupildilation.timestamps<exwin(2))));
subsamplewin = maxtime(1)+[-0.5 0.5].*winsize;

mintime = pupildilation.timestamps(...
    (pupildilation.timestamps>exwin(1) & pupildilation.timestamps<exwin(2)) & pupildilation.data==...
    min(pupildilation.data(pupildilation.timestamps>exwin(1) & pupildilation.timestamps<exwin(2))));
subsamplewin2 = mintime(1)+[-0.5 0.5].*winsize;

figure
subplot(4,1,1)
    imagesc(specgram.timestamps,log2(specgram.freqs),specgram.amp)
    axis xy
    xlim(exwin)
    ylabel({'Specgram','f (Hz)'})
    LogScale('y',2)
    
subplot(6,1,3)
    plot(PSS.timestamps,PSS.data,'k')
    axis tight
    xlim(exwin)
    box off
    xlabel('t (s)')
    ylabel('PSS')
subplot(6,1,4)
    plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2)
    hold on
    plot(pupildilation.timestamps(1:end-1),pupildilation.dpdt,'k')
    axis tight
    xlim(exwin)
    box off
    ylabel('Pupil')
        plot(exwin,[0 0],'k--')
        plot(subsamplewin,3.*ones(size(subsamplewin)),'r','linewidth',2)
        plot(subsamplewin2,3.*ones(size(subsamplewin2)),'r','linewidth',2)
        
        
    

        
        
    subplot(6,2,11)
        bz_MultiLFPPlot(lfp,'timewin',subsamplewin)
        hold on
       % plot(slow.timestamps,10000.*slow.data,'g','linewidth',1)   
        ylabel('High Pupil')
        xlim(subsamplewin)
        
    subplot(6,2,12)
        bz_MultiLFPPlot(lfp,'timewin',subsamplewin2)
        hold on
       % plot(slow.timestamps,10000.*slow.data,'g','linewidth',1)    
        ylabel('Low Pupil')
        xlim(subsamplewin2)

subplot(6,2,9)
    plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2)
    hold on
    plot(pupildilation.timestamps(1:end-1),pupildilation.dpdt,'k')
    plot(PSS.timestamps,PSS.data)
    axis tight
    xlim(subsamplewin)
    box off
    ylabel('Pupil')
        plot(subsamplewin,[0 0],'k--')
        
subplot(6,2,10)
    plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2)
    hold on
    plot(pupildilation.timestamps(1:end-1),pupildilation.dpdt,'k')
    plot(PSS.timestamps,PSS.data)
    axis tight
    xlim(subsamplewin2)
    box off
    ylabel('Pupil')
        plot(subsamplewin2,[0 0],'k--')
        
NiceSave('PSSexample',figfolder,baseName,'tiff')

%%
figure
subplot(3,3,1)
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

subplot(6,3,12)
plot(PSShist.bins,PSShist.hist,'k','linewidth',2)
box off
xlabel('PSS')
    axis tight
    xlim([-1.75 -0.25])
    
subplot(3,2,2)
for ss = 1:2
    plot(SlowWaves.PSS.(updown{ss}),log10(SlowWaves.dur.(updown{ss})),'.','color',UDcolor{ss},'markersize',3)
    hold on
end
xlabel('PSS');ylabel('Dur (s)')
axis tight
box off
LogScale('y',10)
legend(updown{:},'location','eastoutside')


subplot(3,3,4)
    plot(log10(PSS.pupilsize),PSS.data,'k.','markersize',1)
    xlabel('Pupil Area (med^-^1)');ylabel('PSS')
    box off
    axis tight
subplot(3,3,5)
    plot((PSS.dpdt),PSS.data,'k.','markersize',1)
    xlabel('dpdt (med^-^1s^-^1)');ylabel('PSS')
    box off
    axis tight

NiceSave('PSSandUPDOWN',figfolder,baseName,'tiff')


    %%
samplewin = bz_RandomWindowInIntervals(lfp.timestamps([1 end])',300);
    

winsize = 10;
maxtime = pupildilation.timestamps(...
    (pupildilation.timestamps>samplewin(1) & pupildilation.timestamps<samplewin(2)) & pupildilation.data==...
    max(pupildilation.data(pupildilation.timestamps>samplewin(1) & pupildilation.timestamps<samplewin(2))));
subsamplewin = maxtime(1)+[-0.5 0.5].*winsize;

mintime = pupildilation.timestamps(...
    (pupildilation.timestamps>samplewin(1) & pupildilation.timestamps<samplewin(2)) & pupildilation.data==...
    min(pupildilation.data(pupildilation.timestamps>samplewin(1) & pupildilation.timestamps<samplewin(2))));
subsamplewin2 = mintime(1)+[-0.5 0.5].*winsize;

figure
% subplot(3,3,1)
%     plot(slopepupilcorr.dpdt,'k','linewidth',0.5)
%     hold on
%     plot(slopepupilcorr.p,'k','linewidth',2)
%     plot(repchans(1),slopepupilcorr.dpdt(repchans(1)),'o')
%     plot(repchans(2),slopepupilcorr.p(repchans(2)),'o')
%     xlabel('Channel (sorted by depth)')
%     axis tight
%     box off

for cc = 1%:2    
subplot(6,6,12+cc)
    plot(log10(PSS.pupilsize),PSS.data(:,cc),'k.','markersize',3)
    xlabel('Pupil Area (med^-^1)');ylabel('Spectrum Slope')
    box off
    axis tight
subplot(6,6,18+cc)
    plot((PSS.dpdt),PSS.data(:,cc),'k.','markersize',1)
    xlabel('dpdt (med^-^1s^-^1)');ylabel('Spectrum Slope')
    box off
    axis tight
end
    
% subplot(6,1,4)
%     bz_MultiLFPPlot(lfp,'timewin',samplewin)
subplot(6,6,9:12)
    plot(PSS.timestamps,PSS.data)
    axis tight
    xlim(samplewin)
    box off
    xlabel('t (s)')
    ylabel('PSS')
    
subplot(6,6,3:6)
    plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2)
    hold on
    plot(pupildilation.timestamps(1:end-1),pupildilation.dpdt,'k')
    axis tight
    xlim(samplewin)
    box off
    ylabel('Pupil')
        plot(samplewin,[0 0],'k--')
        plot(subsamplewin,3.*ones(size(subsamplewin)),'r','linewidth',2)
        plot(subsamplewin2,3.*ones(size(subsamplewin2)),'r','linewidth',2)
        
        
    subplot(6,3,11:12)
        bz_MultiLFPPlot(lfp,'timewin',samplewin)
        hold on
       % plot(slow.timestamps,10000.*slow.data,'g','linewidth',1)   
        ylabel('High Pupil')
        xlim(subsamplewin)
        
    subplot(6,3,17:18)
        bz_MultiLFPPlot(lfp,'timewin',subsamplewin2)
        hold on
       % plot(slow.timestamps,10000.*slow.data,'g','linewidth',1)    
        ylabel('Low Pupil')
        xlim(subsamplewin2)

subplot(6,3,8:9)
    plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2)
    hold on
    plot(pupildilation.timestamps(1:end-1),pupildilation.dpdt,'k')
    plot(PSS.timestamps,PSS.data)
    axis tight
    xlim(subsamplewin)
    box off
    ylabel('Pupil')
        plot(subsamplewin,[0 0],'k--')
        
        
subplot(6,3,14:15)
    plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2)
    hold on
    plot(pupildilation.timestamps(1:end-1),pupildilation.dpdt,'k')
    plot(PSS.timestamps,PSS.data)
    axis tight
    xlim(subsamplewin2)
    box off
    ylabel('Pupil')
        plot(subsamplewin2,[0 0],'k--')
        
NiceSave('PSSandPupil',figfolder,baseName)
end

