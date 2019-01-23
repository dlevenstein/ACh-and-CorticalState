function [ pupcycleUPDOWN,pupphaseUD,islowhist,pupbyislow ] = UPDOWNandPupilAnalysis(basePath,figfolder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%
%basePath = '/mnt/proraidDL/Database/WMProbeData/180213_WT_M1M3_LFP_Layers_Pupil_EMG_Pole/180213_WT_M1M3_LFP_Layers_Pupil_EMG_180213_113045';
%basePath = '/home/dlevenstein/ProjectRepos/ACh-and-CorticalState/Dataset/180605_WT_M1M3_LFP_Layers_Pupil_EMG_180605_121846';
%basePath = pwd;
%figfolder = '/mnt/data1/Dropbox/research/Current Projects/S1State/AnalysisScripts/figures/UPDOWNandPupilAnalysis';
%figfolder = '/home/dlevenstein/ProjectRepos/ACh-and-CorticalState/AnalysisScripts/AnalysisFigs/UPDOWNandPupilAnalysis';

%%
baseName = bz_BasenameFromBasepath(basePath);
recparms = bz_getSessionInfo(basePath,'noPrompts',true);

%% Detect Slow Waves
%CTXChans = recparms.SpkGrps.Channels(23:46);
[SlowWaves] = DetectSlowWaves(basePath,'noSpikes',true,'noPrompts',true,...
    'NREMInts',[0 Inf]);%,'CTXChans',CTXChans);

%%
pupildilation = bz_LoadBehavior(basePath,'pupildiameter');

%pupildilation.dpdt = 
smoothwin =2;%s
pupildilation.dpdt = diff(smooth(pupildilation.data,smoothwin.*pupildilation.samplingRate,'moving')).*pupildilation.samplingRate;
pupildilation.dpdt = smooth(pupildilation.dpdt,smoothwin.*pupildilation.samplingRate,'moving');

nantimes = isnan(pupildilation.data);
pupildilation.interpdata = interp1(pupildilation.timestamps(~nantimes),...
    pupildilation.data(~nantimes),pupildilation.timestamps);
%% Get the pupil diameter at each UP/DOWN state
states = {'UP','DOWN'};
for ss = 1:length(states)
    SlowWaves.dur.(states{ss}) = diff(SlowWaves.ints.(states{ss}),1,2);
    SlowWaves.midpoint.(states{ss}) = mean(SlowWaves.ints.(states{ss}),2);
    SlowWaves.pupil.(states{ss}) = interp1(pupildilation.timestamps,pupildilation.data,SlowWaves.midpoint.(states{ss}),'nearest');
    SlowWaves.dpdt.(states{ss}) = interp1(pupildilation.timestamps(1:end-1),pupildilation.dpdt,SlowWaves.midpoint.(states{ss}),'nearest');
end


%% DOWN state probability as function of pupil
numbins = 8;
bins = linspace(-0.5,0.5,numbins+1);
pupcycleUPDOWN.bincenters = bins(1:end-1)+0.5.*diff(bins([1 2]));
bins([1 end])=[-Inf Inf];
pDOWN = histcounts2(log10(SlowWaves.pupil.DOWN),SlowWaves.dpdt.DOWN,...
    bins,bins);
pPUP = histcounts2(log10(pupildilation.data(1:end-1)),pupildilation.dpdt,...
    bins,bins)./pupildilation.samplingRate;
pupcycleUPDOWN.pDOWN = pDOWN./pPUP;

puptimethresh = 1; %s
pupcycleUPDOWN.pDOWN(pPUP<puptimethresh)=nan;

%% UP state duration as function of pupil

[N,~,~,BINX,BINY] = histcounts2(log10(SlowWaves.pupil.UP),SlowWaves.dpdt.UP,...
    bins,bins);
pupcycleUPDOWN.logUPdur = zeros(size(pDOWN));
for xx = 1:numbins
    for yy = 1:numbins
        pupcycleUPDOWN.logUPdur(xx,yy) = nanmean(log10(SlowWaves.dur.UP(BINX==xx & BINY==yy)));
    end
end    

nUPthresh = 2;
pupcycleUPDOWN.logUPdur(N<nUPthresh) = nan;

%%
winsize = 300;
samplewin = bz_RandomWindowInIntervals(pupildilation.timestamps([1 end])',winsize);

figure
subplot(3,2,1)
plot(log10(SlowWaves.pupil.DOWN),log10(SlowWaves.dur.DOWN),'b.')
hold on
plot(log10(SlowWaves.pupil.UP),log10(SlowWaves.dur.UP),'r.')
xlabel('Pupil Area (med^-^1)');ylabel('Duration (s)')
LogScale('xy',10)
box off 
axis tight

subplot(3,2,3)
plot((SlowWaves.dpdt.DOWN),log10(SlowWaves.dur.DOWN),'b.')
hold on
plot((SlowWaves.dpdt.UP),log10(SlowWaves.dur.UP),'r.')

xlabel('dpdt (med^-^1s^-^1)');ylabel('Duration (s)')
LogScale('y',10)
box off
axis tight
plot([0 0],get(gca,'ylim'),'k--')

subplot(6,1,6)
plot(pupildilation.timestamps,(pupildilation.data),'k')
hold on
plot(pupildilation.timestamps(1:end-1),pupildilation.dpdt,'r')
plot(samplewin,[0 0],'k--')
xlim(samplewin)

subplot(6,1,5)
plot(SlowWaves.midpoint.DOWN ,log10(SlowWaves.dur.DOWN),'b.')
hold on
plot(SlowWaves.midpoint.UP ,log10(SlowWaves.dur.UP),'r.')
xlim(samplewin)
box off


subplot(3,3,6)
h = imagesc(pupcycleUPDOWN.bincenters,pupcycleUPDOWN.bincenters,pupcycleUPDOWN.pDOWN');
set(h,'AlphaData',~isnan(pupcycleUPDOWN.pDOWN'));
hold on
plot(log10(SlowWaves.pupil.DOWN),SlowWaves.dpdt.DOWN,'k.','markersize',2)
plot(pupcycleUPDOWN.bincenters([1 end]),[0 0],'k--')
axis xy
colorbar
caxis([0 2.5])
xlabel('Pupil Area (med^-^1)')
ylabel('dp/dt')
LogScale('x',10)
title('DOWN rate')

subplot(3,3,3)
h = imagesc(pupcycleUPDOWN.bincenters,pupcycleUPDOWN.bincenters,pupcycleUPDOWN.logUPdur');
set(h,'AlphaData',~isnan(pupcycleUPDOWN.logUPdur'));
hold on
plot(pupcycleUPDOWN.bincenters([1 end]),[0 0],'k--')
LogScale('x',10)
axis xy
colorbar
xlabel('Pupil Area (med^-^1)')
ylabel('dp/dt')
LogScale('x',10)
title('(log) UP duration')

NiceSave('PupilandDpdt',figfolder,baseName)

%% DOWN-Pupil Phase Coupling
frange = [0.001 1];
pupilwave = pupildilation;
pupilwave.data = log10(pupilwave.interpdata);
pupilwavespec = bz_WaveSpec(pupilwave,'frange',[0.001 1],'ncyc',3);
%% Wavelet Coupling with DOWN incidence
[pupphaseUD.freqs,synchcoupling,ratepowercorr,...
    pupphaseUD.phaseDN.mag,spikephaseangle,popcellind,cellpopidx,...
    pupphaseUD.phaseDN.sig,ratepowersig]...
    = GenSpikeLFPCoupling({SlowWaves.timestamps},(pupilwave.data),...
    'sf_LFP',pupilwave.samplingRate,'frange',frange,'ncyc',3,...
    'jittersig',true);

%% Filter in Infraslow
infraslowpup = [0.005 0.05];
islow = bz_Filter(pupilwave,'passband',infraslowpup,...
    'filter','fir1','order',1);
islow.normamp = zscore(log10(islow.amp));

%% DOWN state by phase/power in infraslow
SlowWaves.islowpup.phase.DOWN = interp1(islow.timestamps,islow.phase,...
    SlowWaves.midpoint.DOWN,'nearest');
SlowWaves.islowpup.amp.DOWN = interp1(islow.timestamps,islow.normamp,...
    SlowWaves.midpoint.DOWN,'nearest');

numbins = 6;
islowhist.frange = infraslowpup;
islowhist.phasebins = linspace(-pi,pi,numbins+1);
islowhist.phasebins = islowhist.phasebins(1:end-1)+0.5.*diff(islowhist.phasebins([1 2]));
islowhist.ampbins = linspace(-1.5,1.5,numbins);

islowhist.pDOWN = hist3([SlowWaves.islowpup.phase.DOWN,(SlowWaves.islowpup.amp.DOWN)],...
    {islowhist.phasebins,islowhist.ampbins});
islowhist.pPUP = hist3([islow.phase,islow.normamp],...
    {islowhist.phasebins,islowhist.ampbins})./islow.samplingRate;
islowhist.pDOWNPUP = islowhist.pDOWN./islowhist.pPUP;

islowhist.pDOWN_marg = hist(SlowWaves.islowpup.phase.DOWN,islowhist.phasebins);
islowhist.pPUP_marg = hist(islow.phase,islowhist.phasebins)./islow.samplingRate;
islowhist.pDOWNPUP_marg = islowhist.pDOWN_marg./islowhist.pPUP_marg;

puptimethresh = 1;
islowhist.pDOWNPUP(islowhist.pPUP<puptimethresh) = nan;

%% Mean Pupil as f'n of phase
nphasebins = 25;
pupbyislow.phasebins = linspace(-pi,pi,nphasebins);
%pupbyislow.phasebins = pupbyislow.phasebins(1:end-1)+0.5.*diff(pupbyislow.phasebins([1 2]));
pupbyislow.ampbins = islowhist.ampbins;

islow.nearestphasebin = interp1(pupbyislow.phasebins,pupbyislow.phasebins,...
    islow.phase,'nearest');
islow.nearestampbin = interp1(pupbyislow.ampbins,pupbyislow.ampbins,...
    islow.normamp,'nearest');

pupbyislow.meanpupbyphase = zeros(size(pupbyislow.phasebins));
pupbyislow.meandpdtbyphase = zeros(size(pupbyislow.phasebins));
pupbyislow.stdpupbyphase = zeros(size(pupbyislow.phasebins));
pupbyislow.meanpupbyphaseamp = zeros(length(pupbyislow.phasebins),...
    length(pupbyislow.ampbins));
for pp = 1:nphasebins
    pupbyislow.meanpupbyphase(pp) =...
        nanmean(pupildilation.data(islow.nearestphasebin==pupbyislow.phasebins(pp)));
    pupbyislow.stdpupbyphase(pp) =...
        nanstd(pupildilation.data(islow.nearestphasebin==pupbyislow.phasebins(pp)));

    pupbyislow.meandpdtbyphase(pp) =...
        nanmean(pupildilation.dpdt(islow.nearestphasebin(1:end-1)==pupbyislow.phasebins(pp)));
    
    for aa = 1:length(pupbyislow.ampbins)
        pupbyislow.meanpupbyphaseamp(pp,aa) =...
            nanmean(pupildilation.data...
            (islow.nearestphasebin==pupbyislow.phasebins(pp) & ...
            islow.nearestampbin==pupbyislow.ampbins(aa)));
    end
    
end

%% Coupling with UP duration

SlowWaves.pupphase.UP = zeros(length(SlowWaves.midpoint.UP),length(pupilwavespec.freqs));
pupphaseUD.phaseUPdur.corr = zeros(size(pupilwavespec.freqs));
pupphaseUD.phaseUPdur.pval = zeros(size(pupilwavespec.freqs));
%pupphaseUD.freqs = pupilwavespec.freqs;

for ff = 1:length(pupilwavespec.freqs)
    SlowWaves.pupphase.UP(:,ff) = interp1(pupilwavespec.timestamps,angle(pupilwavespec.data(:,ff)),...
        SlowWaves.midpoint.UP,'nearest');
    [pupphaseUD.phaseUPdur.corr(ff), pupphaseUD.phaseUPdur.pval(ff)] = circ_corrcl(SlowWaves.pupphase.UP(:,ff), log10(SlowWaves.dur.UP));
end

%% Slow Range
slowpup = [0.05 0.5];
slow = bz_Filter(pupilwave,'passband',slowpup,...
    'filter','fir1','order',2);
slow.normamp = zscore(log10(slow.amp));

%%
for ss = 1:length(states)
    SlowWaves.slowpup.phase.(states{ss}) = interp1(slow.timestamps,slow.phase,...
        SlowWaves.midpoint.(states{ss}),'nearest');
    SlowWaves.slowpup.amp.(states{ss}) = interp1(slow.timestamps,slow.normamp,...
        SlowWaves.midpoint.(states{ss}),'nearest');
end
%% Slow Figure
statecolors = {'r','b'};

winsize = 100;
samplewin = bz_RandomWindowInIntervals(pupildilation.timestamps([1 end])',winsize);

figure
subplot(3,2,1)
    plot(log10(pupphaseUD.freqs),pupphaseUD.phaseUPdur.corr,'k')
    hold on
    plot(log10(pupphaseUD.freqs(pupphaseUD.phaseUPdur.pval<0.05)),pupphaseUD.phaseUPdur.corr(pupphaseUD.phaseUPdur.pval<0.05),...
        'ro','markersize',4)
    axis tight
    box off
    LogScale('x',10)
    xlabel('Pupil f (Hz)');ylabel('Phase-Dur_U_P Corr')
    
%
subplot(3,2,2)
for ss = 2:-1:1
    plot(SlowWaves.slowpup.phase.(states{ss}),log10(SlowWaves.dur.(states{ss})),'.','color',statecolors{ss})
    hold on
    plot(SlowWaves.slowpup.phase.(states{ss})+2.*pi,log10(SlowWaves.dur.(states{ss})),'.','color',statecolors{ss})
    xlabel(['Phase: (',num2str(slowpup(1)),'-',num2str(slowpup(2)),' Hz)'])
    ylabel('Duration')
end
axis tight
box off
LogScale('y',10)
%end

subplot(4,1,4)
    plot(pupildilation.timestamps,pupildilation.data,'k')
    hold on
    plot(slow.timestamps,5.*slow.data,'r')
    plot(SlowWaves.timestamps,ones(size(SlowWaves.timestamps)),'b.')
    xlim(samplewin)
NiceSave('SlowPupil',figfolder,baseName)
%% iSlow Figure
winsize = 500;
samplewin = bz_RandomWindowInIntervals(pupildilation.timestamps([1 end])',winsize);

figure
subplot(3,2,2)
    h = imagesc(islowhist.phasebins,islowhist.ampbins,islowhist.pDOWNPUP');
    set(h,'AlphaData',islowhist.pPUP'>1);
    hold on
    h = imagesc(islowhist.phasebins+2.*pi,islowhist.ampbins,islowhist.pDOWNPUP');
    set(h,'AlphaData',islowhist.pPUP'>1);
    plot(SlowWaves.islowpup.phase.DOWN,(SlowWaves.islowpup.amp.DOWN),'k.','markersize',2)
    hold on
    plot(SlowWaves.islowpup.phase.DOWN+2*pi,(SlowWaves.islowpup.amp.DOWN),'k.','markersize',2)
    axis xy
    %colorbar
    clim([0 2])
    xlim(pi.*[-1 3])
    ylabel(['Power: (',num2str(infraslowpup(1)),'-',num2str(infraslowpup(2)),' Hz)'])
%caxis([0 3])

subplot(8,2,8)
    bar(islowhist.phasebins,islowhist.pDOWNPUP_marg,'facecolor','b')
    hold on
    bar(islowhist.phasebins+2.*pi,islowhist.pDOWNPUP_marg,'facecolor','b')
    axis tight
    xlim(pi.*[-1 3])
    ylabel('DOWN rate (s^-^1)')
    box off
    
subplot(8,2,10)
    %plot(linspace(-pi,3.*pi,100),cos(linspace(-pi,3.*pi,100)),'k')
    plot(pupbyislow.phasebins,pupbyislow.meanpupbyphase,'k')
    hold on
    plot(pupbyislow.phasebins+2.*pi,pupbyislow.meanpupbyphase,'k')
    errorshade(pupbyislow.phasebins,pupbyislow.meanpupbyphase,...
        pupbyislow.stdpupbyphase,pupbyislow.stdpupbyphase,'k','scalar')
    errorshade(pupbyislow.phasebins+2.*pi,pupbyislow.meanpupbyphase,...
        pupbyislow.stdpupbyphase,pupbyislow.stdpupbyphase,'k','scalar')
    
    %plot(pupbyislow.phasebins+2.*pi,pupbyislow.meandpdtbyphase,'k')
    
   % plot(pupbyislow.phasebins,pupbyislow.meanpupbyphaseamp(:,2:end-1),'k')
   box off
   axis tight
    xlim(pi.*[-1 3])
    xlabel(['Phase: (',num2str(infraslowpup(1)),'-',num2str(infraslowpup(2)),' Hz)'])
    ylabel('Pupil Area (med^-^1)')
    
subplot(3,2,1)
    plot(log10(pupphaseUD.freqs),pupphaseUD.phaseDN.mag,'k')
    hold on
    plot(log10(pupphaseUD.freqs(pupphaseUD.phaseDN.sig>3)),pupphaseUD.phaseDN.mag(pupphaseUD.phaseDN.sig>3),'bo','markersize',4)
    axis tight
    box off
    LogScale('x',10)
    xlabel('Pupil f (Hz)');ylabel('Pupil-DOWN state coupling')

subplot(4,1,4)
    plot(pupildilation.timestamps,pupildilation.data,'k')
    hold on
    plot(islow.timestamps,islow.data,'m')
    plot(SlowWaves.timestamps,ones(size(SlowWaves.timestamps)),'b.')
    xlim(samplewin)

NiceSave('InfraSlowPupil',figfolder,baseName)
%% Example Window: 
repchans = recparms.AnatGrps.Channels(5:5:45);
lfp = bz_GetLFP(repchans,...
    'basepath',basePath,'noPrompts',true);

%%
winsize = 300;
samplewin = bz_RandomWindowInIntervals(pupildilation.timestamps([1 end])',winsize);

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
    subplot(3,4,5:7)
        bz_MultiLFPPlot(lfp,'timewin',samplewin)
        hold on
        plot(SlowWaves.timestamps,zeros(size(SlowWaves.timestamps)),'ob','markersize',5)
       % plot(slow.timestamps,10000.*slow.data,'g','linewidth',1)   
        ylabel('High Pupil')
        xlim(subsamplewin)
        
    subplot(3,4,9:11)
        bz_MultiLFPPlot(lfp,'timewin',subsamplewin2)
        hold on
        plot(SlowWaves.timestamps,zeros(size(SlowWaves.timestamps)),'ob','markersize',5)
       % plot(slow.timestamps,10000.*slow.data,'g','linewidth',1)    
        ylabel('Low Pupil')
        xlim(subsamplewin2)
        
    subplot(4,1,1)
        plot(pupildilation.timestamps,pupildilation.data,'k','LineWidth',2)
        hold on
        plot(pupildilation.timestamps(1:end-1),pupildilation.dpdt,'k')
        plot(islow.timestamps,islow.data-0.5,'m')
        plot(slow.timestamps,slow.data-1,'g')    
        plot(SlowWaves.timestamps,-2.*ones(size(SlowWaves.timestamps)),'b.')
        plot(samplewin,[0 0],'k--')
        plot(subsamplewin,3.*ones(size(subsamplewin)),'r','linewidth',2)
        plot(subsamplewin2,3.*ones(size(subsamplewin2)),'r','linewidth',2)
        axis tight
        xlim(samplewin)
        legend('Pupil Area','dp/dt','Filtered iSlow','Filtered Slow','DOWN states','location','eastoutside')
        box off
        ylabel('Pupil')
NiceSave('examples',figfolder,baseName)
%%

end

