function [ output_args ] = UPDOWNandPupilAnalysis(basePath,figfolder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%
basePath = '/mnt/proraidDL/Database/WMProbeData/180213_WT_M1M3_LFP_Layers_Pupil_EMG_Pole/180213_WT_M1M3_LFP_Layers_Pupil_EMG_180213_113045';
figfolder = '/mnt/data1/Dropbox/research/Current Projects/S1State/AnalysisScripts/figures/UPDOWNandPupilAnalysis';
figfolder = '/home/dlevenstein/ProjectRepos/ACh-and-CorticalState/AnalysisScripts/AnalysisFigs/UPDOWNandPupilAnalysis';

%%
baseName = bz_BasenameFromBasepath(basePath);
recparms = bz_getSessionInfo(basePath,'noPrompts',true);

%% Detect Slow Waves
CTXChans = recparms.SpkGrps.Channels(23:46);
[SlowWaves] = DetectSlowWaves(basePath,'noSpikes',true,'noPrompts',true,...
    'NREMInts',[0 Inf],'CTXChans',CTXChans);

%%
pupildilation = bz_LoadBehavior(basePath,'pupildiameter');

%pupildilation.dpdt = 
smoothwin =1;%s
pupildilation.dpdt = diff(smooth(pupildilation.data,smoothwin.*pupildilation.samplingRate,'moving')).*pupildilation.samplingRate;
pupildilation.dpdt = smooth(pupildt,smoothwin.*pupildilation.samplingRate,'moving');

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
numbins = 10;
bins = linspace(-0.5,0.5,numbins);
pDOWN = hist3([log10(SlowWaves.pupil.DOWN),SlowWaves.dpdt.DOWN],...
    {bins,bins});
pPUP = hist3([log10(pupildilation.data(1:end-1)),pupildilation.dpdt],...
    {bins,bins})./pupildilation.samplingRate;
pDOWNPUP = pDOWN./pPUP;

%%
winsize = 300;
samplewin = bz_RandomWindowInIntervals(pupildilation.timestamps([1 end])',winsize);

figure
subplot(2,2,1)
plot(log10(SlowWaves.pupil.DOWN),log10(SlowWaves.dur.DOWN),'b.')
hold on
plot(log10(SlowWaves.pupil.UP),log10(SlowWaves.dur.UP),'r.')
xlabel('Pupil Area (med^-^1)');ylabel('Duration (s)')
LogScale('xy',10)

subplot(2,2,2)
plot((SlowWaves.dpdt.DOWN),log10(SlowWaves.dur.DOWN),'b.')
hold on
plot((SlowWaves.dpdt.UP),log10(SlowWaves.dur.UP),'r.')
xlabel('dpdt (med^-^1s^-^1)');ylabel('Duration (s)')
LogScale('y',10)

subplot(4,1,3)
plot(pupildilation.timestamps,log10(pupildilation.data),'k')
hold on
plot(pupildilation.timestamps(1:end-1),pupildilation.dpdt,'r')
xlim(samplewin)

subplot(4,1,4)
plot(SlowWaves.midpoint.DOWN ,log10(SlowWaves.dur.DOWN),'b.')
hold on
plot(SlowWaves.midpoint.UP ,log10(SlowWaves.dur.UP),'r.')
xlim(samplewin)




%%

figure
h = imagesc(bins,bins,pDOWNPUP');
set(h,'AlphaData',pPUP'>1);
hold on
plot(log10(SlowWaves.pupil.DOWN),SlowWaves.dpdt.DOWN,'k.','markersize',2)
axis xy
colorbar
caxis([0 3])
xlabel('Pupil Area (med^-^1)')
ylabel('dp/dt')
LogScale('x',10)

%% DOWN-Pupil Phase Coupling
pupilwave = pupildilation;
pupilwave.data = (pupilwave.interpdata);
%pupilwave = bz_WaveSpec(pupilwave);
%%
[freqs,synchcoupling,ratepowercorr,...
    spikephasemag,spikephaseangle,popcellind,cellpopidx,...
    spikephasesig,ratepowersig]...
    = GenSpikeLFPCoupling({SlowWaves.timestamps},log10(pupilwave.data),...
    'sf_LFP',pupilwave.samplingRate,'frange',[0.001 1],'ncyc',3);

%%
infraslowpup = [0.005 0.05];
filtered = bz_Filter(pupilwave,'passband',infraslowpup,...
    'filter','fir1','order',1);
filtered.normamp = zscore(log10(filtered.amp));

%% DOWN state by phase/power in infraslow
SlowWaves.islowpup.phase.DOWN = interp1(filtered.timestamps,filtered.phase,...
    SlowWaves.midpoint.DOWN,'nearest');
SlowWaves.islowpup.amp.DOWN = interp1(filtered.timestamps,filtered.normamp,...
    SlowWaves.midpoint.DOWN,'nearest');

numbins = 6;
islowhist.phasebins = linspace(-pi,pi,numbins+1);
islowhist.phasebins = islowhist.phasebins(1:end-1)+0.5.*diff(islowhist.phasebins([1 2]));
islowhist.ampbins = linspace(-1.5,1.5,numbins);

islowhist.pDOWN = hist3([SlowWaves.islowpup.phase.DOWN,(SlowWaves.islowpup.amp.DOWN)],...
    {islowhist.phasebins,islowhist.ampbins});
islowhist.pPUP = hist3([filtered.phase,filtered.normamp],...
    {islowhist.phasebins,islowhist.ampbins})./filtered.samplingRate;
islowhist.pDOWNPUP = islowhist.pDOWN./islowhist.pPUP;

islowhist.pDOWN_marg = hist(SlowWaves.islowpup.phase.DOWN,islowhist.phasebins);
islowhist.pPUP_marg = hist(filtered.phase,islowhist.phasebins)./filtered.samplingRate;
islowhist.pDOWNPUP_marg = islowhist.pDOWN_marg./islowhist.pPUP_marg;

%% Mean Pupil as f'n of phase
nphasebins = 25;
pupbyislow.phasebins = linspace(-pi,pi,nphasebins);
%pupbyislow.phasebins = pupbyislow.phasebins(1:end-1)+0.5.*diff(pupbyislow.phasebins([1 2]));
pupbyislow.ampbins = islowhist.ampbins;

filtered.nearestphasebin = interp1(pupbyislow.phasebins,pupbyislow.phasebins,...
    filtered.phase,'nearest');
filtered.nearestampbin = interp1(pupbyislow.ampbins,pupbyislow.ampbins,...
    filtered.normamp,'nearest');

pupbyislow.meanpupbyphase = zeros(size(pupbyislow.phasebins));
pupbyislow.meandpdtbyphase = zeros(size(pupbyislow.phasebins));
pupbyislow.stdpupbyphase = zeros(size(pupbyislow.phasebins));
pupbyislow.meanpupbyphaseamp = zeros(length(pupbyislow.phasebins),...
    length(pupbyislow.ampbins));
for pp = 1:nphasebins
    pupbyislow.meanpupbyphase(pp) =...
        nanmean(pupildilation.data(filtered.nearestphasebin==pupbyislow.phasebins(pp)));
    pupbyislow.stdpupbyphase(pp) =...
        nanstd(pupildilation.data(filtered.nearestphasebin==pupbyislow.phasebins(pp)));

    pupbyislow.meandpdtbyphase(pp) =...
        nanmean(pupildilation.dpdt(filtered.nearestphasebin(1:end-1)==pupbyislow.phasebins(pp)));
    
    for aa = 1:length(pupbyislow.ampbins)
        pupbyislow.meanpupbyphaseamp(pp,aa) =...
            nanmean(pupildilation.data...
            (filtered.nearestphasebin==pupbyislow.phasebins(pp) & ...
            filtered.nearestampbin==pupbyislow.ampbins(aa)));
    end
    
end

%%
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

subplot(6,2,6)
    bar(islowhist.phasebins,islowhist.pDOWNPUP_marg)
    hold on
    bar(islowhist.phasebins+2.*pi,islowhist.pDOWNPUP_marg)
    xlim(pi.*[-1 3])
    ylabel('DOWN rate (s^-^1)')
    
subplot(6,2,8)
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
    xlim(pi.*[-1 3])
    xlabel(['Phase: (',num2str(infraslowpup(1)),'-',num2str(infraslowpup(2)),' Hz)'])
    ylabel('Pupil Area (med^-^1)')
    
subplot(3,2,1)
    plot(log10(freqs),spikephasemag)
    hold on
    plot(log10(freqs(spikephasesig>3)),spikephasemag(spikephasesig>3),'o')
    LogScale('x',10)
    xlabel('f (Hz)');ylabel('Pupil-DOWN state coupling')

subplot(4,1,4)
    plot(pupildilation.timestamps,pupildilation.data,'k')
    hold on
    plot(filtered.timestamps,filtered.data,'r')
    plot(SlowWaves.timestamps,ones(size(SlowWaves.timestamps)),'b.')
    xlim(samplewin)

NiceSave('InfraSlowPupil',figfolder,baseName)
%% Example Window: 
repchans = recparms.SpkGrps.Channels(5:5:45);
lfp = bz_GetLFP(repchans,...
    'basepath',basePath,'noPrompts',true);

%%
winsize = 400;
samplewin = bz_RandomWindowInIntervals(pupildilation.timestamps([1 end])',winsize);

winsize = 15;
maxtime = pupildilation.timestamps(...
    (pupildilation.timestamps>samplewin(1) & pupildilation.timestamps<samplewin(2)) & pupildilation.data==...
    max(pupildilation.data(pupildilation.timestamps>samplewin(1) & pupildilation.timestamps<samplewin(2))));
subsamplewin = maxtime+[-0.5 0.5].*winsize;

mintime = pupildilation.timestamps(...
    (pupildilation.timestamps>samplewin(1) & pupildilation.timestamps<samplewin(2)) & pupildilation.data==...
    min(pupildilation.data(pupildilation.timestamps>samplewin(1) & pupildilation.timestamps<samplewin(2))));
subsamplewin2 = mintime(1)+[-0.5 0.5].*winsize;


figure
    subplot(3,1,2)
        bz_MultiLFPPlot(lfp,'timewin',samplewin)
        hold on
        plot(SlowWaves.timestamps,zeros(size(SlowWaves.timestamps)),'ob')
        ylabel('High Pupil')
        xlim(subsamplewin)
        
    subplot(3,1,3)
        bz_MultiLFPPlot(lfp,'timewin',subsamplewin2)
        hold on
        plot(SlowWaves.timestamps,zeros(size(SlowWaves.timestamps)),'ob')
        ylabel('Low Pupil')
        xlim(subsamplewin2)
        
    subplot(4,1,1)
        plot(pupildilation.timestamps,pupildilation.data,'k','LineWidth',2)
        hold on
        plot(filtered.timestamps,filtered.data,'r')
        plot(pupildilation.timestamps(1:end-1),pupildilation.dpdt,'k')
        plot(SlowWaves.timestamps,ones(size(SlowWaves.timestamps)),'b.')
        plot(subsamplewin,3.*ones(size(subsamplewin)),'r','linewidth',2)
        plot(subsamplewin2,3.*ones(size(subsamplewin2)),'r','linewidth',2)
        xlim(samplewin)
        legend('Pupil Area','Filtered','dp/dt','DOWN states','location','northwest')
        ylabel('Pupil Metrics')
NiceSave('examples',figfolder,baseName)
%%

end

