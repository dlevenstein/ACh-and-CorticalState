function [ pupcycleUPDOWN,pupphaseUD,islowhist,pupbyislow ] = UPDOWNandPupilAnalysis(basePath,figfolder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%
%basePath = '/mnt/proraidDL/Database/WMProbeData/180213_WT_M1M3_LFP_Layers_Pupil_EMG_Pole/180213_WT_M1M3_LFP_Layers_Pupil_EMG_180213_113045';
%basePath = '/home/dlevenstein/ProjectRepos/ACh-and-CorticalState/Dataset/180605_WT_M1M3_LFP_Layers_Pupil_EMG_180605_121846';
basePath = '/Users/dlevenstein/Desktop/180703_WT_EM1M3';
%figfolder = '/mnt/data1/Dropbox/research/Current Projects/S1State/AnalysisScripts/figures/UPDOWNandPupilAnalysis';
%figfolder = '/home/dlevenstein/ProjectRepos/ACh-and-CorticalState/AnalysisScripts/AnalysisFigs/UPDOWNandPupilAnalysis';

%%
baseName = bz_BasenameFromBasepath(basePath);
recparms = bz_getSessionInfo(basePath,'noPrompts',true);

%% Detect Slow Waves
%CTXChans = recparms.SpkGrps.Channels(23:46);
[SlowWaves] = DetectSlowWaves(basePath,'noSpikes',true,'noPrompts',true,...
    'NREMInts',[0 Inf]);%,'CTXChans',CTXChans);

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


%% UP/DOWN durations
UPDOWN = {'UP','DOWN'};
UDcolor = {'r','b'};
for uu = 1:2
    SlowWaves.durs.(UPDOWN{uu}) = diff(SlowWaves.ints.(UPDOWN{uu}),1,2);
    SlowWaves.midpoints.(UPDOWN{uu}) = mean(SlowWaves.ints.(UPDOWN{uu}),2);
    SlowWaves.phase.(UPDOWN{uu}) = interp1(pupilcycle.timestamps,pupilcycle.phase,...
        SlowWaves.midpoints.(UPDOWN{uu}),'nearest');
    SlowWaves.amp.(UPDOWN{uu}) = interp1(pupilcycle.timestamps,log10(pupilcycle.amp),...
        SlowWaves.midpoints.(UPDOWN{uu}),'nearest');
    SlowWaves.highpupil.(UPDOWN{uu}) = SlowWaves.amp.(UPDOWN{uu})>pupilcycle.pupthresh;
    
    SlowWaves.meanEMG.(UPDOWN{uu}) = zeros(size(SlowWaves.midpoints.(UPDOWN{uu})));
    for dd = 1:length(SlowWaves.ints.(UPDOWN{uu})(:,1))
        dd
        EMGin = EMGwhisk.timestamps>SlowWaves.ints.(UPDOWN{uu})(dd,1) & EMGwhisk.timestamps<SlowWaves.ints.(UPDOWN{uu})(dd,2);
        SlowWaves.meanEMG.(UPDOWN{uu})(dd) = mean(EMGwhisk.EMGenvelope(EMGin));
    end
end


%% Figure: UP/DOWN duration return maps
figure
for uu = 1:2
    subplot(2,2,uu)
        plot(log10(SlowWaves.durs.(UPDOWN{uu})(1:end-1)),...
            log10(SlowWaves.durs.(UPDOWN{uu})(2:end)),...
            '.','color',UDcolor{uu})
        xlim([-2 1.5]);ylim([-2 1.5])
        LogScale('xy',10)
        title(UPDOWN{uu})
        
        
    subplot(2,2,uu+2)
        plot(log10(SlowWaves.meanEMG.(UPDOWN{uu})),...
            log10(SlowWaves.durs.(UPDOWN{uu})),...
            '.','color',UDcolor{uu})
        xlim([-2 1.5]);%ylim([-2 1.5])
        LogScale('xy',10)
        title(UPDOWN{uu})
        ylabel('Dur (s)');xlabel('mean EMG')
        
      
end
%%
figure
for uu = 1:2
    subplot(2,2,uu)
        plot(SlowWaves.phase.(UPDOWN{uu})(SlowWaves.highpupil.(UPDOWN{uu})),...
            log10(SlowWaves.durs.(UPDOWN{uu})(SlowWaves.highpupil.(UPDOWN{uu}))),...
            '.','color',UDcolor{uu})
        hold on
        plot(2*pi+SlowWaves.phase.(UPDOWN{uu})(SlowWaves.highpupil.(UPDOWN{uu})),...
            log10(SlowWaves.durs.(UPDOWN{uu})(SlowWaves.highpupil.(UPDOWN{uu}))),...
            '.','color',UDcolor{uu})
        xlim([-pi 3*pi]);ylim([-2 1.5])
        LogScale('y',10)
        title(UPDOWN{uu})
        
    subplot(2,2,uu+2)
        plot(SlowWaves.phase.(UPDOWN{uu})(~SlowWaves.highpupil.(UPDOWN{uu})),...
            log10(SlowWaves.durs.(UPDOWN{uu})(~SlowWaves.highpupil.(UPDOWN{uu}))),...
            '.','color',UDcolor{uu})
        hold on
        plot(2*pi+SlowWaves.phase.(UPDOWN{uu})(~SlowWaves.highpupil.(UPDOWN{uu})),...
            log10(SlowWaves.durs.(UPDOWN{uu})(~SlowWaves.highpupil.(UPDOWN{uu}))),...
            '.','color',UDcolor{uu})
        xlim([-pi 3*pi]);ylim([-2 1.5])
        LogScale('y',10)
        title(UPDOWN{uu})  
end



end

