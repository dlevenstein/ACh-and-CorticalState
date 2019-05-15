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

%% Get rid of recording start/stop artifact and restrict to spontaneous behavior

%Specifying SPONT whisking
load(fullfile(basePath,[baseName,'.MergePoints.events.mat']),'MergePoints');
sidx = find(startsWith(MergePoints.foldernames,"Spont"));
sponttimes = [MergePoints.timestamps(sidx(1),1) MergePoints.timestamps(sidx(end),2)];


maxtimejump = 1; %s
pupilcycle.amp = NanPadJumps( pupilcycle.timestamps,pupilcycle.amp,maxtimejump );
pupilcycle.phase = NanPadJumps( pupilcycle.timestamps,pupilcycle.phase,maxtimejump );

%% When do whisks happen with respect to UP/DOWN transitions

[ccg,t] = CCG({EMGwhisk.ints.Wh(:,1),EMGwhisk.ints.Wh(:,2),...
    SlowWaves.ints.UP(:,1),SlowWaves.ints.UP(:,2)},...
    [],'norm','rate','binSize',0.02 );

%%

figure
subplot(2,2,1)
plot(t,ccg(:,3,1),'LineWidth',2)
hold on
plot(t,ccg(:,3,2),'LineWidth',2)
plot([0 0],[0 0.35],'k')
xlabel('t lag (s) - DOWN->UP');ylabel('Rate (Wh/s)')
legend('Whisk Onsets','Wh Offsets')

subplot(2,2,2)
plot(t,ccg(:,4,1),'LineWidth',2)
hold on
plot(t,ccg(:,4,2),'LineWidth',2)
plot([0 0],[0 0.35],'k')
xlabel('t lag (s) - UP->DOWN');ylabel('Rate (Wh/s)')
legend('Whisk Onsets','Wh Offsets')
%% UP/DOWN durations
UPDOWN = {'UP','DOWN'};
UDcolor = {'r','b'};
for uu = 1:2
    SlowWaves.durs.(UPDOWN{uu}) = diff(SlowWaves.ints.(UPDOWN{uu}),1,2);
    SlowWaves.durs_np1.(UPDOWN{uu}) = [SlowWaves.durs.(UPDOWN{uu})(2:end);nan];
    SlowWaves.midpoints.(UPDOWN{uu}) = mean(SlowWaves.ints.(UPDOWN{uu}),2);
    SlowWaves.phase.(UPDOWN{uu}) = interp1(pupilcycle.timestamps,pupilcycle.phase,...
        SlowWaves.midpoints.(UPDOWN{uu}),'nearest');
    SlowWaves.amp.(UPDOWN{uu}) = interp1(pupilcycle.timestamps,log10(pupilcycle.amp),...
        SlowWaves.midpoints.(UPDOWN{uu}),'nearest');
    SlowWaves.hipup.(UPDOWN{uu}) = SlowWaves.amp.(UPDOWN{uu})>pupilcycle.pupthresh;
    SlowWaves.lopup.(UPDOWN{uu}) = ~SlowWaves.hipup.(UPDOWN{uu});
    
    [~,SlowWaves.nwh.(UPDOWN{uu})] = ExcludeIntervals(SlowWaves.ints.(UPDOWN{uu}),EMGwhisk.ints.Wh);
    SlowWaves.wh.(UPDOWN{uu}) = ~SlowWaves.nwh.(UPDOWN{uu});
    [~,SlowWaves.nwh_doubles.(UPDOWN{uu})] = ExcludeIntervals([[SlowWaves.ints.(UPDOWN{uu})(1:end-1,1);nan],...
        [SlowWaves.ints.(UPDOWN{uu})(2:end,2);nan]],EMGwhisk.ints.Wh);
%     SlowWaves.meanEMG.(UPDOWN{uu}) = zeros(size(SlowWaves.midpoints.(UPDOWN{uu})));
%     for dd = 1:length(SlowWaves.ints.(UPDOWN{uu})(:,1))
%         dd
%         EMGin = EMGwhisk.timestamps>SlowWaves.ints.(UPDOWN{uu})(dd,1) & EMGwhisk.timestamps<SlowWaves.ints.(UPDOWN{uu})(dd,2);
%         SlowWaves.meanEMG.(UPDOWN{uu})(dd) = mean(EMGwhisk.EMGenvelope(EMGin));
%     end
end


%% Histograms
WHNWH = {'wh','nwh'};
HILO = {'hipup','lopup'};

for uu = 1:2
    for ww = 1:2
        for pp = 1:2
            incondition = SlowWaves.(HILO{pp}).(UPDOWN{uu}) & ...
                SlowWaves.(WHNWH{ww}).(UPDOWN{uu});
            [ConditionalUPDOWN.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww})] = ...
                ConditionalHist(SlowWaves.phase.(UPDOWN{uu})(incondition),...
                log10(SlowWaves.durs.(UPDOWN{uu})(incondition)),...
                'Xbounds',[-pi pi],'Ybounds',[-1.5 1.5],'minX',5,'numXbins',20);
        end
    end
end


%% Figure: durations conditioned on pupil cycle
UDcmap = {makeColorMap([1 1 1],[0.8 0 0]),makeColorMap([1 1 1],[0 0 0.8])};

cosx = linspace(-pi,3*pi,100);
cospamp = [0.6 0.05];
figure
for uu = 1:2; for ww = 1:2; for pp = 1:2
    subplot(4,4,4*(uu-1)+pp+2*(ww-1))
    colormap(gca,UDcmap{uu})
        imagesc(ConditionalUPDOWN.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).Xbins,...
            ConditionalUPDOWN.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).Ybins,...
            ConditionalUPDOWN.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).pYX')
        hold on; axis xy; box off
        imagesc(ConditionalUPDOWN.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).Xbins+2*pi,...
            ConditionalUPDOWN.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).Ybins,...
            ConditionalUPDOWN.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).pYX')
        plot(cosx,(cos(cosx)+1).*cospamp(pp)-1.5,'k')
        LogScale('y',10)
         xlim([-pi 3*pi])
         caxis([0 0.15])
         %colorbar
%         if cc == 1 & uu ==1
%             switch ww
%                 case 1; ylabel('Non-Whisk')
%                 case 2; ylabel('Whisky')
%         elseif cc == 3 & uu ==1
%             ylabel('Whisky')
%         end
        
        
    incondition = SlowWaves.(HILO{pp}).(UPDOWN{uu}) & ...
        SlowWaves.(WHNWH{ww}).(UPDOWN{uu});
    
    subplot(4,4,4*(pp-1)+uu+2*(ww-1)+8)
        plot(SlowWaves.phase.(UPDOWN{uu})(incondition),...
            log10(SlowWaves.durs.(UPDOWN{uu})(incondition)),...
            '.','color',UDcolor{uu})
        hold on
        plot(SlowWaves.phase.(UPDOWN{uu})(incondition)+2*pi,...
            log10(SlowWaves.durs.(UPDOWN{uu})(incondition)),...
            '.','color',UDcolor{uu})
        plot(cosx,(cos(cosx)+1).*cospamp(pp)-1.8,'k')
        xlim([-pi 3*pi]);ylim([-2 1.5])
        LogScale('y',10)
        
end;end;end
%% Figure: UP/DOWN duration return maps
figure
for uu = 1:2
    subplot(2,2,uu)
        plot(log10(SlowWaves.durs.(UPDOWN{uu})(SlowWaves.nwh_doubles.(UPDOWN{uu}))),...
            log10(SlowWaves.durs_np1.(UPDOWN{uu})(SlowWaves.nwh_doubles.(UPDOWN{uu}))),...
            '.','color',UDcolor{uu})
        
        xlim([-2 1.5]);ylim([-2 1.5])
        LogScale('xy',10)
        title(UPDOWN{uu})
        
        
%     subplot(2,2,uu+2)
%         plot(log10(SlowWaves.meanEMG.(UPDOWN{uu})),...
%             log10(SlowWaves.durs.(UPDOWN{uu})),...
%             '.','color',UDcolor{uu})
%         xlim([-2 1.5]);%ylim([-2 1.5])
%         LogScale('xy',10)
%         title(UPDOWN{uu})
%         ylabel('Dur (s)');xlabel('mean EMG')
        
      
end
%%
cosx = linspace(-pi,3*pi,100);
cosmag = [1 0.1 1 0.1];

figure
for uu = 1:2
    
    conditions = {[SlowWaves.hipup.(UPDOWN{uu}) & SlowWaves.nwh.(UPDOWN{uu})],...
        [~SlowWaves.hipup.(UPDOWN{uu}) & SlowWaves.nwh.(UPDOWN{uu})],...
        [SlowWaves.hipup.(UPDOWN{uu}) & ~SlowWaves.nwh.(UPDOWN{uu})],...
        [~SlowWaves.hipup.(UPDOWN{uu}) & ~SlowWaves.nwh.(UPDOWN{uu})]};
    
    for cc = 1:length(conditions)
    subplot(4,4,uu+2*(cc-1))
        plot(SlowWaves.phase.(UPDOWN{uu})(conditions{cc}),...
            log10(SlowWaves.durs.(UPDOWN{uu})(conditions{cc})),...
            '.','color',UDcolor{uu})
        hold on
        plot(2*pi+SlowWaves.phase.(UPDOWN{uu})(conditions{cc}),...
            log10(SlowWaves.durs.(UPDOWN{uu})(conditions{cc})),...
            '.','color',UDcolor{uu})
        plot(cosx,cos(cosx).*cosmag(cc)-1,'k')
        xlim([-pi 3*pi]);ylim([-2 1.5])
        LogScale('y',10)
        title(UPDOWN{uu})
        if cc == 1 & uu ==1
            ylabel('Non-Whisk')
        elseif cc == 3 & uu ==1
            ylabel('Whisky')
        end
    end
end



end

