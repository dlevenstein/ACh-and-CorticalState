function [PSShist,ConditionalUPDOWN,ReturnHist  ] = UPDOWNDurationAnalysis(basePath,figfolder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%
%basePath = '/mnt/proraidDL/Database/WMProbeData/180213_WT_M1M3_LFP_Layers_Pupil_EMG_Pole/180213_WT_M1M3_LFP_Layers_Pupil_EMG_180213_113045';
%basePath = '/home/dlevenstein/ProjectRepos/ACh-and-CorticalState/Dataset/180605_WT_M1M3_LFP_Layers_Pupil_EMG_180605_121846';
%basePath = '/Users/dlevenstein/Desktop/171206_WT_EM1M3';
%basePath = '/Users/dlevenstein/Desktop/180704_KO_EM1M3';
%basePath = pwd;
%figfolder = '/Users/dlevenstein/Project Repos/ACh-and-CorticalState/AnalysisScripts/AnalysisFigs/UPDOWNDurationAnalysis';
%figfolder = '/home/dlevenstein/ProjectRepos/ACh-and-CorticalState/AnalysisScripts/AnalysisFigs/UPDOWNandPupilAnalysis';
%figfolder = fullfile(basePath,'AnalysisFigures');
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
PSS.winsize = PSpecSlope.Shortwin.movingwin(1);

%%
UPDOWN = {'UP','DOWN'};
UDcolor = {'r','b'};
WHNWH = {'wh','nwh'};
HILO = {'hipup','lopup'};
%% Get rid of recording start/stop artifact and restrict to spontaneous behavior

%Restricting SPONT UP/DOWNs
load(fullfile(basePath,[baseName,'.MergePoints.events.mat']),'MergePoints');
sidx = find(startsWith(MergePoints.foldernames,"Spont"));
sponttimes = [MergePoints.timestamps(sidx(1),1) MergePoints.timestamps(sidx(end),2)];

for uu = 1:2
    inspont = InIntervals(SlowWaves.ints.(UPDOWN{uu})(:,1),sponttimes);
    SlowWaves.ints.(UPDOWN{uu})(~inspont,:) = [];
end
inspont = InIntervals(PSS.timestamps,sponttimes);
PSS.timestamps = PSS.timestamps(inspont);
PSS.data = PSS.data(inspont);

maxtimejump = 1; %s
pupilcycle.amp = NanPadJumps( pupilcycle.timestamps,pupilcycle.amp,maxtimejump );
pupilcycle.phase = NanPadJumps( pupilcycle.timestamps,pupilcycle.phase,maxtimejump );
pupildilation.data = NanPadJumps( pupildilation.timestamps,pupildilation.data,maxtimejump );

PSS.wh = InIntervals(PSS.timestamps,EMGwhisk.ints.Wh);
PSS.nwh = ~InIntervals(PSS.timestamps,EMGwhisk.ints.Wh);
PSS.hipup = interp1(pupilcycle.timestamps,single(pupilcycle.highpup),PSS.timestamps,'nearest')==1;
PSS.lopup = interp1(pupilcycle.timestamps,single(~pupilcycle.highpup),PSS.timestamps,'nearest')==1;
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
%% Get all the features of UP/DOWN states

for uu = 1:2
    %Time/Dur Stuff
    SlowWaves.durs.(UPDOWN{uu}) = diff(SlowWaves.ints.(UPDOWN{uu}),1,2);
    SlowWaves.durs_np1.(UPDOWN{uu}) = [SlowWaves.durs.(UPDOWN{uu})(2:end);nan];
    SlowWaves.midpoints.(UPDOWN{uu}) = mean(SlowWaves.ints.(UPDOWN{uu}),2);
    
    %Pupil Stuff
    SlowWaves.phase.(UPDOWN{uu}) = interp1(pupilcycle.timestamps,pupilcycle.phase,...
        SlowWaves.midpoints.(UPDOWN{uu}),'nearest');
    SlowWaves.amp.(UPDOWN{uu}) = interp1(pupilcycle.timestamps,log10(pupilcycle.amp),...
        SlowWaves.midpoints.(UPDOWN{uu}),'nearest');
    SlowWaves.pupsize.(UPDOWN{uu}) = interp1(pupildilation.timestamps,log10(pupildilation.data),...
        SlowWaves.midpoints.(UPDOWN{uu}),'nearest');
    SlowWaves.hipup.(UPDOWN{uu}) = SlowWaves.amp.(UPDOWN{uu})>pupilcycle.pupthresh;
    SlowWaves.lopup.(UPDOWN{uu}) = ~SlowWaves.hipup.(UPDOWN{uu});
    
    %Whisk Stuff
    [~,SlowWaves.nwh.(UPDOWN{uu})] = ExcludeIntervals(SlowWaves.ints.(UPDOWN{uu}),EMGwhisk.ints.Wh);
    SlowWaves.wh.(UPDOWN{uu}) = ~SlowWaves.nwh.(UPDOWN{uu});
    [~,SlowWaves.doubles.nwh.(UPDOWN{uu})] = ExcludeIntervals([[SlowWaves.ints.(UPDOWN{uu})(1:end-1,1);nan],...
        [SlowWaves.ints.(UPDOWN{uu})(2:end,2);nan]],EMGwhisk.ints.Wh);
    SlowWaves.doubles.wh.(UPDOWN{uu}) = ~SlowWaves.doubles.nwh.(UPDOWN{uu});
%     SlowWaves.meanEMG.(UPDOWN{uu}) = zeros(size(SlowWaves.midpoints.(UPDOWN{uu})));
%     for dd = 1:length(SlowWaves.ints.(UPDOWN{uu})(:,1))
%         dd
%         EMGin = EMGwhisk.timestamps>SlowWaves.ints.(UPDOWN{uu})(dd,1) & EMGwhisk.timestamps<SlowWaves.ints.(UPDOWN{uu})(dd,2);
%         SlowWaves.meanEMG.(UPDOWN{uu})(dd) = mean(EMGwhisk.EMGenvelope(EMGin));
%     end

    %PSS Stuff
    SlowWaves.PSS.(UPDOWN{uu}) = interp1(PSS.timestamps,PSS.data,...
        SlowWaves.midpoints.(UPDOWN{uu}),'nearest');
end


%% Histograms


for uu = 1:2
    for ww = 1:2
        
        %UD Duration by PSS
        [ConditionalUPDOWN.(UPDOWN{uu}).PSS.(WHNWH{ww})] = ...
            ConditionalHist((SlowWaves.PSS.(UPDOWN{uu})(SlowWaves.(WHNWH{ww}).(UPDOWN{uu}))),...
            log10(SlowWaves.durs.(UPDOWN{uu})(SlowWaves.(WHNWH{ww}).(UPDOWN{uu}))),...
            'Xbounds',[-3.25 -0.75],'Ybounds',[-1.5 1],'minX',15,'numXbins',20);
        
        %PSS Histogram by wh/nwh
        PSShist.bins = ConditionalUPDOWN.(UPDOWN{uu}).PSS.(WHNWH{ww}).Xbins;
        PSShist.allpup.(WHNWH{ww}) = hist(PSS.data(PSS.(WHNWH{ww})),PSShist.bins);
        PSShist.allpup.(WHNWH{ww}) = PSShist.allpup.(WHNWH{ww})./sum(PSShist.allpup.(WHNWH{ww}));
        
        %Hi/lo pupil amplitude split
        for pp = 1:2
            incondition = SlowWaves.(HILO{pp}).(UPDOWN{uu}) & ...
                SlowWaves.(WHNWH{ww}).(UPDOWN{uu});
            %UP/DOWN durations by pupil phase
            [ConditionalUPDOWN.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww})] = ...
                ConditionalHist(SlowWaves.phase.(UPDOWN{uu})(incondition),...
                log10(SlowWaves.durs.(UPDOWN{uu})(incondition)),...
                'Xbounds',[-pi pi],'Ybounds',[-1.5 1],'minX',1,'numXbins',20);
            %U/D Duration Return maps split by pupil/whisk groups
            ReturnHist.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}) = hist3(...
                [log10(SlowWaves.durs.(UPDOWN{uu})(incondition)),...
                log10(SlowWaves.durs_np1.(UPDOWN{uu})(incondition))],...
                {ConditionalUPDOWN.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).Ybins,...
                ConditionalUPDOWN.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).Ybins});
            ReturnHist.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}) = ...
                ReturnHist.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww})./...
                sum(ReturnHist.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}));
            
            %UD Duration by PSS - split by hi/lo pup
            [ConditionalUPDOWN.(UPDOWN{uu}).PSS.(WHNWH{ww}).(HILO{pp})] = ...
                ConditionalHist((SlowWaves.PSS.(UPDOWN{uu})(incondition)),...
                log10(SlowWaves.durs.(UPDOWN{uu})(incondition)),...
                'Xbounds',[-3.25 -0.75],'Ybounds',[-1.5 1],'minX',15,'numXbins',20);
            
            PSShist.(HILO{pp}).(WHNWH{ww}) = hist(PSS.data(PSS.(WHNWH{ww})&PSS.(HILO{pp})),PSShist.bins);
            PSShist.(HILO{pp}).(WHNWH{ww}) = PSShist.(HILO{pp}).(WHNWH{ww})./sum(PSShist.(HILO{pp}).(WHNWH{ww}));
        end
        
        %UD Duration Return maps split by whisk 
        ReturnHist.(UPDOWN{uu}).allpup.(WHNWH{ww}) = hist3(...
            [log10(SlowWaves.durs.(UPDOWN{uu})(SlowWaves.doubles.(WHNWH{ww}).(UPDOWN{uu}))),...
            log10(SlowWaves.durs_np1.(UPDOWN{uu})(SlowWaves.doubles.(WHNWH{ww}).(UPDOWN{uu})))],...
            {ConditionalUPDOWN.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).Ybins,...
            ConditionalUPDOWN.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).Ybins});
        ReturnHist.(UPDOWN{uu}).allpup.(WHNWH{ww}) = ...
            ReturnHist.(UPDOWN{uu}).allpup.(WHNWH{ww})./...
            sum(ReturnHist.(UPDOWN{uu}).allpup.(WHNWH{ww}));
        
        %UD Duation by pupil size
        [ConditionalUPDOWN.(UPDOWN{uu}).pupmag.(WHNWH{ww})] = ...
            ConditionalHist((SlowWaves.pupsize.(UPDOWN{uu})(SlowWaves.(WHNWH{ww}).(UPDOWN{uu}))),...
            log10(SlowWaves.durs.(UPDOWN{uu})(SlowWaves.(WHNWH{ww}).(UPDOWN{uu}))),...
            'Xbounds',[-0.2 0.2],'Ybounds',[-1.5 1],'minX',1,'numXbins',10);
       

        
    end

end

%%
UDcmap = {makeColorMap([1 1 1],[0.8 0 0]),makeColorMap([1 1 1],[0 0 0.8])};

figure
for uu = 1:2; for ww = 1:2;
    subplot(3,2,1+(uu-1)*2+(ww-1))
    colormap(gca,UDcmap{uu})
    imagesc(ConditionalUPDOWN.(UPDOWN{uu}).PSS.(WHNWH{ww}).Xbins,...
        ConditionalUPDOWN.(UPDOWN{uu}).PSS.(WHNWH{ww}).Ybins,...
        ConditionalUPDOWN.(UPDOWN{uu}).PSS.(WHNWH{ww}).pYX')
    hold on
    plot(ConditionalUPDOWN.(UPDOWN{uu}).PSS.(WHNWH{ww}).Xbins([1 end]),...
        log10(PSS.winsize*[1 1]),'k--')
    plot(PSShist.bins,PSShist.allpup.(WHNWH{ww})*8+...
        ConditionalUPDOWN.(UPDOWN{uu}).PSS.(WHNWH{ww}).Ybins(1)-0.5,'k','linewidth',1)
    axis xy; box off;axis tight
   % xlim([-1.5 1])
    LogScale('y',10)
    colorbar
    caxis([0 0.125])
    if uu==2
    xlabel('PSS')
    else
        title((WHNWH{ww}))
    end
    if ww==1
        ylabel('Duration (s)')
    end
       
end;end
NiceSave('UDbyPSS',figfolder,baseName)

%%
figure
for uu = 1:2; for ww = 1:2; for pp=1:2;
    subplot(4,2,1+(uu-1)*2+(ww-1)+(pp-1).*4)
    colormap(gca,UDcmap{uu})
    imagesc(ConditionalUPDOWN.(UPDOWN{uu}).PSS.(WHNWH{ww}).(HILO{pp}).Xbins,...
        ConditionalUPDOWN.(UPDOWN{uu}).PSS.(WHNWH{ww}).(HILO{pp}).Ybins,...
        ConditionalUPDOWN.(UPDOWN{uu}).PSS.(WHNWH{ww}).(HILO{pp}).pYX')
    hold on
    plot(ConditionalUPDOWN.(UPDOWN{uu}).PSS.(WHNWH{ww}).(HILO{pp}).Xbins([1 end]),...
        log10(PSS.winsize*[1 1]),'k--')
    plot(PSShist.bins,PSShist.(HILO{pp}).(WHNWH{ww})*10+...
        ConditionalUPDOWN.(UPDOWN{uu}).PSS.(WHNWH{ww}).(HILO{pp}).Ybins(1))
    axis xy; box off
   % xlim([-1.5 1])
    LogScale('y',10)
    if uu==2
    xlabel('PSS')
    end
NiceSave('UDbyPSS_hilopup',figfolder,baseName)
       
end;end;end

%% Figure: durations conditioned on pupil cycle

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
             if uu==2
                xlabel('Pup Phase')
             end
    if ww == 1 &pp==1
        ylabel('Dur')
    end
         %colorbar
%         if cc == 1 & uu ==1
%             switch ww
%                 case 1; ylabel('Non-Whisk')
%                 case 2; ylabel('Whisky')
%         elseif cc == 3 & uu ==1
%             ylabel('Whisky')
%         end
        
        
%     incondition = SlowWaves.(HILO{pp}).(UPDOWN{uu}) & ...
%         SlowWaves.(WHNWH{ww}).(UPDOWN{uu});
%     
%     subplot(4,4,4*(uu-1)+pp+2*(ww-1)+8)
%         plot(SlowWaves.phase.(UPDOWN{uu})(incondition),...
%             log10(SlowWaves.durs.(UPDOWN{uu})(incondition)),...
%             '.','color',UDcolor{uu})
%         hold on
%         plot(SlowWaves.phase.(UPDOWN{uu})(incondition)+2*pi,...
%             log10(SlowWaves.durs.(UPDOWN{uu})(incondition)),...
%             '.','color',UDcolor{uu})
%         plot(cosx,(cos(cosx)+1).*cospamp(pp)-1.8,'k')
%         xlim([-pi 3*pi]);ylim([-2 1.5])
%         LogScale('y',10)'

    subplot(4,4,9+(uu-1)*4+(ww-1)*2)
    colormap(gca,UDcmap{uu})
    imagesc(ConditionalUPDOWN.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).Ybins,...
        ConditionalUPDOWN.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).Ybins,...
        ReturnHist.(UPDOWN{uu}).allpup.(WHNWH{ww}))
    axis xy
    LogScale('xy',10)
    if uu==2
        xlabel('Dur (n-1)')
    end
    if ww == 1
        ylabel('Dur')
    end
    
    subplot(4,4,10+(uu-1)*4+(ww-1)*2)
    colormap(gca,UDcmap{uu})
    imagesc(ConditionalUPDOWN.(UPDOWN{uu}).pupmag.(WHNWH{ww}).Xbins,...
        ConditionalUPDOWN.(UPDOWN{uu}).pupmag.(WHNWH{ww}).Ybins,...
        ConditionalUPDOWN.(UPDOWN{uu}).pupmag.(WHNWH{ww}).pYX')
    axis xy
   % xlim([-1.5 1])
    LogScale('y',10)
    if uu==2
    xlabel('Pupil')
    end
        
end;end;end

NiceSave('UDbyPupilCycle',figfolder,baseName)



%% Figure: UP/DOWN duration return maps
figure
for uu = 1:2
    subplot(2,2,uu)
        plot(log10(SlowWaves.durs.(UPDOWN{uu})(SlowWaves.doubles.nwh.(UPDOWN{uu}))),...
            log10(SlowWaves.durs_np1.(UPDOWN{uu})(SlowWaves.doubles.nwh.(UPDOWN{uu}))),...
            '.','color',UDcolor{uu})
        
        xlim([-2 1.5]);ylim([-2 1.5])
        LogScale('xy',10)
        title(UPDOWN{uu})
end    

for uu = 1:2; for ww = 1:2; for pp = 1:2
    subplot(4,4,4*(uu-1)+pp+2*(ww-1))
        imagesc(ConditionalUPDOWN.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).Ybins,...
            ConditionalUPDOWN.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).Ybins,...
            ReturnHist.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}))
        axis xy
        LogScale('xy',10)
end;end;end

for uu = 1:2
subplot(4,4,8+uu)
colormap(gca,UDcmap{uu})
imagesc(ConditionalUPDOWN.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).Ybins,...
    ConditionalUPDOWN.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).Ybins,...
    ReturnHist.(UPDOWN{uu}).allpup.nwh)
axis xy
LogScale('xy',10)
end
%     subplot(2,2,uu+2)
%         plot(log10(SlowWaves.meanEMG.(UPDOWN{uu})),...
%             log10(SlowWaves.durs.(UPDOWN{uu})),...
%             '.','color',UDcolor{uu})
%         xlim([-2 1.5]);%ylim([-2 1.5])
%         LogScale('xy',10)
%         title(UPDOWN{uu})
%         ylabel('Dur (s)');xlabel('mean EMG')
        
      
NiceSave('ReturnMap',figfolder,baseName)




end

