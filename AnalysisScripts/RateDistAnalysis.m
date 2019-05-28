function [  ] = RateDistAnalysis(basePath,figfolder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%
%basePath = '/mnt/proraidDL/Database/WMProbeData/180213_WT_M1M3_LFP_Layers_Pupil_EMG_Pole/180213_WT_M1M3_LFP_Layers_Pupil_EMG_180213_113045';
%basePath = '/home/dlevenstein/ProjectRepos/ACh-and-CorticalState/Dataset/180605_WT_M1M3_LFP_Layers_Pupil_EMG_180605_121846';
basePath = '/Users/dlevenstein/Desktop/171206_WT_EM1M3';
%basePath = '/Users/dlevenstein/Desktop/180704_KO_EM1M3';
%basePath = pwd;
figfolder = '/Users/dlevenstein/Project Repos/ACh-and-CorticalState/AnalysisScripts/AnalysisFigs/RateDistAnalysis';
%figfolder = '/home/dlevenstein/ProjectRepos/ACh-and-CorticalState/AnalysisScripts/AnalysisFigs/UPDOWNandPupilAnalysis';
%figfolder = fullfile(basePath,'AnalysisFigures');
%%
baseName = bz_BasenameFromBasepath(basePath);
recparms = bz_getSessionInfo(basePath,'noPrompts',true);

%% Detect Slow Waves
%CTXChans = recparms.SpkGrps.Channels(23:46);
spikes = bz_GetSpikes('basepath',basePath);

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

%%
WHNWH = {'Wh','NWh'};
HILO = {'lopup','hipup'};
%% Get rid of recording start/stop artifact and restrict to spontaneous behavior

%Restricting SPONT UP/DOWNs
load(fullfile(basePath,[baseName,'.MergePoints.events.mat']),'MergePoints');
sidx = find(startsWith(MergePoints.foldernames,"Spont"));
sponttimes = [MergePoints.timestamps(sidx(1),1) MergePoints.timestamps(sidx(end),2)];


maxtimejump = 1; %s
pupilcycle.amp = NanPadJumps( pupilcycle.timestamps,pupilcycle.amp,maxtimejump );
pupilcycle.phase = NanPadJumps( pupilcycle.timestamps,pupilcycle.phase,maxtimejump );
pupildilation.data = NanPadJumps( pupildilation.timestamps,pupildilation.data,maxtimejump );

%% Calculate rate mat

spikemat = bz_SpktToSpkmat(spikes,'dt',0.005,'binsize',0.025,'win',sponttimes);
spikemat.MUA = sum(spikemat.data,2);
spikemat.logMUA = log10(spikemat.MUA);
spikemat.logMUA(isinf(spikemat.logMUA)) = -1;
%%
spikemat.pupphase = interp1(pupilcycle.timestamps,pupilcycle.phase,...
        spikemat.timestamps,'nearest');
spikemat.amp = interp1(pupilcycle.timestamps,log10(pupilcycle.amp),...
    spikemat.timestamps,'nearest');
spikemat.pupsize = interp1(pupildilation.timestamps,log10(pupildilation.data),...
    spikemat.timestamps,'nearest');
spikemat.hipup = spikemat.amp>pupilcycle.pupthresh;
spikemat.lopup = ~spikemat.hipup;
spikemat.Wh = InIntervals(spikemat.timestamps,EMGwhisk.ints.Wh);
spikemat.NWh = InIntervals(spikemat.timestamps,EMGwhisk.ints.NWh);
spikemat.EMG = interp1(EMGwhisk.timestamps,EMGwhisk.EMGsm,...
    spikemat.timestamps,'nearest');

%%

maxMUA = 80;
mintime = 2.5; %s
for ww = 1:2
        [ConditionalRate.pup.(WHNWH{ww})] = ...
            ConditionalHist(spikemat.pupsize(spikemat.(WHNWH{ww})),spikemat.MUA(spikemat.(WHNWH{ww})),...
            'Ybounds',[0 maxMUA],'numYbins',maxMUA+1,...
            'Xbounds',[-0.25 0.25],'minX',mintime./spikemat.dt,'numXbins',30);
        
        [ConditionalRate.EMG] = ...
            ConditionalHist(log10(spikemat.EMG),spikemat.MUA,...
            'Ybounds',[0 maxMUA],'numYbins',maxMUA+1,...
            'Xbounds',[-1.5 0.75],'minX',mintime./spikemat.dt,'numXbins',30);
    for pp = 1:2
        [ConditionalRate.(HILO{pp}).(WHNWH{ww})] = ...
            ConditionalHist(spikemat.pupphase(spikemat.(WHNWH{ww})&spikemat.(HILO{pp})),...
            spikemat.MUA(spikemat.(WHNWH{ww})&spikemat.(HILO{pp})),...
            'Ybounds',[0 maxMUA],'numYbins',maxMUA+1,...
            'Xbounds',[-pi pi],'minX',mintime./spikemat.dt,'numXbins',30);
        
    end
end

% for ww = 1:2
%         [ConditionalRate.pup.(WHNWH{ww})] = ...
%             ConditionalHist(spikemat.pupsize(spikemat.(WHNWH{ww})),spikemat.logMUA(spikemat.(WHNWH{ww})),...
%             'Ybounds',[0 2.5],'numYbins',maxMUA+1,...
%             'Xbounds',[-0.25 0.25],'minX',15,'numXbins',30);
%     for pp = 1:2
%         [ConditionalRate.(HILO{pp}).(WHNWH{ww})] = ...
%             ConditionalHist(spikemat.pupphase(spikemat.(WHNWH{ww})&spikemat.(HILO{pp})),...
%             spikemat.logMUA(spikemat.(WHNWH{ww})&spikemat.(HILO{pp})),...
%             'Ybounds',[0 2.5],'numYbins',maxMUA+1,...
%             'Xbounds',[-pi pi],'minX',15,'numXbins',30);
%     end
% end


%%
cosx = linspace(-pi,pi,100);
cospamp = [2 20];

figure
for ww = 1:2
    subplot(3,3,(ww-1)*3+1)
    imagesc(ConditionalRate.pup.(WHNWH{ww}).Xbins,ConditionalRate.pup.(WHNWH{ww}).Ybins,...
        ConditionalRate.pup.(WHNWH{ww}).pYX')
    axis xy
    xlabel('Pupil (log)')
    ylabel(['Spike Count (',num2str(spikemat.binsize*1000),'ms) bins'])
    
    subplot(3,3,(ww-1)*3+2)
    for pp = 1:2
    imagesc(ConditionalRate.(HILO{pp}).(WHNWH{ww}).Xbins+2*pi*(pp-1),...
        ConditionalRate.(HILO{pp}).(WHNWH{ww}).Ybins,...
        ConditionalRate.(HILO{pp}).(WHNWH{ww}).pYX')
    hold on; axis xy; box off
    plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp)-1.5,'w')
    end
    xlim([-pi 3*pi])
    xlabel('Pupil Phase');
    ylabel(['Spike Count (',num2str(spikemat.binsize*1000),'ms) bins'])
end

	subplot(3,3,3)
    imagesc(ConditionalRate.EMG.Xbins,ConditionalRate.EMG.Ybins,...
        ConditionalRate.EMG.pYX')
    axis xy
    xlabel('EMG')
    ylabel(['Spike Count (',num2str(spikemat.binsize*1000),'ms) bins'])
    
NiceSave('RateDistByPupil',figfolder,baseName)





%% COnditional Duration
durrange = [-1.5 1.5]; %range for UP/DOWN durations - log
[ConditionalUPdurbyWhisk.allpup] = ...
    ConditionalHist(SlowWaves.timesinceWH,...
    log10(SlowWaves.durs.UP),...
    'Xbounds',[-1 4],'numYbins',75,'Ybounds',durrange,'minX',15,'numXbins',10);
for pp = 1:2
[ConditionalUPdurbyWhisk.(HILO{pp})] = ...
    ConditionalHist(SlowWaves.timesinceWH(SlowWaves.(HILO{pp}).UP),...
    log10(SlowWaves.durs.UP(SlowWaves.(HILO{pp}).UP)),...
    'Xbounds',[-1 4],'numYbins',75,'Ybounds',durrange,'minX',15,'numXbins',10);
end
%%
durrange = [-1.5 1.5]; %range for UP/DOWN durations - log
[ConditionalUPdurbyWhisk.allpup] = ...
    ConditionalHist(SlowWaves.timebeforeWH,...
    log10(SlowWaves.durs.UP),...
    'Xbounds',[-5 0],'numYbins',75,'Ybounds',durrange,'minX',15,'numXbins',10);
for pp = 1:2
[ConditionalUPdurbyWhisk.(HILO{pp})] = ...
    ConditionalHist(SlowWaves.timebeforeWH(SlowWaves.(HILO{pp}).UP),...
    log10(SlowWaves.durs.UP(SlowWaves.(HILO{pp}).UP)),...
    'Xbounds',[-5 0],'numYbins',75,'Ybounds',durrange,'minX',15,'numXbins',10);
end

%%
UDcmap = {makeColorMap([1 1 1],[0.8 0 0]),makeColorMap([1 1 1],[0 0 0.8])};

figure
subplot(2,2,1)
scatter(SlowWaves.timesinceWH,(SlowWaves.pupsize.UP),3,log10(SlowWaves.durs.UP))
xlim([-1 4])
ColorbarWithAxis([-1.25 1],'UP Duration (s)','location','east','labelloc','top')
LogScale('c',10)
xlabel('Time Since WhUP Offset')
%LogScale('y',10)
ylabel('Pupil Size (log)')

subplot(2,2,3)
    colormap(gca,UDcmap{1})
    imagesc(ConditionalUPdurbyWhisk.allpup.Xbins,...
        ConditionalUPdurbyWhisk.allpup.Ybins,...
        ConditionalUPdurbyWhisk.allpup.pYX')
    axis xy
    LogScale('y',10)
    xlabel('Time Since WhUP Offset');ylabel('UP druation')
    
for pp = 1:2
subplot(2,2,pp*2)
    colormap(gca,UDcmap{1})
    imagesc(ConditionalUPdurbyWhisk.(HILO{pp}).Xbins,...
        ConditionalUPdurbyWhisk.(HILO{pp}).Ybins,...
        ConditionalUPdurbyWhisk.(HILO{pp}).pYX')
    axis xy
    title(HILO{pp})
    LogScale('y',10)
    xlabel('Time Since WhUP Offset');ylabel('UP druation')
end
    

%%
figure
subplot(2,2,1)
scatter(SlowWaves.timebeforeWH,log10(SlowWaves.durNWH),2,log10(SlowWaves.durs.UP))
xlim([-5 0])
ColorbarWithAxis([-1.25 1],'UP Duration (s)','location','northoutside','labelloc','top')
LogScale('c',10)
xlabel('Time Before Wh Onset')
LogScale('y',10)
ylabel('NWh Dur (s)')

subplot(2,2,3)
    colormap(gca,UDcmap{1})
    imagesc(ConditionalUPdurbyWhisk.allpup.Xbins,...
        ConditionalUPdurbyWhisk.allpup.Ybins,...
        ConditionalUPdurbyWhisk.allpup.pYX')
    axis xy
    LogScale('y',10)
    xlabel('Time Before Wh Onset');ylabel('UP druation')
    
for pp = 1:2
subplot(2,2,pp*2)
    colormap(gca,UDcmap{1})
    imagesc(ConditionalUPdurbyWhisk.(HILO{pp}).Xbins,...
        ConditionalUPdurbyWhisk.(HILO{pp}).Ybins,...
        ConditionalUPdurbyWhisk.(HILO{pp}).pYX')
    axis xy
    title(HILO{pp})
    LogScale('y',10)
    xlabel('Time Before Wh Onset');ylabel('UP druation')
end
%% Histograms

durrange = [-1.5 1.5]; %range for UP/DOWN durations - log

for uu = 1:2
    for ww = 1:2
        
        %UD Duration by PSS
        [ConditionalUPDOWN.(UPDOWN{uu}).PSS.(WHNWH{ww})] = ...
            ConditionalHist((SlowWaves.PSS.(UPDOWN{uu})(SlowWaves.(WHNWH{ww}).(UPDOWN{uu}))),...
            log10(SlowWaves.durs.(UPDOWN{uu})(SlowWaves.(WHNWH{ww}).(UPDOWN{uu}))),...
            'Xbounds',[-3.25 -0.75],'numYbins',75,'Ybounds',durrange,'minX',15,'numXbins',15);
        
        %PSS Histogram by wh/nwh
        PSShist.winsize = PSS.winsize;
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
                'Xbounds',[-pi pi],'numYbins',75,'Ybounds',durrange,'minX',15,'numXbins',20);
            %U/D Duration Return maps split by pupil/whisk groups
            ReturnHist.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}) = hist3(...
                [log10(SlowWaves.durs.(UPDOWN{uu})(incondition)),...
                log10(SlowWaves.durs_np1.(UPDOWN{uu})(incondition))],...
                {ConditionalUPDOWN.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).Ybins,...
                ConditionalUPDOWN.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).Ybins});
            ReturnHist.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}) = ...
                ReturnHist.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww})./...
                sum(ReturnHist.(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww})(:));
            
            %UD Duration by PSS - split by hi/lo pup
            [ConditionalUPDOWN.(UPDOWN{uu}).PSS.(WHNWH{ww}).(HILO{pp})] = ...
                ConditionalHist((SlowWaves.PSS.(UPDOWN{uu})(incondition)),...
                log10(SlowWaves.durs.(UPDOWN{uu})(incondition)),...
                'Xbounds',[-3.25 -0.75],'numYbins',75,'Ybounds',durrange,'minX',15,'numXbins',15);
            
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
            sum(ReturnHist.(UPDOWN{uu}).allpup.(WHNWH{ww})(:));
        
        %UD Duation by pupil size
        [ConditionalUPDOWN.(UPDOWN{uu}).pupmag.(WHNWH{ww})] = ...
            ConditionalHist((SlowWaves.pupsize.(UPDOWN{uu})(SlowWaves.(WHNWH{ww}).(UPDOWN{uu}))),...
            log10(SlowWaves.durs.(UPDOWN{uu})(SlowWaves.(WHNWH{ww}).(UPDOWN{uu}))),...
            'Xbounds',[-0.2 0.2],'numYbins',75,'Ybounds',durrange,'minX',15,'numXbins',15);
       

        
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

