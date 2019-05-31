function [ SPECdepth,OSCdepth, speccorr] = LFPSpecbyDepthAnalysis(basePath,figfolder)
% Date XX/XX/20XX
%
%Question: 
%
%Plots
%-
%-
%
%% Load Header
%Initiate Paths
%reporoot = '/home/dlevenstein/ProjectRepos/ACh-and-CorticalState/';
%reporoot = '/Users/dlevenstein/Project Repos/ACh-and-CorticalState/';
%basePath = '/mnt/proraidDL/Database/WMData/AChPupil/171209_WT_EM1M3/';
%basePath = '/mnt/proraidDL/Database/WMData/AChPupil/180706_WT_EM1M3/';
%basePath = pwd;
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath,'noPrompts',true);

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
%Restricting SPONT UP/DOWNs
load(fullfile(basePath,[baseName,'.MergePoints.events.mat']),'MergePoints');
sidx = find(startsWith(MergePoints.foldernames,"Spont"));
sponttimes = [MergePoints.timestamps(sidx(1),1) MergePoints.timestamps(sidx(end),2)];

%% Get the depth info

depthinfo = rescaleCx(basePath);
inCTX = find(~isnan(depthinfo.ndepth));
CTXchans = depthinfo.channels(inCTX);
CTXdepth = -depthinfo.ndepth(inCTX);



%% Load only the cortex channels in the spontaneous time window
downsamplefactor = 2;
lfp = bz_GetLFP(CTXchans,...
    'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor,...
    'intervals',sponttimes);
lfp.chanlayers = depthinfo.layer(inCTX);
%% Calculate Spectrogram on all channels
clear spec
dt = 0.1;
spec.winsize = 1;
spec.frange = [0.5 256]; %Frequency lower than can be assessed for window because IRASA... but maybe this is bad for IRASA too
spec.nfreqs = 200;

noverlap = spec.winsize-dt;
spec.freqs = logspace(log10(spec.frange(1)),log10(spec.frange(2)),spec.nfreqs);
winsize_sf = round(spec.winsize .*lfp.samplingRate);
noverlap_sf = round(noverlap.*lfp.samplingRate);

spec.channels = lfp.channels;
for cc =1:length(spec.channels)
    cc
    [temp,~,spec.timestamps] = spectrogram(single(lfp.data(:,cc)),winsize_sf,noverlap_sf,spec.freqs,lfp.samplingRate);
    spec.data(:,:,cc) = log10(abs(temp))';
    spec.timestamps = spec.timestamps';
end

%% Correlation power-pupil,power-whisk




%% Take Mean Specgram by layer and calculate irasa

LAYERS = depthinfo.lnames;

for ll = 1:length(LAYERS)
    layerchans = strcmp(LAYERS{ll},lfp.chanlayers);
    spec.Layer(:,:,ll) = mean(spec.data(:,:,layerchans),3);
    [~,spec.osci(:,:,ll),spec.oscifreqs] = WaveIRASA(spec.Layer(:,:,ll),...
        'logamp',true,'freqs',spec.freqs);
end


%% Align Specgram by behavior
maxtimejump = 1; %s
pupilcycle.amp = NanPadJumps( pupilcycle.timestamps,pupilcycle.amp,maxtimejump );
pupilcycle.phase = NanPadJumps( pupilcycle.timestamps,pupilcycle.phase,maxtimejump );
pupildilation.data = NanPadJumps( pupildilation.timestamps,pupildilation.data,maxtimejump );
EMGwhisk.EMGsm = NanPadJumps( EMGwhisk.timestamps,EMGwhisk.EMGsm,maxtimejump );

spec.Wh = InIntervals(spec.timestamps,EMGwhisk.ints.Wh);
EMGwhisk.ints.ExpWh = bsxfun(@plus,EMGwhisk.ints.Wh,[-0.5 0.5]*spec.winsize);
spec.NWh = ~InIntervals(spec.timestamps,EMGwhisk.ints.ExpWh);
spec.hipup = interp1(pupilcycle.timestamps,single(pupilcycle.highpup),spec.timestamps,'nearest')==1;
spec.lopup = interp1(pupilcycle.timestamps,single(~pupilcycle.highpup),spec.timestamps,'nearest')==1;


spec.Wh = InIntervals(spec.timestamps,EMGwhisk.ints.Wh);
EMGwhisk.ints.ExpWh = bsxfun(@plus,EMGwhisk.ints.Wh,[-0.5 0.5]*spec.winsize);
spec.NWh = ~InIntervals(spec.timestamps,EMGwhisk.ints.ExpWh);
spec.hipup = interp1(pupilcycle.timestamps,single(pupilcycle.highpup),spec.timestamps,'nearest')==1;
spec.lopup = interp1(pupilcycle.timestamps,single(~pupilcycle.highpup),spec.timestamps,'nearest')==1;

spec.pupphase = interp1(pupilcycle.timestamps,pupilcycle.phase,spec.timestamps,'nearest');
spec.pup = interp1(pupildilation.timestamps,pupildilation.data,spec.timestamps,'nearest');
spec.EMG = interp1(EMGwhisk.timestamps,EMGwhisk.EMGsm,spec.timestamps,'nearest');

%% Whisk phase and duration

EMGwhisk.pupphase = interp1(pupilcycle.timestamps,pupilcycle.phase,EMGwhisk.ints.Wh(:,1),'nearest');
EMGwhisk.pupamp = interp1(pupilcycle.timestamps,pupilcycle.amp,EMGwhisk.ints.Wh(:,1),'nearest');
EMGwhisk.pup = interp1(pupildilation.timestamps,pupildilation.data,EMGwhisk.ints.Wh(:,1),'nearest');
EMGwhisk.dur = diff(EMGwhisk.ints.Wh,[],2);
EMGwhisk.hipup = log10(EMGwhisk.pupamp)>pupilcycle.pupthresh;
EMGwhisk.lopup = ~EMGwhisk.hipup;

%% Get relative time to nearest Wh onset
ONOFF = {'WhOn','WhOFF'};
WHNWH = {'Wh','NWh'};
HILO = {'lopup','hipup'};
LONGSHORT = {'long','short'};

%% Get Whisk On/Offset aligned spec, separated by long/short whisks

window = 5; %s
durthresh = 1; %
spec.whtime.WhOn = nan(size(spec.timestamps));
spec.whtime.WhOFF = nan(size(spec.timestamps));
spec.long.WhOn = false(size(spec.timestamps));
spec.long.WhOFF = false(size(spec.timestamps));
spec.short.WhOn = false(size(spec.timestamps));
spec.short.WhOFF = false(size(spec.timestamps));
for wh = 2:(size(EMGwhisk.ints.Wh,1)-1)
    %wh
    longwhisk = EMGwhisk.dur(wh)>durthresh;
    
    for oo = 1:2
        if oo == 1
            inwintimes = spec.timestamps > EMGwhisk.ints.ExpWh(wh-1,2) & ...
                spec.timestamps <= EMGwhisk.ints.Wh(wh,2);  
        elseif oo==2
            inwintimes = spec.timestamps > EMGwhisk.ints.Wh(wh,1) & ...
                spec.timestamps <= EMGwhisk.ints.ExpWh(wh+1,1); 
        end
        spec.whtime.(ONOFF{oo})(inwintimes) = spec.timestamps(inwintimes)-EMGwhisk.ints.Wh(wh,oo);
        spec.long.(ONOFF{oo})(inwintimes) = longwhisk;
        spec.short.(ONOFF{oo})(inwintimes) = ~longwhisk;
    end
    
end

%% Correlation with pupil
 speccorr.freqs = spec.freqs;
 speccorr.channels = spec.channels;
for cc = 1:length(spec.channels)
    specvarcorr = bz_LFPSpecToExternalVar( spec.data(:,:,cc),...
        log10(spec.pup),'specparms','input',...
        'figparms',true,'numvarbins',20,'varlim',[-0.25 0.25]);
    speccorr.pup(:,cc)  = specvarcorr.corr;
    
    specvarcorr = bz_LFPSpecToExternalVar( spec.data(:,:,cc),...
        log10(spec.EMG),'specparms','input',...
        'figparms',true,'numvarbins',20,'varlim',[-0.25 0.25]);
    speccorr.EMG(:,cc)  = specvarcorr.corr;
    
    close all
end
close all

% Interpolate PSS to normalized depth
speccorr.interpdepth = linspace(-1,0,100);
speccorr.pupinterp = interp1(CTXdepth',speccorr.pup',speccorr.interpdepth')';
speccorr.EMGinterp = interp1(CTXdepth',speccorr.EMG',speccorr.interpdepth')';
%%
figure

subplot(2,2,1)
imagesc(log10(speccorr.freqs),speccorr.interpdepth,speccorr.pupinterp')
hold on
plot(log10(speccorr.freqs([1 end])),-depthinfo.boundaries'*[1 1],'k')
crameri vik
LogScale('x',10)
axis xy
ColorbarWithAxis([-0.4 0.4],'Pupil Corr.')
%clim([-0.35 0.35])
xlabel('f (Hz)');ylabel('Depth')

subplot(2,2,2)
imagesc(log10(speccorr.freqs),speccorr.interpdepth,speccorr.EMGinterp')
hold on
plot(log10(speccorr.freqs([1 end])),-depthinfo.boundaries'*[1 1],'k')
LogScale('x',10)
axis xy
ColorbarWithAxis([-0.4 0.4],'EMG Corr.')
crameri vik
xlabel('f (Hz)');ylabel('Depth')

NiceSave('DepthFreqCorr',figfolder,baseName)
%%  Mean layer spec by pupil size, phase and whisking
%prepare for LFPspec....
SPECdepth.freqs = spec.freqs;
for dd = 1:length(LAYERS)
for ww = 1:2
[ ~,SPECdepth.(LAYERS{dd}).pup.(WHNWH{ww}) ] = bz_LFPSpecToExternalVar( spec.Layer(spec.(WHNWH{ww}),:,dd),...
    log10(spec.pup(spec.(WHNWH{ww}),:)),'specparms','input',...
    'figparms',true,'numvarbins',20,'varlim',[-0.25 0.25]);
    [~,SPECdepth.(LAYERS{dd}).pup.(WHNWH{ww}).mean_osc,SPECdepth.oscfreqs] = ...
        WaveIRASA(SPECdepth.(LAYERS{dd}).pup.(WHNWH{ww}).mean','logamp',true,'freqs',spec.freqs);


    for pp= 1:2
    [ ~,SPECdepth.(LAYERS{dd}).(HILO{pp}).(WHNWH{ww}) ] = bz_LFPSpecToExternalVar(...
        spec.Layer(spec.(WHNWH{ww})&spec.(HILO{pp}),:,dd),...
        spec.pupphase(spec.(WHNWH{ww})&spec.(HILO{pp})),'specparms','input',...
        'figparms',true,'numvarbins',20,'varlim',[-pi pi]);
    [~,SPECdepth.(LAYERS{dd}).(HILO{pp}).(WHNWH{ww}).mean_osc,SPECdepth.oscfreqs] = ...
        WaveIRASA(SPECdepth.(LAYERS{dd}).(HILO{pp}).(WHNWH{ww}).mean','logamp',true,'freqs',spec.freqs);

    end
end

for oo = 1:2
    for ll = 1:2
        [ ~,SPECdepth.(LAYERS{dd}).(ONOFF{oo}).(LONGSHORT{ll}) ] = bz_LFPSpecToExternalVar(...
            spec.Layer(spec.(LONGSHORT{ll}).(ONOFF{oo}),:,dd),...
            spec.whtime.(ONOFF{oo})(spec.(LONGSHORT{ll}).(ONOFF{oo}),:),'specparms','input',...
            'figparms',true,'numvarbins',40,'varlim',[-window window]);
        [~,SPECdepth.(LAYERS{dd}).(ONOFF{oo}).(LONGSHORT{ll}).mean_osc,SPECdepth.oscfreqs] = ...
            WaveIRASA(SPECdepth.(LAYERS{dd}).(ONOFF{oo}).(LONGSHORT{ll}).mean','logamp',true,'freqs',spec.freqs);
    end
    
    [ ~,SPECdepth.(LAYERS{dd}).(ONOFF{oo}).all ] = bz_LFPSpecToExternalVar(...
        spec.Layer(:,:,dd),...
        spec.whtime.(ONOFF{oo}),'specparms','input',...
        'figparms',true,'numvarbins',40,'varlim',[-window window]);
    [~,SPECdepth.(LAYERS{dd}).(ONOFF{oo}).all.mean_osc,SPECdepth.oscfreqs] = ...
        WaveIRASA(SPECdepth.(LAYERS{dd}).(ONOFF{oo}).all.mean','logamp',true,'freqs',spec.freqs);
end

    
[ ~,SPECdepth.(LAYERS{dd}).EMG ] = bz_LFPSpecToExternalVar( spec.Layer(:,:,dd),...
    log10(spec.EMG),'specparms','input',...
    'figparms',true,'numvarbins',40,'varlim',[-1.7 0.9]);
[~,SPECdepth.(LAYERS{dd}).EMG.mean_osc,SPECdepth.oscfreqs] = ...
    WaveIRASA(SPECdepth.(LAYERS{dd}).EMG.mean','logamp',true,'freqs',spec.freqs);
end


  %%
  
  cosx = linspace(-pi,pi,100);
cospamp = [0.025 0.3]*2;




figure
for dd = 1:6
for ww = 1:2
    subplot(6,4,(dd-1)*4+ww)
        for pp = 1:2
        imagesc( SPECdepth.(LAYERS{dd}).(HILO{pp}).(WHNWH{ww}).varbins+2*pi*(pp-1),...
            log10(SPECdepth.freqs),...
            SPECdepth.(LAYERS{dd}).(HILO{pp}).(WHNWH{ww}).mean)
        hold on; axis xy; box off
        plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp),'w')
        end   
        LogScale('y',10)
        
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        clim([3 5])
        xlim([-pi 3*pi])
        if dd == 6
        xlabel('Pupil Phase');
        end
        if ww == 1
            ylabel({LAYERS{dd},'Freq'})
        end
        if dd == 1
        title((WHNWH{ww}))
        end

        
        
    subplot(6,4,(dd-1)*4+ww+2)
        imagesc( SPECdepth.(LAYERS{dd}).pup.(WHNWH{ww}).varbins,...
            log10(SPECdepth.freqs),...
            SPECdepth.(LAYERS{dd}).pup.(WHNWH{ww}).mean)
        hold on; axis xy; box off
        LogScale('y',10)
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        clim([3 5])
        if dd == 6
        xlabel('Pupil Size');
        end
        if ww == 1
            ylabel('Freq')
        end
        if dd == 1
        title((WHNWH{ww}))
        end
        
       
        
end
end

NiceSave('DepthSPECandPup',figfolder,baseName)

%%

figure
for dd = 1:6
for oo = 1:2

subplot(6,4,(dd-1)*4+oo)
        imagesc( SPECdepth.(LAYERS{dd}).(ONOFF{oo}).all.varbins,...
            log10(SPECdepth.freqs),...
            SPECdepth.(LAYERS{dd}).(ONOFF{oo}).all.mean)
        hold on; axis xy; box off
        plot([0 0],[0 max(SPECdepth.freqs)],'w')
        clim([3 5])
        LogScale('y',10)
        if dd == 6
        xlabel(['t - ',(ONOFF{oo})])
        end
        if oo == 1
            ylabel({LAYERS{dd},'Freq'})
        end
        if dd == 1
        title((ONOFF{oo}))
        end
end
 
subplot(6,4,(dd-1)*4+4)
        imagesc( SPECdepth.(LAYERS{dd}).EMG.varbins,...
            log10(SPECdepth.freqs),...
            SPECdepth.(LAYERS{dd}).EMG.mean)
        hold on; axis xy; box off
        plot(SPECdepth.(LAYERS{dd}).EMG.varbins,SPECdepth.(LAYERS{dd}).EMG.vardist*1000,'w')
        plot(log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],[0 max(SPECdepth.freqs)],'k--')
        clim([3 5])
        ylabel('Freq');
        LogScale('y',10)
                if dd == 6
        xlabel('EMG')
        end
        
end
NiceSave('DepthSPECandWhisk',figfolder,baseName)


  %%
  
  cosx = linspace(-pi,pi,100);
cospamp = [0.025 0.3]*2;




figure
for dd = 1:6
for ww = 1:2
    subplot(6,4,(dd-1)*4+ww)
        for pp = 1:2
        imagesc( SPECdepth.(LAYERS{dd}).(HILO{pp}).(WHNWH{ww}).varbins+2*pi*(pp-1),...
            log10(SPECdepth.oscfreqs),...
            SPECdepth.(LAYERS{dd}).(HILO{pp}).(WHNWH{ww}).mean_osc')
        hold on; axis xy; box off
        plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp),'w')
        end   
        LogScale('y',10)
        
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        clim([-0.2 0.2])
        xlim([-pi 3*pi])
        if dd == 6
        xlabel('Pupil Phase');
        end
        if ww == 1
            ylabel({LAYERS{dd},'Freq'})
        end
        if dd == 1
        title((WHNWH{ww}))
        end

        
        
    subplot(6,4,(dd-1)*4+ww+2)
        imagesc( SPECdepth.(LAYERS{dd}).pup.(WHNWH{ww}).varbins,...
            log10(SPECdepth.oscfreqs),...
            SPECdepth.(LAYERS{dd}).pup.(WHNWH{ww}).mean_osc')
        hold on; axis xy; box off
        LogScale('y',10)
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        clim([-0.2 0.2])
        if dd == 6
        xlabel('Pupil Size');
        end
        if ww == 1
            ylabel('Freq')
        end
        if dd == 1
        title((WHNWH{ww}))
        end
        
       
        
end
end

NiceSave('DepthSPECremandPup',figfolder,baseName)

%%

figure
for dd = 1:6
for oo = 1:2

subplot(6,4,(dd-1)*4+oo)
        imagesc( SPECdepth.(LAYERS{dd}).(ONOFF{oo}).all.varbins,...
            log10(SPECdepth.oscfreqs),...
            SPECdepth.(LAYERS{dd}).(ONOFF{oo}).all.mean_osc')
        hold on; axis xy; box off
        plot([0 0],[0 max(SPECdepth.oscfreqs)],'w')
        clim([-0.2 0.2])
        LogScale('y',10)
        if dd == 6
        xlabel(['t - ',(ONOFF{oo})])
        end
        if oo == 1
            ylabel({LAYERS{dd},'Freq'})
        end
        if dd == 1
        title((ONOFF{oo}))
        end
end
 
subplot(6,4,(dd-1)*4+4)
        imagesc( SPECdepth.(LAYERS{dd}).EMG.varbins,...
            log10(SPECdepth.oscfreqs),...
            SPECdepth.(LAYERS{dd}).EMG.mean_osc')
        hold on; axis xy; box off
        plot(SPECdepth.(LAYERS{dd}).EMG.varbins,SPECdepth.(LAYERS{dd}).EMG.vardist*1000,'w')
        plot(log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],[0 max(SPECdepth.oscfreqs)],'k--')
        clim([-0.2 0.2])
        LogScale('y',10)
        ylabel('Freq');
                if dd == 6
        xlabel('EMG')
        end
        
end
NiceSave('DepthSPECremandWhisk',figfolder,baseName)


%%  Mean depth spec by pupil size, phase and whisking
%prepare for LFPspec....
OSCdepth.freqs = spec.oscifreqs;
for dd = 1:length(LAYERS)
for ww = 1:2
[ ~,OSCdepth.(LAYERS{dd}).pup.(WHNWH{ww}) ] = bz_LFPSpecToExternalVar( spec.osci(spec.(WHNWH{ww}),:,dd),...
    log10(spec.pup(spec.(WHNWH{ww}),:)),'specparms','input',...
    'figparms',true,'numvarbins',20,'varlim',[-0.25 0.25]);

    for pp= 1:2
    [ ~,OSCdepth.(LAYERS{dd}).(HILO{pp}).(WHNWH{ww}) ] = bz_LFPSpecToExternalVar(...
        spec.osci(spec.(WHNWH{ww})&spec.(HILO{pp}),:,dd),...
        spec.pupphase(spec.(WHNWH{ww})&spec.(HILO{pp})),'specparms','input',...
        'figparms',true,'numvarbins',20,'varlim',[-pi pi]);

    end
end

for oo = 1:2
    for ll = 1:2
        [ ~,OSCdepth.(LAYERS{dd}).(ONOFF{oo}).(LONGSHORT{ll}) ] = bz_LFPSpecToExternalVar(...
            spec.osci(spec.(LONGSHORT{ll}).(ONOFF{oo}),:,dd),...
            spec.whtime.(ONOFF{oo})(spec.(LONGSHORT{ll}).(ONOFF{oo}),:),'specparms','input',...
            'figparms',true,'numvarbins',40,'varlim',[-window window]);
    end
    
    [ ~,OSCdepth.(LAYERS{dd}).(ONOFF{oo}).all ] = bz_LFPSpecToExternalVar(...
        spec.osci(:,:,dd),...
        spec.whtime.(ONOFF{oo}),'specparms','input',...
        'figparms',true,'numvarbins',40,'varlim',[-window window]);
end

    
[ ~,OSCdepth.(LAYERS{dd}).EMG ] = bz_LFPSpecToExternalVar( spec.osci(:,:,dd),...
    log10(spec.EMG),'specparms','input',...
    'figparms',true,'numvarbins',40,'varlim',[-1.7 0.9]);
end

  %%
  
  cosx = linspace(-pi,pi,100);
cospamp = [0.025 0.3]*2;


figure
for dd = 1:6
for ww = 1:2
    subplot(6,4,(dd-1)*4+ww)
        for pp = 1:2
        imagesc( OSCdepth.(LAYERS{dd}).(HILO{pp}).(WHNWH{ww}).varbins+2*pi*(pp-1),...
            log10(OSCdepth.freqs),...
            OSCdepth.(LAYERS{dd}).(HILO{pp}).(WHNWH{ww}).mean)
        hold on; axis xy; box off
        plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp),'w')
        end   
        LogScale('y',10)
        
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        clim([-0.2 0.2])
        %colorbar
        xlim([-pi 3*pi])
        if dd == 6
        xlabel('Pupil Phase');
        end
        if ww == 1
            ylabel({LAYERS{dd},'Freq'})
        end
        if dd == 1
        title((WHNWH{ww}))
        end

        
        
    subplot(6,4,(dd-1)*4+ww+2)
        imagesc( OSCdepth.(LAYERS{dd}).pup.(WHNWH{ww}).varbins,...
            log10(OSCdepth.freqs),...
            OSCdepth.(LAYERS{dd}).pup.(WHNWH{ww}).mean)
        hold on; axis xy; box off
        LogScale('y',10)
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        clim([-0.2 0.2])
        if dd == 6
        xlabel('Pupil Size');
        end
        if ww == 1
            ylabel('Freq')
        end
        if dd == 1
        title((WHNWH{ww}))
        end
        
       
        
end
end

NiceSave('DepthOSCandPup',figfolder,baseName)

%%

figure
for dd = 1:6
for oo = 1:2

subplot(6,4,(dd-1)*4+oo)
        imagesc( OSCdepth.(LAYERS{dd}).(ONOFF{oo}).all.varbins,...
            log10(OSCdepth.freqs),...
            OSCdepth.(LAYERS{dd}).(ONOFF{oo}).all.mean)
        hold on; axis xy; box off
        plot([0 0],[0 max(OSCdepth.freqs)],'w')
        clim([-0.2 0.2])
         LogScale('y',10)
        if dd == 6
        xlabel(['t - ',(ONOFF{oo})])
        end
        if oo == 1
            ylabel({LAYERS{dd},'Freq'})
        end
        if dd == 1
        title((ONOFF{oo}))
        end
end
 
subplot(6,4,(dd-1)*4+4)
        imagesc( OSCdepth.(LAYERS{dd}).EMG.varbins,...
            log10(OSCdepth.freqs),...
            OSCdepth.(LAYERS{dd}).EMG.mean)
        hold on; axis xy; box off
        plot(OSCdepth.(LAYERS{dd}).EMG.varbins,OSCdepth.(LAYERS{dd}).EMG.vardist*1000,'w')
        plot(log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],[0 max(OSCdepth.freqs)],'k--')
        clim([-0.2 0.2])
        LogScale('y',10)
        ylabel('Freq');
                if dd == 6
        xlabel('EMG')
        end
        
end
NiceSave('DepthOSCandWhisk',figfolder,baseName)



