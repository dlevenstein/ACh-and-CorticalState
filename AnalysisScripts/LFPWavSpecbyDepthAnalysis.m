function [ SPECdepth,OSCdepth] = LFPWavSpecbyDepthAnalysis(basePath,figfolder)
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
%pupildilation = bz_LoadBehavior(basePath,'pupildiameter');
[ pupilcycle ] = ExtractPupilCycle( basePath );


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
% %dt = 0.1;
% spec.winsize = 1;
 spec.frange = [2 128]; %Frequency lower than can be assessed for window because IRASA... but maybe this is bad for IRASA too
 spec.nfreqs = 100;
% 
% % noverlap = spec.winsize-dt;
% % spec.freqs = logspace(log10(spec.frange(1)),log10(spec.frange(2)),spec.nfreqs);
% % winsize_sf = round(spec.winsize .*lfp.samplingRate);
% % noverlap_sf = round(noverlap.*lfp.samplingRate);
% 
% spec.channels = lfp.channels;
spec.channels = CTXchans;
for cc =1:length(spec.channels)
    bz_Counter(cc,length(spec.channels),'Channel')
    %FFT
    %[temp,~,spec.timestamps] = spectrogram(single(lfp.data(:,cc)),winsize_sf,noverlap_sf,spec.freqs,lfp.samplingRate);
%     spec.data(:,:,cc) = log10(abs(temp))';
%     spec.timestamps = spec.timestamps';

    %Wavelets - recalculate
%     [wavespec] = bz_WaveSpec(lfp,'frange',spec.frange,'nfreqs',spec.nfreqs,...
%         'chanID',lfp.channels(cc),'ncyc',12); 
%     spec.data(:,:,cc) = log10(abs(downsample(wavespec.data,5)));
%     spec.timestamps = downsample(wavespec.timestamps,5);
%     spec.freqs = wavespec.freqs; 

    %Loaded Wavelets
%     load(fullfile(basePath,'WaveSpec_Downsampled',[baseName,'.',num2str(spec.channels(cc)),'.WaveSpec.lfp.mat']));
%     inspont = InIntervals(wavespec.timestamps,sponttimes);
%     spec.data(:,:,cc) = log10(abs(wavespec.data(inspont,:)));
%     spec.timestamps = wavespec.timestamps(inspont,:);
%     spec.freqs = wavespec.freqs; 
%%
    [specslope,~] = bz_PowerSpectrumSlope(lfp,10,0.005,'channels',spec.channels(cc),...
        'frange',spec.frange,'spectype','wavelet','nfreqs',spec.nfreqs,'ints',sponttimes);
    spec.data(:,:,cc) = specslope.specgram;
    spec.osci(:,:,cc) = specslope.resid;
    spec.PSS(:,cc) = specslope.data;
    spec.timestamps = specslope.timestamps;
    spec.freqs = specslope.freqs; 
    clear specslope
end
%%
spec.winsize = 1;
spec.chanlayers = depthinfo.layer(inCTX);
spec.osci(spec.osci<0) = 0;
%% Take Mean Specgram by layer and calculate irasa

LAYERS = depthinfo.lnames;

for ll = 1:length(LAYERS)
    layerchans = strcmp(LAYERS{ll},spec.chanlayers);
    repchan(ll) = round(median(find(layerchans)));
    spec.Layer(:,:,ll) = mean(spec.data(:,:,layerchans),3);
    spec.LayerOsci(:,:,ll) = mean(spec.osci(:,:,layerchans),3);
    spec.LayerPSS(:,ll) = mean(spec.PSS(:,layerchans),2);
    % Median Normalize Spectrogram
    spec.Layer(:,:,ll) = NormToInt(spec.Layer(:,:,ll),'median');
end
spec = rmfield(spec,'data');
spec = rmfield(spec,'osci');
%%

%%
xwin = bz_RandomWindowInIntervals(spec.timestamps([1 end]),10);
figure
subplot(2,1,1)
imagesc(spec.timestamps,log2(spec.freqs),spec.LayerOsci(:,:,4)')
hold on
plot(lfp.timestamps,bz_NormToRange(single(lfp.data(:,repchan(4)))),'w')
axis xy
colorbar
caxis([0 0.7])
LogScale('y',2)
xlim(xwin)
%% Get PSS by layer
% load(fullfile(basePath,[baseName,'.PowerSpectrumSlope.lfp.mat']));
% 
% [~,~,PSS.CTXchans] = intersect(CTXchans,PSpecSlope.Shortwin.channels,'stable');
% PSS.data = PSpecSlope.Shortwin.PSS(:,PSS.CTXchans);
% PSS.timestamps = PSpecSlope.Shortwin.timestamps;
% PSS.samplingRate = 1/mean(diff(PSS.timestamps));
% PSS.winsize = PSpecSlope.Shortwin.movingwin(1);
% PSS.depth = CTXdepth;
% 
% PSS.osci = PSpecSlope.Shortwin.OSCI(:,:,PSS.CTXchans);
% PSS.osci = shiftdim(PSS.osci,1);
% PSS.freqs = PSpecSlope.Shortwin.freqs;
% PSS.chanlayers = depthinfo.layer(inCTX);
% 
% inspont = InIntervals(PSS.timestamps,sponttimes);
% PSS.timestamps = PSS.timestamps(inspont);
% PSS.data = PSS.data(inspont,:);
% 
% for ll = 1:length(LAYERS)
%     PSS.Lchans.(LAYERS{ll}) = strcmp(LAYERS{ll},PSS.chanlayers);
%     PSS.Layer(:,ll) = (mean(PSS.data(:,PSS.Lchans.(LAYERS{ll})),2));
%     
%     %PSS at spec times
%     spec.LayerPSS(:,ll) = interp1(PSS.timestamps, PSS.Layer(:,ll),spec.timestamps,'nearest');
% end
% 
% clear PSpecSlope

%% Align Specgram by behavior
maxtimejump = 1; %s
pupilcycle.amp = NanPadJumps( pupilcycle.timestamps,pupilcycle.amp,maxtimejump );
pupilcycle.phase = NanPadJumps( pupilcycle.timestamps,pupilcycle.phase,maxtimejump );
pupilcycle.data = NanPadJumps( pupilcycle.timestamps,pupilcycle.data,maxtimejump );
EMGwhisk.EMGsm = NanPadJumps( EMGwhisk.timestamps,EMGwhisk.EMGsm,maxtimejump );

spec.Wh = InIntervals(spec.timestamps,EMGwhisk.ints.Wh);
EMGwhisk.ints.ExpWh = bsxfun(@plus,EMGwhisk.ints.Wh,[-0.5 0.5]*spec.winsize);
spec.NWh = ~InIntervals(spec.timestamps,EMGwhisk.ints.ExpWh);
spec.hipup = interp1(pupilcycle.timestamps,single(pupilcycle.states==2),spec.timestamps,'nearest')==1;
spec.lopup = interp1(pupilcycle.timestamps,single(pupilcycle.states==1),spec.timestamps,'nearest')==1;


spec.Wh = InIntervals(spec.timestamps,EMGwhisk.ints.Wh);
EMGwhisk.ints.ExpWh = bsxfun(@plus,EMGwhisk.ints.Wh,[-0.5 0.5]*spec.winsize);
spec.NWh = ~InIntervals(spec.timestamps,EMGwhisk.ints.ExpWh);
spec.hipup = interp1(pupilcycle.timestamps,single(pupilcycle.states==2),spec.timestamps,'nearest')==1;
spec.lopup = interp1(pupilcycle.timestamps,single(pupilcycle.states==1),spec.timestamps,'nearest')==1;

spec.pupphase = interp1(pupilcycle.timestamps,pupilcycle.phase,spec.timestamps,'nearest');
spec.pup = interp1(pupilcycle.timestamps,pupilcycle.data,spec.timestamps,'nearest');
spec.EMG = interp1(EMGwhisk.timestamps,EMGwhisk.EMGsm,spec.timestamps,'nearest');


%% Whisk phase and duration

EMGwhisk.pupphase = interp1(pupilcycle.timestamps,pupilcycle.phase,EMGwhisk.ints.Wh(:,1),'nearest');
EMGwhisk.pupamp = interp1(pupilcycle.timestamps,pupilcycle.amp,EMGwhisk.ints.Wh(:,1),'nearest');
EMGwhisk.pup = interp1(pupilcycle.timestamps,pupilcycle.data,EMGwhisk.ints.Wh(:,1),'nearest');
EMGwhisk.dur = diff(EMGwhisk.ints.Wh,[],2);
EMGwhisk.hipup = log10(EMGwhisk.pupamp)>pupilcycle.detectionparms.pupthresh;
EMGwhisk.lopup = ~EMGwhisk.hipup;

%% Get relative time to nearest Wh onset
ONOFF = {'WhOn','WhOFF'};
WHNWH = {'Wh','NWh'};
HILO = {'lopup','hipup'};
LONGSHORT = {'long','short'};

%% Get Whisk On/Offset aligned spec, separated by long/short whisks

window = 2; %s
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
%  speccorr.freqs = spec.freqs;
%  speccorr.channels = spec.channels;
% for cc = 1:length(spec.channels)
%     specvarcorr = bz_LFPSpecToExternalVar( spec.data(:,:,cc),...
%         log10(spec.pup),'specparms','input',...
%         'figparms',[],'numvarbins',20,'varlim',[-0.25 0.25]);
%     speccorr.pup(:,cc)  = specvarcorr.corr;
%     
%     specvarcorr = bz_LFPSpecToExternalVar( spec.data(:,:,cc),...
%         log10(spec.EMG),'specparms','input',...
%         'figparms',[],'numvarbins',20,'varlim',[-0.25 0.25]);
%     speccorr.EMG(:,cc)  = specvarcorr.corr;
%     
%     
% end
% %close all
% spec = rmfield(spec,'data');
% % Interpolate PSS to normalized depth
% speccorr.interpdepth = linspace(-1,0,100);
% speccorr.pupinterp = interp1(CTXdepth',speccorr.pup',speccorr.interpdepth')';
% speccorr.EMGinterp = interp1(CTXdepth',speccorr.EMG',speccorr.interpdepth')';
% %% Figure: Correlation
% figure
% 
% subplot(2,2,1)
% imagesc(log10(speccorr.freqs),speccorr.interpdepth,speccorr.pupinterp')
% hold on
% plot(log10(speccorr.freqs([1 end])),-depthinfo.boundaries'*[1 1],'k')
% crameri vik
% LogScale('x',10)
% axis xy
% ColorbarWithAxis([-0.4 0.4],'Pupil Corr.')
% %clim([-0.35 0.35])
% xlabel('f (Hz)');ylabel('Depth')
% 
% subplot(2,2,2)
% imagesc(log10(speccorr.freqs),speccorr.interpdepth,speccorr.EMGinterp')
% hold on
% plot(log10(speccorr.freqs([1 end])),-depthinfo.boundaries'*[1 1],'k')
% LogScale('x',10)
% axis xy
% ColorbarWithAxis([-0.4 0.4],'EMG Corr.')
% crameri vik
% xlabel('f (Hz)');ylabel('Depth')
% 
% NiceSave('DepthFreqCorr',figfolder,baseName)

%%  Mean depth spec by pupil size, phase and whisking
%prepare for LFPspec....
SPECdepth.freqs = spec.freqs;
for dd = 1:length(LAYERS)
for ww = 1:2
[ ~,SPECdepth.(LAYERS{dd}).pup.(WHNWH{ww}) ] = bz_LFPSpecToExternalVar( spec.Layer(spec.(WHNWH{ww}),:,dd),...
    log10(spec.pup(spec.(WHNWH{ww}),:)),'specparms','input',...
    'figparms',[],'numvarbins',20,'varlim',[-0.25 0.25],'minX',500);
%     [~,SPECdepth.(LAYERS{dd}).pup.(WHNWH{ww}).mean_osc,SPECdepth.oscfreqs] = ...
%         WaveIRASA(SPECdepth.(LAYERS{dd}).pup.(WHNWH{ww}).mean','logamp',true,'freqs',spec.freqs);


    for pp= 1:2
    [ ~,SPECdepth.(LAYERS{dd}).(HILO{pp}).(WHNWH{ww}) ] = bz_LFPSpecToExternalVar(...
        spec.Layer(spec.(WHNWH{ww})&spec.(HILO{pp}),:,dd),...
        spec.pupphase(spec.(WHNWH{ww})&spec.(HILO{pp})),'specparms','input',...
        'figparms',[],'numvarbins',25,'varlim',[-pi pi],'minX',500);
%     [~,SPECdepth.(LAYERS{dd}).(HILO{pp}).(WHNWH{ww}).mean_osc,SPECdepth.oscfreqs] = ...
%         WaveIRASA(SPECdepth.(LAYERS{dd}).(HILO{pp}).(WHNWH{ww}).mean','logamp',true,'freqs',spec.freqs);

    end
end

for oo = 1:2
    for ll = 1:2
        [ ~,SPECdepth.(LAYERS{dd}).(ONOFF{oo}).(LONGSHORT{ll}) ] = bz_LFPSpecToExternalVar(...
            spec.Layer(spec.(LONGSHORT{ll}).(ONOFF{oo}),:,dd),...
            spec.whtime.(ONOFF{oo})(spec.(LONGSHORT{ll}).(ONOFF{oo}),:),'specparms','input',...
            'figparms',[],'numvarbins',80,'varlim',[-window window],'minX',500);
%         [~,SPECdepth.(LAYERS{dd}).(ONOFF{oo}).(LONGSHORT{ll}).mean_osc,SPECdepth.oscfreqs] = ...
%             WaveIRASA(SPECdepth.(LAYERS{dd}).(ONOFF{oo}).(LONGSHORT{ll}).mean','logamp',true,'freqs',spec.freqs);
    end
    
    [ ~,SPECdepth.(LAYERS{dd}).(ONOFF{oo}).all ] = bz_LFPSpecToExternalVar(...
        spec.Layer(:,:,dd),...
        spec.whtime.(ONOFF{oo}),'specparms','input',...
        'figparms',[],'numvarbins',400,'varlim',[-window window],'minX',500);
%     [~,SPECdepth.(LAYERS{dd}).(ONOFF{oo}).all.mean_osc,SPECdepth.oscfreqs] = ...
%         WaveIRASA(SPECdepth.(LAYERS{dd}).(ONOFF{oo}).all.mean','logamp',true,'freqs',spec.freqs);
end

    
[ ~,SPECdepth.(LAYERS{dd}).EMG ] = bz_LFPSpecToExternalVar( spec.Layer(:,:,dd),...
    log10(spec.EMG),'specparms','input',...
    'figparms',[],'numvarbins',40,'varlim',[-1.7 0.9],'minX',500);
% [~,SPECdepth.(LAYERS{dd}).EMG.mean_osc,SPECdepth.oscfreqs] = ...
%     WaveIRASA(SPECdepth.(LAYERS{dd}).EMG.mean','logamp',true,'freqs',spec.freqs);

[ ~,SPECdepth.(LAYERS{dd}).PSS ] = bz_LFPSpecToExternalVar( spec.Layer(:,:,dd),...
    spec.LayerPSS(:,dd),'specparms','input',...
    'figparms',[],'numvarbins',40,'varlim',[-3 -0.5],'minX',500);
end


  %%
  
  cosx = linspace(-pi,pi,100);
cospamp = [0.025 0.3]*2;


speclim = [2 4];
speclim = [0.9 1.1];

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
        ylim([0 2.5])
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        clim(speclim)
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
        imagesc( SPECdepth.(LAYERS{dd}).pup.(WHNWH{ww}).varbins,...
            log10(SPECdepth.freqs),...
            SPECdepth.(LAYERS{dd}).pup.(WHNWH{ww}).mean)
        hold on; axis xy; box off
        LogScale('y',10)
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        clim(speclim)
        ylim([0 2.5])
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
        clim(speclim)
        ylim([0 2.5])
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
        %plot(log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],[0 max(SPECdepth.freqs)],'k--')
        clim(speclim)
        ylim([0 2.5])
        ylabel('Freq');
        LogScale('y',10)
                if dd == 6
        xlabel('EMG')
        end
        
end
NiceSave('DepthSPECandWhisk',figfolder,baseName)

%%
%   %%
%   
%   cosx = linspace(-pi,pi,100);
% cospamp = [0.025 0.3]*2;
% 
% 
% 
% 
% figure
% for dd = 1:6
% for ww = 1:2
%     subplot(6,4,(dd-1)*4+ww)
%         for pp = 1:2
%         imagesc( SPECdepth.(LAYERS{dd}).(HILO{pp}).(WHNWH{ww}).varbins+2*pi*(pp-1),...
%             log10(SPECdepth.oscfreqs),...
%             SPECdepth.(LAYERS{dd}).(HILO{pp}).(WHNWH{ww}).mean_osc')
%         hold on; axis xy; box off
%         plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp),'w')
%         end   
%         LogScale('y',10)
%         
%         %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
%         clim([-0.2 0.2])
%         xlim([-pi 3*pi])
%         if dd == 6
%         xlabel('Pupil Phase');
%         end
%         if ww == 1
%             ylabel({LAYERS{dd},'Freq'})
%         end
%         if dd == 1
%         title((WHNWH{ww}))
%         end
% 
%         
%         
%     subplot(6,4,(dd-1)*4+ww+2)
%         imagesc( SPECdepth.(LAYERS{dd}).pup.(WHNWH{ww}).varbins,...
%             log10(SPECdepth.oscfreqs),...
%             SPECdepth.(LAYERS{dd}).pup.(WHNWH{ww}).mean_osc')
%         hold on; axis xy; box off
%         LogScale('y',10)
%         %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
%         clim([-0.2 0.2])
%         if dd == 6
%         xlabel('Pupil Size');
%         end
%         if ww == 1
%             ylabel('Freq')
%         end
%         if dd == 1
%         title((WHNWH{ww}))
%         end
%         
%        
%         
% end
% end
% 
% NiceSave('DepthSPECremandPup',figfolder,baseName)
% 
% %%
% 
% figure
% for dd = 1:6
% for oo = 1:2
% 
% subplot(6,4,(dd-1)*4+oo)
%         imagesc( SPECdepth.(LAYERS{dd}).(ONOFF{oo}).all.varbins,...
%             log10(SPECdepth.oscfreqs),...
%             SPECdepth.(LAYERS{dd}).(ONOFF{oo}).all.mean_osc')
%         hold on; axis xy; box off
%         plot([0 0],[0 max(SPECdepth.oscfreqs)],'w')
%         clim([-0.2 0.2])
%         LogScale('y',10)
%         if dd == 6
%         xlabel(['t - ',(ONOFF{oo})])
%         end
%         if oo == 1
%             ylabel({LAYERS{dd},'Freq'})
%         end
%         if dd == 1
%         title((ONOFF{oo}))
%         end
% end
%  
% subplot(6,4,(dd-1)*4+4)
%         imagesc( SPECdepth.(LAYERS{dd}).EMG.varbins,...
%             log10(SPECdepth.oscfreqs),...
%             SPECdepth.(LAYERS{dd}).EMG.mean_osc')
%         hold on; axis xy; box off
%         plot(SPECdepth.(LAYERS{dd}).EMG.varbins,SPECdepth.(LAYERS{dd}).EMG.vardist*1000,'w')
%         plot(log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],[0 max(SPECdepth.oscfreqs)],'k--')
%         clim([-0.2 0.2])
%         LogScale('y',10)
%         ylabel('Freq');
%                 if dd == 6
%         xlabel('EMG')
%         end
%         
% end
% NiceSave('DepthSPECremandWhisk',figfolder,baseName)
% 

%%  Mean depth spec by pupil size, phase and whisking
%prepare for LFPspec....
OSCdepth.freqs = spec.freqs;
for dd = 1:length(LAYERS)
for ww = 1:2
[ ~,OSCdepth.(LAYERS{dd}).pup.(WHNWH{ww}) ] = bz_LFPSpecToExternalVar( spec.LayerOsci(spec.(WHNWH{ww}),:,dd),...
    log10(spec.pup(spec.(WHNWH{ww}),:)),'specparms','input',...
    'figparms',[],'numvarbins',20,'varlim',[-0.25 0.25],'minX',500);

    for pp= 1:2
    [ ~,OSCdepth.(LAYERS{dd}).(HILO{pp}).(WHNWH{ww}) ] = bz_LFPSpecToExternalVar(...
        spec.LayerOsci(spec.(WHNWH{ww})&spec.(HILO{pp}),:,dd),...
        spec.pupphase(spec.(WHNWH{ww})&spec.(HILO{pp})),'specparms','input',...
        'figparms',[],'numvarbins',20,'varlim',[-pi pi],'minX',500);

    end
end

for oo = 1:2
    for ll = 1:2
        [ ~,OSCdepth.(LAYERS{dd}).(ONOFF{oo}).(LONGSHORT{ll}) ] = bz_LFPSpecToExternalVar(...
            spec.LayerOsci(spec.(LONGSHORT{ll}).(ONOFF{oo}),:,dd),...
            spec.whtime.(ONOFF{oo})(spec.(LONGSHORT{ll}).(ONOFF{oo}),:),'specparms','input',...
            'figparms',[],'numvarbins',40,'varlim',[-window window],'minX',500);
    end
    
    [ ~,OSCdepth.(LAYERS{dd}).(ONOFF{oo}).all ] = bz_LFPSpecToExternalVar(...
        spec.LayerOsci(:,:,dd),...
        spec.whtime.(ONOFF{oo}),'specparms','input',...
        'figparms',[],'numvarbins',400,'varlim',[-window window],'minX',500);
end

    
[ ~,OSCdepth.(LAYERS{dd}).EMG ] = bz_LFPSpecToExternalVar( spec.LayerOsci(:,:,dd),...
    log10(spec.EMG),'specparms','input',...
    'figparms',[],'numvarbins',40,'varlim',[-1.7 0.9],'minX',500);

[ ~,OSCdepth.(LAYERS{dd}).PSS ] = bz_LFPSpecToExternalVar( spec.LayerOsci(:,:,dd),...
    spec.LayerPSS(:,dd),'specparms','input',...
    'figparms',[],'numvarbins',40,'varlim',[-3 -0.5],'minX',500);
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
            log2(OSCdepth.freqs),...
            OSCdepth.(LAYERS{dd}).(HILO{pp}).(WHNWH{ww}).mean)
        hold on; axis xy; box off
        plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp),'w')
        end   
        LogScale('y',2)
        
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        clim([0 0.3])
        colorbar
        ylim([1 7])
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
            log2(OSCdepth.freqs),...
            OSCdepth.(LAYERS{dd}).pup.(WHNWH{ww}).mean)
        hold on; axis xy; box off
        LogScale('y',2)
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        clim([0 0.3])
        %colorbar
        ylim([1 7])
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
        clim([0 0.3])
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
        clim([0 0.3])
        LogScale('y',10)
        ylabel('Freq');
                if dd == 6
        xlabel('EMG')
        end
        
end
NiceSave('DepthOSCandWhisk',figfolder,baseName)


%%
%%
speclim = [0.8 1.2];
figure
for dd = 1:6
   subplot(6,3,(dd-1)*3+1)
        a = imagesc( SPECdepth.(LAYERS{dd}).PSS.varbins,...
            log10(SPECdepth.freqs),...
            SPECdepth.(LAYERS{dd}).PSS.mean);
        alpha(a,single(~isnan(SPECdepth.(LAYERS{dd}).PSS.mean)))
        hold on; axis xy; box off
        plot(SPECdepth.(LAYERS{dd}).PSS.varbins,SPECdepth.(LAYERS{dd}).PSS.vardist*10,'k','linewidth',2)
        %plot(log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],[0 max(SPECdepth.freqs)],'k--')
        crameri vik
        ColorbarWithAxis(speclim,'Power (med^-^1)')
        ylim([0 2.5])
        ylabel('Freq');
        title('Spec')
        LogScale('y',10)
                if dd == 6
        xlabel('PSS')
                end 
                ylabel({LAYERS{dd},'f (Hz)'})
        
                
   subplot(6,3,(dd-1)*3+2)
        a = imagesc( OSCdepth.(LAYERS{dd}).PSS.varbins,...
            log10(OSCdepth.freqs),...
            OSCdepth.(LAYERS{dd}).PSS.mean);
        alpha(a,single(~isnan(OSCdepth.(LAYERS{dd}).PSS.mean)))
        hold on; axis xy; box off
         plot(SPECdepth.(LAYERS{dd}).PSS.varbins,SPECdepth.(LAYERS{dd}).PSS.vardist*10,'k','linewidth',2)
        %plot(log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],[0 max(SPECdepth.freqs)],'k--')
        ColorbarWithAxis([-0.2 0.2],'PSS-subtract')
        title('OSC')
        crameri vik
        ylim([0 2.5])
        ylabel('Freq');
        LogScale('y',10)
                if dd == 6
        xlabel('PSS')
        end 
end
NiceSave('DepthSPECandPSS',figfolder,baseName)
