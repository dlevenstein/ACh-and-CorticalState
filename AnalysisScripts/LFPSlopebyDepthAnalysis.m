function [PSScomponents,PSSdepth,PSSdist ] = LFPSlopebyDepthAnalysis(basePath,figfolder)
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
lowfilter = [0.01 0.1]; %old. order 3
lowfilter = [0.02 0.2]; %new: EMG coupled. order 1
%highfilter = [0.3 0.8];

pupil4filter = pupildilation;
pupilcycle = bz_Filter(pupil4filter,'passband',lowfilter,'filter' ,'fir1','order',1);
%highpupildata = bz_Filter(pupil4filter,'passband',highfilter,'filter' ,'fir1');
pupilcycle.pupthresh = -0.8;
pupilcycle.highpup = log10(pupilcycle.amp)>pupilcycle.pupthresh; 
pupilcycle.lowpup = log10(pupilcycle.amp)<=pupilcycle.pupthresh; 

% EMG
EMGwhisk = bz_LoadBehavior(basePath,'EMGwhisk');

%% Load PSS and aligning behavior
%load(fullfile(basePath,[baseName,'.PowerSpectrumSlope.mat']));
load(fullfile(basePath,[baseName,'.PowerSpectrumSlope.lfp.mat']));
depthinfo = rescaleCx(basePath);

%%
inCTX = find(~isnan(depthinfo.ndepth));
CTXchans = depthinfo.channels(inCTX);
CTXdepth = -depthinfo.ndepth(inCTX);
clear PSS
%[~,~,PSS.CTXchans] = intersect(CTXchans,PSpecSlope.Shortwin.channels,'stable');
[~,~,PSS.CTXchans] = intersect(CTXchans,specslope.channels,'stable');
PSS.data = specslope.data(:,PSS.CTXchans);
PSS.timestamps = specslope.timestamps;
PSS.samplingRate = 1/mean(diff(PSS.timestamps));
PSS.winsize = specslope.detectionparms.winsize;
PSS.depth = CTXdepth;

PSS.osci = specslope.resid(:,:,PSS.CTXchans);
%PSS.osci = shiftdim(PSS.osci,1);
PSS.freqs = specslope.freqs;
PSS.chanlayers = depthinfo.layer(inCTX);

LAYERS = depthinfo.lnames;

for ll = 1:length(LAYERS)
    PSS.Lchans.(LAYERS{ll}) = strcmp(LAYERS{ll},PSS.chanlayers);
    PSS.oscLayer(:,:,ll) = squeeze(mean(PSS.osci(:,PSS.Lchans.(LAYERS{ll}),:),2));
end

% Interpolate PSS to normalized depth
PSS.interpdepth = linspace(-1,0,100);
PSS.depthinterp = interp1(PSS.depth',PSS.data',PSS.interpdepth')';

specslope = []; %clear memory
%% Get rid of recording start/stop artifact and restrict to spontaneous behavior

%Restricting SPONT UP/DOWNs
load(fullfile(basePath,[baseName,'.MergePoints.events.mat']),'MergePoints');
sidx = find(startsWith(MergePoints.foldernames,"Spont"));
sponttimes = [MergePoints.timestamps(sidx(1),1) MergePoints.timestamps(sidx(end),2)];

inspont = InIntervals(PSS.timestamps,sponttimes);
PSS.timestamps = PSS.timestamps(inspont);
PSS.data = PSS.data(inspont,:);
PSS.oscLayer = PSS.oscLayer(inspont,:,:);
PSS.osci = PSS.osci(inspont,:,:);

maxtimejump = 1; %s
pupilcycle.amp = NanPadJumps( pupilcycle.timestamps,pupilcycle.amp,maxtimejump );
pupilcycle.phase = NanPadJumps( pupilcycle.timestamps,pupilcycle.phase,maxtimejump );
pupildilation.data = NanPadJumps( pupildilation.timestamps,pupildilation.data,maxtimejump );
EMGwhisk.EMGsm = NanPadJumps( EMGwhisk.timestamps,EMGwhisk.EMGsm,maxtimejump );

PSS.Wh = InIntervals(PSS.timestamps,EMGwhisk.ints.Wh);
EMGwhisk.ints.ExpWh = bsxfun(@plus,EMGwhisk.ints.Wh,[-0.5 0.5]*PSS.winsize);
PSS.NWh = ~InIntervals(PSS.timestamps,EMGwhisk.ints.ExpWh);
PSS.hipup = interp1(pupilcycle.timestamps,single(pupilcycle.highpup),PSS.timestamps,'nearest')==1;
PSS.lopup = interp1(pupilcycle.timestamps,single(pupilcycle.lowpup),PSS.timestamps,'nearest')==1;


PSS.pupphase = interp1(pupilcycle.timestamps,pupilcycle.phase,PSS.timestamps,'nearest');
PSS.pup = interp1(pupildilation.timestamps,pupildilation.data,PSS.timestamps,'nearest');
PSS.EMG = interp1(EMGwhisk.timestamps,EMGwhisk.EMGsm,PSS.timestamps,'nearest');

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

%% Get Whisk On/Offset aligned PSS, separated by long/short whisks

window = 5; %s
durthresh = 1; %
PSS.whtime.WhOn = nan(size(PSS.timestamps));
PSS.whtime.WhOFF = nan(size(PSS.timestamps));
PSS.long.WhOn = false(size(PSS.timestamps));
PSS.long.WhOFF = false(size(PSS.timestamps));
PSS.short.WhOn = false(size(PSS.timestamps));
PSS.short.WhOFF = false(size(PSS.timestamps));
for wh = 2:(size(EMGwhisk.ints.Wh,1)-1)
    %wh
    longwhisk = EMGwhisk.dur(wh)>durthresh;
    
    for oo = 1:2
        if oo == 1
            inwintimes = PSS.timestamps > EMGwhisk.ints.ExpWh(wh-1,2) & ...
                PSS.timestamps <= EMGwhisk.ints.Wh(wh,2);  
        elseif oo==2
            inwintimes = PSS.timestamps > EMGwhisk.ints.Wh(wh,1) & ...
                PSS.timestamps <= EMGwhisk.ints.ExpWh(wh+1,1); 
        end
        PSS.whtime.(ONOFF{oo})(inwintimes) = PSS.timestamps(inwintimes)-EMGwhisk.ints.Wh(wh,oo);
        PSS.long.(ONOFF{oo})(inwintimes) = longwhisk;
        PSS.short.(ONOFF{oo})(inwintimes) = ~longwhisk;
    end
    
end



%% Identifying components: 
%PSS correlation by depth
PSScomponents.corr = corr(PSS.depthinterp);
PSScomponents.depth = PSS.interpdepth;

%PCA
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(PSS.data);

%ICA
% [weights,sphere,compvars,bias,signs,lrates,activations] ...
%                               = runica(PSS.data','pca',10,'sphering','off');
%                           
nPC = 10;
PSScomponents.PCAcoeff = interp1(PSS.depth',COEFF(:,1:nPC),PSS.interpdepth');
PSScomponents.EV = EXPLAINED(1:nPC);

%%
% figure
% subplot(2,2,1)
% plot(compvars,'o-')
% xlim([0 6])
% xlabel('IC');ylabel('% EV')
% subplot(2,2,2)
% plot(activations(1,:),activations(2,:),'.')
% xlabel('IC 1');ylabel('IC 2')
% 
% subplot(2,2,3)
% plot(PSS.depth,weights(1:3,:))
% hold on
% plot(PSS.depth([1 end]),[0 0],'k--')                

%%
% figure

% for pp = 1:2
%     plot(EMGwhisk.pupphase(EMGwhisk.(HILO{pp}))+2*pi*(pp-1),log10(EMGwhisk.dur(EMGwhisk.(HILO{pp}))),'k.')
%     hold on; axis xy; box off
%     plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp)-1,'k')
% end 
% LogScale('y',10)
% xlabel('Pupil Phase');ylabel('Duration')
% xlim([-pi 3*pi])



%% Mean depth activation by pupil size, phase and whisking
%prepare for LFPspec....
PSSdepth.depth = PSS.interpdepth;
for ww = 1:2
[ ~,PSSdepth.pup.(WHNWH{ww}) ] = bz_LFPSpecToExternalVar( PSS.data(PSS.(WHNWH{ww}),:),...
    log10(PSS.pup(PSS.(WHNWH{ww}),:)),'specparms','input',...
    'figparms',[],'numvarbins',20,'varlim',[-0.25 0.25]);
    PSSdepth.pup.(WHNWH{ww}).mean_interp = interp1(PSS.depth',PSSdepth.pup.(WHNWH{ww}).mean,PSS.interpdepth');
    PSSdepth.pup.(WHNWH{ww}).std_interp = interp1(PSS.depth',PSSdepth.pup.(WHNWH{ww}).std,PSS.interpdepth');
    
    for pp= 1:2
    [ ~,PSSdepth.(HILO{pp}).(WHNWH{ww}) ] = bz_LFPSpecToExternalVar(...
        PSS.data(PSS.(WHNWH{ww})&PSS.(HILO{pp}),:),...
        PSS.pupphase(PSS.(WHNWH{ww})&PSS.(HILO{pp}),:),'specparms','input',...
        'figparms',[],'numvarbins',20,'varlim',[-pi pi]);
    PSSdepth.(HILO{pp}).(WHNWH{ww}).mean_interp = interp1(PSS.depth',PSSdepth.(HILO{pp}).(WHNWH{ww}).mean,PSS.interpdepth');
    PSSdepth.(HILO{pp}).(WHNWH{ww}).std_interp = interp1(PSS.depth',PSSdepth.(HILO{pp}).(WHNWH{ww}).std,PSS.interpdepth');

    end
end

for oo = 1:2
    for ll = 1:2
        [ ~,PSSdepth.(ONOFF{oo}).(LONGSHORT{ll}) ] = bz_LFPSpecToExternalVar(...
            PSS.data(PSS.(LONGSHORT{ll}).(ONOFF{oo}),:),...
            PSS.whtime.(ONOFF{oo})(PSS.(LONGSHORT{ll}).(ONOFF{oo}),:),'specparms','input',...
            'figparms',[],'numvarbins',40,'varlim',[-window window]);
            PSSdepth.(ONOFF{oo}).(LONGSHORT{ll}).mean_interp = interp1(PSS.depth',PSSdepth.(ONOFF{oo}).(LONGSHORT{ll}).mean,PSS.interpdepth');
            PSSdepth.(ONOFF{oo}).(LONGSHORT{ll}).std_interp = interp1(PSS.depth',PSSdepth.(ONOFF{oo}).(LONGSHORT{ll}).std,PSS.interpdepth');
    end
    
    [ ~,PSSdepth.(ONOFF{oo}).all ] = bz_LFPSpecToExternalVar(...
        PSS.data,...
        PSS.whtime.(ONOFF{oo}),'specparms','input',...
        'figparms',[],'numvarbins',40,'varlim',[-window window]);
        PSSdepth.(ONOFF{oo}).all.mean_interp = interp1(PSS.depth',PSSdepth.(ONOFF{oo}).all.mean,PSS.interpdepth');
        PSSdepth.(ONOFF{oo}).all.std_interp = interp1(PSS.depth',PSSdepth.(ONOFF{oo}).all.std,PSS.interpdepth');
end

    
[ ~,PSSdepth.EMG ] = bz_LFPSpecToExternalVar( PSS.data,...
    log10(PSS.EMG),'specparms','input',...
    'figparms',[],'numvarbins',40,'varlim',[-1.7 0.9]);
    PSSdepth.EMG.mean_interp = interp1(PSS.depth',PSSdepth.EMG.mean,PSS.interpdepth');
    PSSdepth.EMG.std_interp = interp1(PSS.depth',PSSdepth.EMG.std,PSS.interpdepth');

    

    
%%
xwinsize = 80;
xwin = bz_RandomWindowInIntervals(PSS.timestamps([1 end]),xwinsize);
figure
subplot(2,1,1)
    imagesc(PSS.timestamps,PSS.interpdepth,PSS.depthinterp')
    hold on
    plot(pupildilation.timestamps,bz_NormToRange(pupildilation.data,[-1 -0.5]),'w','linewidth',2)
    plot(EMGwhisk.timestamps,bz_NormToRange(EMGwhisk.EMG,[-1.3 -0.8]),'k')
    xlim(xwin)
    axis xy
    ylim([-1.1 0])
    ColorbarWithAxis([-2 -0.25],'PSS')
    xlabel('t (s)');ylabel('Depth')
    
    subplot(2,2,3)
    imagesc(PSScomponents.depth,PSScomponents.depth,PSScomponents.corr)
    colorbar
    axis xy
    xlabel('Depth');ylabel('Depth')
    title('PSS Corr')
    
%     subplot(2,2,4)
%         imagesc(PSSdepth.byPup.varbins,PSS.interpdepth,PSSdepth.byPup.mean)
%         axis xy
subplot(4,2,6)
    plot(log10(PSScomponents.EV),'o-')
    xlim([1 6])
    xlabel('PC');ylabel('% EV')
    LogScale('y',10)

subplot(4,2,8)
    plot(PSScomponents.depth,PSScomponents.PCAcoeff(:,1:3))
    hold on
    plot(PSScomponents.depth([1 end]),[0 0],'k--')
    legend('PC1','PC2','PC3','location','eastoutside')
    xlabel('Deptah');ylabel('Weight')
    
        NiceSave('DepthPSSEx',figfolder,baseName)


%%
cosx = linspace(-pi,pi,100);
cospamp = [0.025 0.3];

meanPSSrange = [-1 -0.5];

figure
for ww = 1:2
    subplot(4,2,(ww-1)*2+1)
        for pp = 1:2
        imagesc( PSSdepth.(HILO{pp}).(WHNWH{ww}).varbins+2*pi*(pp-1),...
            PSS.interpdepth,...
            PSSdepth.(HILO{pp}).(WHNWH{ww}).mean_interp)
        hold on; axis xy; box off
        plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp)-1,'k')
        end   
        ColorbarWithAxis(meanPSSrange,'Mean PSS')
        xlim([-pi 3*pi])
        xlabel('Pupil Phase');ylabel({WHNWH{ww},'Depth'})

        
        
    subplot(4,2,(ww-1)*2+2)
        imagesc( PSSdepth.pup.(WHNWH{ww}).varbins,...
            PSS.interpdepth,...
            PSSdepth.pup.(WHNWH{ww}).mean_interp)
        hold on; axis xy; box off
        ColorbarWithAxis(meanPSSrange,'Mean PSS')
        xlabel('Pupil Size');ylabel('Depth')
        
        
    subplot(4,2,(ww-1)*2+1+4)
        for pp = 1:2
        imagesc( PSSdepth.(HILO{pp}).(WHNWH{ww}).varbins+2*pi*(pp-1),...
            PSS.interpdepth,...
            PSSdepth.(HILO{pp}).(WHNWH{ww}).std_interp)
        hold on; axis xy; box off
        plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp)-1,'k')
        end   
        ColorbarWithAxis([0 0.5],'STD PSS')

        xlim([-pi 3*pi])
        xlabel('Pupil Phase');ylabel({WHNWH{ww},'Depth'})

        
        
    subplot(4,2,(ww-1)*2+2+4)
        imagesc( PSSdepth.pup.(WHNWH{ww}).varbins,...
            PSS.interpdepth,...
            PSSdepth.pup.(WHNWH{ww}).std_interp)
        hold on; axis xy; box off
        ColorbarWithAxis([0 0.5],'STD PSS')
        xlabel('Pupil Size');ylabel('Depth')
        
end


NiceSave('DepthPSSandPup',figfolder,baseName)


%%

figure

for oo = 1:2
    for ll= 1:2
subplot(4,3,(ll-1)*3+oo)
        imagesc( PSSdepth.(ONOFF{oo}).(LONGSHORT{ll}).varbins,...
            PSSdepth.depth,...
            PSSdepth.(ONOFF{oo}).(LONGSHORT{ll}).mean_interp)
        hold on; axis xy; box off
        plot([0 0],[-1 0],'w')
        ColorbarWithAxis(meanPSSrange,'Mean PSS')
        xlabel(['t - ',(ONOFF{oo})]);ylabel('Depth')
    end
end
 
subplot(4,3,6)
        imagesc( PSSdepth.EMG.varbins,...
            PSSdepth.depth,...
            PSSdepth.EMG.mean_interp)
        hold on; axis xy; box off
        plot(PSSdepth.EMG.varbins,PSSdepth.EMG.vardist*10-1,'k')
        plot(log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],[-1 0],'k--')
        ColorbarWithAxis(meanPSSrange,'Mean PSS')
        xlabel('EMG');ylabel('Depth')
        
for oo = 1:2
    for ll= 1:2
subplot(4,3,(ll-1)*3+oo+6)
        imagesc( PSSdepth.(ONOFF{oo}).(LONGSHORT{ll}).varbins,...
            PSSdepth.depth,...
            PSSdepth.(ONOFF{oo}).(LONGSHORT{ll}).std_interp)
        hold on; axis xy; box off
        plot([0 0],[-1 0],'w')
        ColorbarWithAxis([0 0.5],'STD PSS')
        xlabel(['t - ',(ONOFF{oo})]);ylabel('Depth')
    end
end
 
subplot(4,3,9)
        imagesc( PSSdepth.EMG.varbins,...
            PSSdepth.depth,...
            PSSdepth.EMG.std_interp)
        hold on; axis xy; box off
        plot(PSSdepth.EMG.varbins,PSSdepth.EMG.vardist*10-1,'k')
        plot(log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],[-1 0],'k--')
        ColorbarWithAxis([0 0.5],'STD PSS')
        xlabel('EMG');ylabel('Depth')
        
        
NiceSave('DepthPSSandWhisk',figfolder,baseName)

%% PSS distirbution conditioned on pupil

%Get the best l5 channel

PSSdist.depth = PSS.interpdepth;
PSSrange = [-1.5 0];

%for cc = 1:length(PSSdist.depth)

%for cc = 1:length(CTXchans)
%bz_Counter(cc,length(CTXchans),'Conditional Dist, Channel')
for ll = 1:length(LAYERS)
for ww = 1:2
PSSdist.(LAYERS{ll}).pup.(WHNWH{ww})  = ConditionalHist( log10(PSS.pup(PSS.(WHNWH{ww}))),...
    PSS.data(PSS.(WHNWH{ww}),PSS.Lchans.(LAYERS{ll})),...
    'numXbins',20,'Xbounds',[-0.25 0.25],'Ybounds',PSSrange,'numYbins',50);
    PSSdist.(LAYERS{ll}).pup.(WHNWH{ww}) =  bz_CollapseStruct(PSSdist.(LAYERS{ll}).pup.(WHNWH{ww}),...
        3,'mean',true);

    for pp= 1:2
    PSSdist.(LAYERS{ll}).(HILO{pp}).(WHNWH{ww})  = ConditionalHist(...
        PSS.pupphase(PSS.(WHNWH{ww})&PSS.(HILO{pp}) ),...
        PSS.data(PSS.(WHNWH{ww})&PSS.(HILO{pp}),PSS.Lchans.(LAYERS{ll})),...
        'Ybounds',PSSrange,'numYbins',50,'numXbins',20,'Xbounds',[-pi pi]);
    PSSdist.(LAYERS{ll}).(HILO{pp}).(WHNWH{ww}) =  bz_CollapseStruct(PSSdist.(LAYERS{ll}).(HILO{pp}).(WHNWH{ww}),...
        3,'mean',true);
    end
end

    
PSSdist.(LAYERS{ll}).EMG  = ConditionalHist( log10(PSS.EMG),PSS.data(:,PSS.Lchans.(LAYERS{ll})),...
    'Ybounds',PSSrange,'numYbins',50,'numXbins',40,'Xbounds',[-1.7 0.9]);
    PSSdist.(LAYERS{ll}).EMG =  bz_CollapseStruct(PSSdist.(LAYERS{ll}).EMG,...
        3,'mean',true);
    
end

    

%%
cosx = linspace(-pi,pi,100);
cospamp = [0.05 0.5];


figure
for ll = 1:length(LAYERS)
for ww = 1:2
    subplot(6,5,(ll-1)*5+ww)
        for pp = 1:2
        imagesc( PSSdist.(LAYERS{ll}).(HILO{pp}).(WHNWH{ww}).Xbins+2*pi*(pp-1),...
            PSSdist.(LAYERS{ll}).(HILO{pp}).(WHNWH{ww}).Ybins,...
            PSSdist.(LAYERS{ll}).(HILO{pp}).(WHNWH{ww}).pYX')
        hold on; axis xy; box off
        plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp)-3.1,'w')
        end   
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        xlim([-pi 3*pi])
        xlabel('Pupil Phase');ylabel({WHNWH{ww},'PSS'})
    crameri bilbao
        
        
    subplot(6,5,(ll-1)*5+ww+2)
        imagesc( PSSdist.(LAYERS{ll}).pup.(WHNWH{ww}).Xbins,...
            PSSdist.(LAYERS{ll}).(HILO{pp}).(WHNWH{ww}).Ybins,...
            PSSdist.(LAYERS{ll}).pup.(WHNWH{ww}).pYX')
        hold on; axis xy; box off
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        xlabel('Pupil Size');ylabel('PSS')
   crameri bilbao
        
        
    subplot(6,5,(ll-1)*5+5)
        imagesc( PSSdist.(LAYERS{ll}).EMG.Xbins,...
            PSSdist.(LAYERS{ll}).EMG.Ybins,...
            PSSdist.(LAYERS{ll}).EMG.pYX')
        hold on; axis xy; box off
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        xlabel('EMG');ylabel('PSS')
        
    crameri bilbao
end
end

NiceSave('LayerPSSandBehavior',figfolder,baseName)
%%
% 
% %% Mean depth osci by pupil size, phase and whisking
% %prepare for LFPspec....
% OSCdepth.freqs = PSS.freqs;
% for dd = 1:length(LAYERS)
% for ww = 1:2
% [ ~,OSCdepth.(LAYERS{dd}).pup.(WHNWH{ww}) ] = bz_LFPSpecToExternalVar( PSS.oscLayer(PSS.(WHNWH{ww}),:,dd),...
%     log10(PSS.pup(PSS.(WHNWH{ww}),:)),'specparms','input',...
%     'figparms',[],'numvarbins',20,'varlim',[-0.25 0.25]);
% 
%     for pp= 1:2
%     [ ~,OSCdepth.(LAYERS{dd}).(HILO{pp}).(WHNWH{ww}) ] = bz_LFPSpecToExternalVar(...
%         PSS.oscLayer(PSS.(WHNWH{ww})&PSS.(HILO{pp}),:,dd),...
%         PSS.pupphase(PSS.(WHNWH{ww})&PSS.(HILO{pp})),'specparms','input',...
%         'figparms',[],'numvarbins',20,'varlim',[-pi pi]);
% 
%     end
% end
% 
% for oo = 1:2
%     for ll = 1:2
%         [ ~,OSCdepth.(LAYERS{dd}).(ONOFF{oo}).(LONGSHORT{ll}) ] = bz_LFPSpecToExternalVar(...
%             PSS.oscLayer(PSS.(LONGSHORT{ll}).(ONOFF{oo}),:,dd),...
%             PSS.whtime.(ONOFF{oo})(PSS.(LONGSHORT{ll}).(ONOFF{oo}),:),'specparms','input',...
%             'figparms',[],'numvarbins',40,'varlim',[-window window]);
%     end
%     
%     [ ~,OSCdepth.(LAYERS{dd}).(ONOFF{oo}).all ] = bz_LFPSpecToExternalVar(...
%         PSS.oscLayer(:,:,dd),...
%         PSS.whtime.(ONOFF{oo}),'specparms','input',...
%         'figparms',[],'numvarbins',40,'varlim',[-window window]);
% end
% 
%     
% [ ~,OSCdepth.(LAYERS{dd}).EMG ] = bz_LFPSpecToExternalVar( PSS.oscLayer(:,:,dd),...
%     log10(PSS.EMG),'specparms','input',...
%     'figparms',[],'numvarbins',40,'varlim',[-1.7 0.9]);
% end
% 
%   %%
%   xwinsize = 30;
% xwin = bz_RandomWindowInIntervals(PSS.timestamps([1 end]),xwinsize);
% figure
%   for cc = 1:6
% 
% subplot(6,1,cc)
%     imagesc(PSS.timestamps,PSS.freqs,(PSS.oscLayer(:,:,cc))')
%     hold on
%     %plot(PSS.timestamps,PSS.depthinterp(:,cc)')
% 
%     plot(pupildilation.timestamps,bz_NormToRange(pupildilation.data,[-1 -0.5]),'w','linewidth',2)
%     plot(EMGwhisk.timestamps,bz_NormToRange(EMGwhisk.EMG,[-1.3 -0.8]),'k')
%     xlim(xwin)
%     axis xy
%     %ylim([-1.1 0])
%     %ColorbarWithAxis([-500 800000],'Osci')
%     xlabel('t (s)');ylabel('Depth')
%     
%   end
%   
%   
%   %%
%   
%   cosx = linspace(-pi,pi,100);
% cospamp = [0.025 0.3].*100;
% 
% 
% figure
% for dd = 1:6
% for ww = 1:2
%     subplot(6,4,(dd-1)*4+ww)
%         for pp = 1:2
%         imagesc( OSCdepth.(LAYERS{dd}).(HILO{pp}).(WHNWH{ww}).varbins+2*pi*(pp-1),...
%             PSS.freqs,...
%             OSCdepth.(LAYERS{dd}).(HILO{pp}).(WHNWH{ww}).mean)
%         hold on; axis xy; box off
%         plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp),'w')
%         end   
%         %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
%         clim([0 5e4])
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
%         imagesc( OSCdepth.(LAYERS{dd}).pup.(WHNWH{ww}).varbins,...
%             PSS.freqs,...
%             OSCdepth.(LAYERS{dd}).pup.(WHNWH{ww}).mean)
%         hold on; axis xy; box off
%         %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
%         clim([0 5e4])
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
% NiceSave('DepthOSCandPup',figfolder,baseName)
% 
% %%
% 
% figure
% for dd = 1:6
% for oo = 1:2
% 
% subplot(6,4,(dd-1)*4+oo)
%         imagesc( OSCdepth.(LAYERS{dd}).(ONOFF{oo}).all.varbins,...
%             OSCdepth.freqs,...
%             OSCdepth.(LAYERS{dd}).(ONOFF{oo}).all.mean)
%         hold on; axis xy; box off
%         plot([0 0],[0 max(PSS.freqs)],'w')
%         clim([0 5e4])
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
%         imagesc( OSCdepth.(LAYERS{dd}).EMG.varbins,...
%             OSCdepth.freqs,...
%             OSCdepth.(LAYERS{dd}).EMG.mean)
%         hold on; axis xy; box off
%         plot(OSCdepth.(LAYERS{dd}).EMG.varbins,OSCdepth.(LAYERS{dd}).EMG.vardist*1000,'w')
%         plot(log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],[0 max(PSS.freqs)],'k--')
%         clim([0 5e4])
%         ylabel('Freq');
%                 if dd == 6
%         xlabel('EMG')
%         end
%         
% end
% NiceSave('DepthOSCandWhisk',figfolder,baseName)
