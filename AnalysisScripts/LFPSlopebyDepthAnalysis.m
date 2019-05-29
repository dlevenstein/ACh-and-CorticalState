function [PSScomponents,PSSdepth ] = LFPSlopebyDepthAnalysis(basePath,figfolder)
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
sessionInfo = bz_getSessionInfo(basePath);


%% Load the LFP if needed

% lfpchan = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID;
% downsamplefactor = 1;
% lfp = bz_GetLFP(lfpchan,...
%     'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);
%Noralize the LFP
%lfp.data = NormToInt(single(lfp.data),'modZ', SleepState.ints.NREMstate,lfp.samplingRate);

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

%%
inCTX = find(~isnan(depthinfo.ndepth));
CTXchans = depthinfo.channels(inCTX);
CTXdepth = -depthinfo.ndepth(inCTX);
clear PSS
[~,~,PSS.CTXchans] = intersect(CTXchans,PSpecSlope.Shortwin.channels,'stable');
PSS.data = PSpecSlope.Shortwin.PSS(:,PSS.CTXchans);
PSS.timestamps = PSpecSlope.Shortwin.timestamps;
PSS.samplingRate = 1/mean(diff(PSS.timestamps));
PSS.winsize = PSpecSlope.Shortwin.movingwin(1);
PSS.depth = CTXdepth;


%% Get rid of recording start/stop artifact and restrict to spontaneous behavior

%Restricting SPONT UP/DOWNs
load(fullfile(basePath,[baseName,'.MergePoints.events.mat']),'MergePoints');
sidx = find(startsWith(MergePoints.foldernames,"Spont"));
sponttimes = [MergePoints.timestamps(sidx(1),1) MergePoints.timestamps(sidx(end),2)];

inspont = InIntervals(PSS.timestamps,sponttimes);
PSS.timestamps = PSS.timestamps(inspont);
PSS.data = PSS.data(inspont,:);

maxtimejump = 1; %s
pupilcycle.amp = NanPadJumps( pupilcycle.timestamps,pupilcycle.amp,maxtimejump );
pupilcycle.phase = NanPadJumps( pupilcycle.timestamps,pupilcycle.phase,maxtimejump );
pupildilation.data = NanPadJumps( pupildilation.timestamps,pupildilation.data,maxtimejump );

PSS.Wh = InIntervals(PSS.timestamps,EMGwhisk.ints.Wh);
EMGwhisk.ints.ExpWh = bsxfun(@plus,EMGwhisk.ints.Wh,[-0.5 0.5]*PSS.winsize);
PSS.NWh = ~InIntervals(PSS.timestamps,EMGwhisk.ints.ExpWh);
PSS.hipup = interp1(pupilcycle.timestamps,single(pupilcycle.highpup),PSS.timestamps,'nearest')==1;
PSS.lopup = interp1(pupilcycle.timestamps,single(~pupilcycle.highpup),PSS.timestamps,'nearest')==1;

%%
PSS.pupphase = interp1(pupilcycle.timestamps,pupilcycle.phase,PSS.timestamps,'nearest');
PSS.pup = interp1(pupildilation.timestamps,pupildilation.data,PSS.timestamps,'nearest');


%% Get relative time to nearest Wh onset
window = 5; %s
PSS.whtime = nan(size(PSS.timestamps));
for wh = 2:size(EMGwhisk.ints.Wh,1)
    wh
    inwintimes = PSS.timestamps > EMGwhisk.ints.Wh(wh-1,2) & ...
        PSS.timestamps <= EMGwhisk.ints.Wh(wh,2);
    
    PSS.whtime(inwintimes) = PSS.timestamps(inwintimes)-EMGwhisk.ints.Wh(wh,1);
end


%% Interpolate PSS to normalized depth
PSS.interpdepth = linspace(-1,0,100);
PSS.depthinterp = interp1(PSS.depth',PSS.data',PSS.interpdepth')';

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
WHNWH = {'Wh','NWh'};
HILO = {'lopup','hipup'};
%% Mean depth activation by pupil size, phase and whisking
%prepare for LFPspec....
PSSdepth.depth = PSS.interpdepth;
for ww = 1:2
[ ~,PSSdepth.pup.(WHNWH{ww}) ] = bz_LFPSpecToExternalVar( PSS.data(PSS.(WHNWH{ww}),:),...
    log10(PSS.pup(PSS.(WHNWH{ww}),:)),'specparms','input',...
    'figparms',true,'numvarbins',20,'varlim',[-0.25 0.25]);
    PSSdepth.pup.(WHNWH{ww}).mean_interp = interp1(PSS.depth',PSSdepth.pup.(WHNWH{ww}).mean,PSS.interpdepth');

    for pp= 1:2
    [ ~,PSSdepth.(HILO{pp}).(WHNWH{ww}) ] = bz_LFPSpecToExternalVar(...
        PSS.data(PSS.(WHNWH{ww})&PSS.(HILO{pp}),:),...
        PSS.pupphase(PSS.(WHNWH{ww})&PSS.(HILO{pp}),:),'specparms','input',...
        'figparms',true,'numvarbins',20,'varlim',[-pi pi]);
    PSSdepth.(HILO{pp}).(WHNWH{ww}).mean_interp = interp1(PSS.depth',PSSdepth.(HILO{pp}).(WHNWH{ww}).mean,PSS.interpdepth');

    end
end

[ ~,PSSdepth.whOn ] = bz_LFPSpecToExternalVar( PSS.data,...
    PSS.whtime,'specparms','input',...
    'figparms',true,'numvarbins',50,'varlim',[-window window]);
    PSSdepth.whOn.mean_interp = interp1(PSS.depth',PSSdepth.whOn.mean,PSS.interpdepth');
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
    ColorbarWithAxis([-2.75 -0.5],'PSS')
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


figure
for ww = 1:2
    subplot(3,2,(ww-1)*2+1)
        for pp = 1:2
        imagesc( PSSdepth.(HILO{pp}).(WHNWH{ww}).varbins+2*pi*(pp-1),...
            PSS.interpdepth,...
            PSSdepth.(HILO{pp}).(WHNWH{ww}).mean_interp)
        hold on; axis xy; box off
        plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp)-1,'k')
        end   
        ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        xlim([-pi 3*pi])
        xlabel('Pupil Phase');ylabel({WHNWH{ww},'Depth'})

        
        
    subplot(3,2,(ww-1)*2+2)
        imagesc( PSSdepth.pup.(WHNWH{ww}).varbins,...
            PSS.interpdepth,...
            PSSdepth.pup.(WHNWH{ww}).mean_interp)
        hold on; axis xy; box off
        ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        xlabel('Pupil Size');ylabel('Depth')
end

subplot(3,2,5.5)
        imagesc( PSSdepth.whOn.varbins,...
            PSS.interpdepth,...
            PSSdepth.whOn.mean_interp)
        hold on; axis xy; box off
        plot([0 0],[-1 0],'w')
        ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        xlabel('t (s - relative to WhOn');ylabel('Depth')
NiceSave('DepthPSSandBeh',figfolder,baseName)

