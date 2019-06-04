function [ ] = FBandsAnalysis(basePath,figfolder)
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
lfp.depth = CTXdepth;




%%

BANDS = {'HiGamma','LoGamma','Theta','Delta','Slow'};
bandranges = [[100 312];[30 55];[6 12];[2 5];[0.5 2]];

banddownsamples = [1 1 5 5 5];



for bb = 1:length(BANDS)
    bb
    
    lfp.(BANDS{bb}) = bz_Filter(bz_DownsampleLFP( lfp, banddownsamples(bb) ),...
        'passband',bandranges(bb,:),'filter','fir1','fast',false);
    lfp.(BANDS{bb}).amp = log10(lfp.(BANDS{bb}).amp);
end


%% Pad for jumps
maxtimejump = 1; %s
pupilcycle.amp = NanPadJumps( pupilcycle.timestamps,pupilcycle.amp,maxtimejump );
pupilcycle.phase = NanPadJumps( pupilcycle.timestamps,pupilcycle.phase,maxtimejump );
pupildilation.data = NanPadJumps( pupildilation.timestamps,pupildilation.data,maxtimejump );
EMGwhisk.EMGsm = NanPadJumps( EMGwhisk.timestamps,EMGwhisk.EMGsm,maxtimejump );

%Interpolate all behavioral variables at band timepoints
for bb = length(BANDS):-1:1
lfp.(BANDS{bb}).Wh = InIntervals(lfp.(BANDS{bb}).timestamps,EMGwhisk.ints.Wh);
EMGwhisk.ints.ExpWh = bsxfun(@plus,EMGwhisk.ints.Wh,[-0.5 0.5]*1);
lfp.(BANDS{bb}).NWh = ~InIntervals(lfp.(BANDS{bb}).timestamps,EMGwhisk.ints.ExpWh);
lfp.(BANDS{bb}).hipup = interp1(pupilcycle.timestamps,single(pupilcycle.highpup),lfp.(BANDS{bb}).timestamps,'nearest')==1;
lfp.(BANDS{bb}).lopup = interp1(pupilcycle.timestamps,single(~pupilcycle.highpup),lfp.(BANDS{bb}).timestamps,'nearest')==1;


lfp.(BANDS{bb}).pupphase = interp1(pupilcycle.timestamps,pupilcycle.phase,lfp.(BANDS{bb}).timestamps,'nearest');
lfp.(BANDS{bb}).pup = interp1(pupildilation.timestamps,pupildilation.data,lfp.(BANDS{bb}).timestamps,'nearest');
lfp.(BANDS{bb}).EMG = interp1(EMGwhisk.timestamps,EMGwhisk.EMGsm,lfp.(BANDS{bb}).timestamps,'nearest');

end




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

%% Get Whisk On/Offset aligned Band data, separated by long/short whisks

window = 2; %s
durthresh = 1; %
for bb = 1:length(BANDS)
lfp.(BANDS{bb}).whtime.WhOn = nan(size(lfp.(BANDS{bb}).timestamps));
lfp.(BANDS{bb}).whtime.WhOFF = nan(size(lfp.(BANDS{bb}).timestamps));
lfp.(BANDS{bb}).long.WhOn = false(size(lfp.(BANDS{bb}).timestamps));
lfp.(BANDS{bb}).long.WhOFF = false(size(lfp.(BANDS{bb}).timestamps));
lfp.(BANDS{bb}).short.WhOn = false(size(lfp.(BANDS{bb}).timestamps));
lfp.(BANDS{bb}).short.WhOFF = false(size(lfp.(BANDS{bb}).timestamps));
for wh = 2:(size(EMGwhisk.ints.Wh,1)-1)
    %wh
    longwhisk = EMGwhisk.dur(wh)>durthresh;
    
    for oo = 1:2
        if oo == 1
            inwintimes = lfp.(BANDS{bb}).timestamps > EMGwhisk.ints.ExpWh(wh-1,2) & ...
                lfp.(BANDS{bb}).timestamps <= EMGwhisk.ints.Wh(wh,2);  
        elseif oo==2
            inwintimes = lfp.(BANDS{bb}).timestamps > EMGwhisk.ints.Wh(wh,1) & ...
                lfp.(BANDS{bb}).timestamps <= EMGwhisk.ints.ExpWh(wh+1,1); 
        end
        lfp.(BANDS{bb}).whtime.(ONOFF{oo})(inwintimes) = lfp.(BANDS{bb}).timestamps(inwintimes)-EMGwhisk.ints.Wh(wh,oo);
        lfp.(BANDS{bb}).long.(ONOFF{oo})(inwintimes) = longwhisk;
        lfp.(BANDS{bb}).short.(ONOFF{oo})(inwintimes) = ~longwhisk;
    end
    
end

end


%%
lfp.interpdepth = linspace(-1,0,100);


%% Mean depth band by pupil size, phase and whisking
%prepare for LFPspec....
for bb = 1:length(BANDS)
BANDdepth.(BANDS{bb}).depth = lfp.interpdepth;
for ww = 1:2
[ ~,BANDdepth.(BANDS{bb}).pup.(WHNWH{ww}) ] = bz_LFPSpecToExternalVar( lfp.(BANDS{bb}).amp(lfp.(BANDS{bb}).(WHNWH{ww}),:),...
    log10(lfp.(BANDS{bb}).pup(lfp.(BANDS{bb}).(WHNWH{ww}),:)),'specparms','input',...
    'figparms',[],'numvarbins',25,'varlim',[-0.25 0.25]);
    BANDdepth.(BANDS{bb}).pup.(WHNWH{ww}).mean_interp = interp1(lfp.depth',BANDdepth.(BANDS{bb}).pup.(WHNWH{ww}).mean,lfp.interpdepth');
    
    for pp= 1:2
    [ ~,BANDdepth.(BANDS{bb}).(HILO{pp}).(WHNWH{ww}) ] = bz_LFPSpecToExternalVar(...
        lfp.(BANDS{bb}).amp(lfp.(BANDS{bb}).(WHNWH{ww})&lfp.(BANDS{bb}).(HILO{pp}),:),...
        lfp.(BANDS{bb}).pupphase(lfp.(BANDS{bb}).(WHNWH{ww})&lfp.(BANDS{bb}).(HILO{pp}),:),'specparms','input',...
        'figparms',[],'numvarbins',25,'varlim',[-pi pi]);
    BANDdepth.(BANDS{bb}).(HILO{pp}).(WHNWH{ww}).mean_interp = interp1(lfp.depth',BANDdepth.(BANDS{bb}).(HILO{pp}).(WHNWH{ww}).mean,lfp.interpdepth');

    end
end

for oo = 1:2
    for ll = 1:2
        [ ~,BANDdepth.(BANDS{bb}).(ONOFF{oo}).(LONGSHORT{ll}) ] = bz_LFPSpecToExternalVar(...
            lfp.(BANDS{bb}).amp(lfp.(BANDS{bb}).(LONGSHORT{ll}).(ONOFF{oo}),:),...
            lfp.(BANDS{bb}).whtime.(ONOFF{oo})(lfp.(BANDS{bb}).(LONGSHORT{ll}).(ONOFF{oo}),:),'specparms','input',...
            'figparms',[],'numvarbins',401,'varlim',[-window window]);
            BANDdepth.(BANDS{bb}).(ONOFF{oo}).(LONGSHORT{ll}).mean_interp = interp1(lfp.depth',BANDdepth.(BANDS{bb}).(ONOFF{oo}).(LONGSHORT{ll}).mean,lfp.interpdepth');
    end
    
    [ ~,BANDdepth.(BANDS{bb}).(ONOFF{oo}).all ] = bz_LFPSpecToExternalVar(...
        lfp.(BANDS{bb}).amp,...
        lfp.(BANDS{bb}).whtime.(ONOFF{oo}),'specparms','input',...
        'figparms',[],'numvarbins',401,'varlim',[-window window]);
        BANDdepth.(BANDS{bb}).(ONOFF{oo}).all.mean_interp = interp1(lfp.depth',BANDdepth.(BANDS{bb}).(ONOFF{oo}).all.mean,lfp.interpdepth');
end

    
[ ~,BANDdepth.(BANDS{bb}).EMG ] = bz_LFPSpecToExternalVar( lfp.(BANDS{bb}).amp,...
    log10(lfp.(BANDS{bb}).EMG),'specparms','input',...
    'figparms',[],'numvarbins',40,'varlim',[-1.7 0.9]);
    BANDdepth.(BANDS{bb}).EMG.mean_interp = interp1(lfp.depth',BANDdepth.(BANDS{bb}).EMG.mean,lfp.interpdepth');

    
    
    
end    
    
    
    
%%
cosx = linspace(-pi,pi,100);
cospamp = [0.025 0.3];

bandranges = [[1.75 2.3];[1.8 2.2];[2.1 2.8];[2.1 2.8];[2.2 2.8]];


figure
for bb = 1:length(BANDS)
for ww = 1:2
    subplot(5,4,(bb-1)*4+(ww-1)+1)
        for pp = 1:2
        imagesc( BANDdepth.(BANDS{bb}).(HILO{pp}).(WHNWH{ww}).varbins+2*pi*(pp-1),...
            lfp.interpdepth,...
            BANDdepth.(BANDS{bb}).(HILO{pp}).(WHNWH{ww}).mean_interp)
        hold on; axis xy; box off
        plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp)-1,'k')
        end   
        ColorbarWithAxis(bandranges(bb,:),'Power (dB)')
        colorbar
        xlim([-pi 3*pi])
        xlabel('Pupil Phase');
        if ww == 1
            ylabel({BANDS{bb},'Depth'})
        end
        if bb == 1
            title(WHNWH{ww})
        end
        
        
    subplot(5,4,(bb-1)*4+(ww-1)+3)
        imagesc( BANDdepth.(BANDS{bb}).pup.(WHNWH{ww}).varbins,...
            lfp.interpdepth,...
            BANDdepth.(BANDS{bb}).pup.(WHNWH{ww}).mean_interp)
        hold on; axis xy; box off
        ColorbarWithAxis(bandranges(bb,:),'Power (dB)')
        colorbar
        xlabel('Pupil Size');ylabel('Depth')
        
        if bb == 1
            title(WHNWH{ww})
        end

        
end

end
NiceSave('FBandsPupil',figfolder,baseName)
%%

figure
for bb = 1:length(BANDS)
for oo = 1:2

subplot(5,4,(bb-1)*4+oo)
        imagesc( BANDdepth.(BANDS{bb}).(ONOFF{oo}).all.varbins,...
            BANDdepth.(BANDS{bb}).depth,...
            BANDdepth.(BANDS{bb}).(ONOFF{oo}).all.mean_interp)
        hold on; axis xy; box off
        plot([0 0],[-1 0],'w')
       ColorbarWithAxis(bandranges(bb,:),'Power (dB)')
        xlabel(['t - ',(ONOFF{oo})]);
        if oo == 1   
            ylabel({(BANDS{bb}),'Depth'})
        end
        xlim([-1 1])

end
 
subplot(5,4,(bb-1)*4+3)
        imagesc( BANDdepth.(BANDS{bb}).EMG.varbins,...
            BANDdepth.(BANDS{bb}).depth,...
            BANDdepth.(BANDS{bb}).EMG.mean_interp)
        hold on; axis xy; box off
        plot(BANDdepth.(BANDS{bb}).EMG.varbins,BANDdepth.(BANDS{bb}).EMG.vardist*10-1,'k')
        plot(log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],[-1 0],'k--')
       ColorbarWithAxis(bandranges(bb,:),'Power (dB)')
        xlabel('EMG');
        
        
end
NiceSave('FBandsWhisk',figfolder,baseName)


    
%%
% xwinsize = 80;
% xwin = bz_RandomWindowInIntervals(lfp.timestamps([1 end]),xwinsize);
% figure
% subplot(2,1,1)
%     imagesc(lfp.timestamps,lfp.interpdepth,lfp.depthinterp')
%     hold on
%     plot(pupildilation.timestamps,bz_NormToRange(pupildilation.data,[-1 -0.5]),'w','linewidth',2)
%     plot(EMGwhisk.timestamps,bz_NormToRange(EMGwhisk.EMG,[-1.3 -0.8]),'k')
%     xlim(xwin)
%     axis xy
%     ylim([-1.1 0])
%     ColorbarWithAxis([-2.75 -0.5],'PSS')
%     xlabel('t (s)');ylabel('Depth')
%     
%     subplot(2,2,3)
%     imagesc(PSScomponents.depth,PSScomponents.depth,PSScomponents.corr)
%     colorbar
%     axis xy
%     xlabel('Depth');ylabel('Depth')
%     title('PSS Corr')
%     
% %     subplot(2,2,4)
% %         imagesc(PSSdepth.byPup.varbins,PSS.interpdepth,PSSdepth.byPup.mean)
% %         axis xy
% subplot(4,2,6)
%     plot(log10(PSScomponents.EV),'o-')
%     xlim([1 6])
%     xlabel('PC');ylabel('% EV')
%     LogScale('y',10)
% 
% subplot(4,2,8)
%     plot(PSScomponents.depth,PSScomponents.PCAcoeff(:,1:3))
%     hold on
%     plot(PSScomponents.depth([1 end]),[0 0],'k--')
%     legend('PC1','PC2','PC3','location','eastoutside')
%     xlabel('Deptah');ylabel('Weight')
%     
%         NiceSave('DepthPSSEx',figfolder,baseName)


%%

