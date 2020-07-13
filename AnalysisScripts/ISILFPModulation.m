function [ISILFPMap] = ISILFPModulation(basePath,figfolder)
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
%reporoot = '/Users/dl2820/Project Repos/ACh-and-CorticalState/';
%reporoot = '/gpfs/data/buzsakilab/DL/ACh-and-CorticalState/';
%basePath = '/mnt/proraidDL/Database/WMData/AChPupil/171209_WT_EM1M3/';
%basePath = '/gpfs/data/rudylab/William/171209_WT_EM1M3';
%basePath = '/mnt/proraidDL/Database/WMData/AChPupil/180706_WT_EM1M3/';
%basePath = '/Users/dl2820/Dropbox/research/Datasets/WMProbeData/171209_WT_EM1M3';
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

%%
%%
%EMGwhisk.ints.Wh = EMGwhisk.ints.Wh();
%EMGwhisk.ints.NWh = EMGwhisk.ints.NWh(InIntervals(EMGwhisk.ints.NWh,sponttimes));
% buffer = 0.5;
% EMGwhisk.ints.ExpWh = bsxfun(@plus,EMGwhisk.ints.Wh,[-1 1]*buffer);
% EMGwhisk.ints.Wh = EMGwhisk.ints.Wh(InIntervals(EMGwhisk.ints.Wh,sponttimes));
%% Get the depth info

depthinfo = rescaleCx(basePath);
inCTX = find(~isnan(depthinfo.ndepth));
CTXchans = depthinfo.channels(inCTX);
CTXdepth = -depthinfo.ndepth(inCTX);


%For Piloting
%CTXchans = CTXchans([1:5]);
%CTXdepth = CTXdepth([1:5]);
%%
MapIntNames = {'Spont','NWh','Wh'};
MapInts.Spont = sponttimes;
MapInts.Wh = EMGwhisk.ints.Wh;
MapInts.NWh =EMGwhisk.ints.NWh;
%%
[ISILFPMap] = bz_ISILFPMap(basePath,'groups',{CTXchans},'ints',MapInts,...
    'figfolder',figfolder,'nfreqs',200,'dt',0.25,'forceRedetect',false);

ISILFPMap.MapIntNames = MapIntNames;
%% Now: interpolate to depth

ISILFPMap.interpdepth = linspace(-1,0,100);

ISILFPMap.interp.Wh = interp1(CTXdepth',ISILFPMap.NA.Wh.AllCells(:,ISILFPMap.NA.SGorder)',ISILFPMap.interpdepth');
ISILFPMap.interp.NWh = interp1(CTXdepth',ISILFPMap.NA.NWh.AllCells(:,ISILFPMap.NA.SGorder)',ISILFPMap.interpdepth');
ISILFPMap.interp.AllTime = interp1(CTXdepth',ISILFPMap.NA.Spont.AllCells(:,ISILFPMap.NA.SGorder)',ISILFPMap.interpdepth');

%%

spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
depthinfo_allchans = rescaleCx(basePath,'BADOUT',false);
for cc = 1:spikes.numcells
    spikes.layer(cc) = depthinfo_allchans.layer(depthinfo_allchans.channels==spikes.maxWaveformCh(cc));
    ISILFPMap.NA.layer(cc) = spikes.layer(cc);
end

%%
LAYERS = depthinfo.lnames;

for ii = 1:length(MapIntNames)
for ll = 1:length(LAYERS)
    layercells = strcmp(LAYERS{ll},spikes.layer);
    Layermap.(MapIntNames{ii})(:,:,ll) = nanmedian(ISILFPMap.NA.(MapIntNames{ii}).AllUnits(:,:,layercells),3);
end
end
%[~,~,spikes.layer] = intersect(depthinfo.channels,spikes.maxWaveformCh);
%%
figure
for ii = 1:length(MapIntNames)
    for ll = 1:length(LAYERS)
        subplot(6,3,ii+(ll-1)*3)
            imagesc(log10(ISILFPMap.freqs),[-1 0],Layermap.(MapIntNames{ii})(:,ISILFPMap.NA.SGorder,ll)')
            
            LogScale('x',10)
            if ll == 1
                title(MapIntNames{ii})
            elseif ll == 6
                xlabel('freq (Hz)')
            end
            
            if ii == 1
                ylabel(LAYERS{ll})
            end
    end 
end
NiceSave('LayerISIMod',figfolder,baseName)
%spikes.layer = depthinfo.layer(spikes.layer);
%% Plot with layer dividers
% figure
% subplot(2,3,1)
% imagesc(ISILFPMap.freqs,ISILFPMap.interpdepth,ISILFPMap.interp.Wh)
% axis xy
% subplot(2,3,2)
% imagesc(ISILFPMap.freqs,ISILFPMap.interpdepth,ISILFPMap.interp.NWh)
% axis xy
% subplot(2,3,3)
% imagesc(ISILFPMap.freqs,ISILFPMap.interpdepth,ISILFPMap.interp.AllTime)
% axis xy
% 
% 
% subplot(2,3,4)
% imagesc(ISILFPMap.freqs,ISILFPMap.interpdepth,ISILFPMap.NA.Wh.AllCells(:,ISILFPMap.NA.SGorder)')
% subplot(2,3,5)
% imagesc(ISILFPMap.freqs,ISILFPMap.interpdepth,ISILFPMap.NA.NWh.AllCells(:,ISILFPMap.NA.SGorder)')
% subplot(2,3,6)
% imagesc(ISILFPMap.freqs,ISILFPMap.interpdepth,ISILFPMap.NA.Spont.AllCells(:,ISILFPMap.NA.SGorder)')
% 





