function [PSSConditionalISI,ConditionalRate] = PSSandSpikesAnalysis(basePath,figfolder)
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
%reporoot = '/Users/dl2820/Project Repos/ACh-and-CorticalState/';
%basePath = '/mnt/proraidDL/Database/WMData/AChPupil/171209_WT_EM1M3/';
%basePath = '/mnt/proraidDL/Database/WMData/AChPupil/180706_WT_EM1M3/';
%basePath = '/Users/dl2820/Dropbox/research/Datasets/WMProbeData/171209_WT_EM1M3';
%basePath = pwd;
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath,'noPrompts',true);


CellClass = bz_LoadCellinfo(basePath,'CellClass');


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


%For Piloting
%CTXchans = CTXchans([1:5]);
%CTXdepth = CTXdepth([1:5]);


%% Calculate Spectrogram on all channels
clear spec
% %dt = 0.1;
% spec.winsize = 1;
 spec.frange = [2 128]; %Frequency lower than can be assessed for window because IRASA... but maybe this is bad for IRASA too
 spec.nfreqs = 150;
% 
% % noverlap = spec.winsize-dt;
% % spec.freqs = logspace(log10(spec.frange(1)),log10(spec.frange(2)),spec.nfreqs);
% % winsize_sf = round(spec.winsize .*lfp.samplingRate);
% % noverlap_sf = round(noverlap.*lfp.samplingRate);
% 
% spec.channels = lfp.channels;
spec.channels = CTXchans;
spec.depth = CTXdepth;
ncycles = 10; %prev 10
for cc =1:length(spec.channels)
    bz_Counter(cc,length(spec.channels),'Channel')
    
    specslope = bz_PowerSpectrumSlope([],ncycles,0.01,'channels',spec.channels(cc),...
        'frange',spec.frange,'spectype','wavelet','nfreqs',spec.nfreqs,'ints',sponttimes,...
        'saveMat',basePath,'saveName',['wav',num2str(spec.channels(cc))],...
        'Redetect',false,'suppressText',true);
    
    try
    spec.PSS(:,cc) = specslope.data; %Issue here: Unable to perform assignment because the size of the left side is 809294-by-1
    spec.timestamps = specslope.timestamps;
    catch
        display(['Issue: channel ',num2str(spec.channels(cc))])
        spec.PSS(:,cc) = nan(size(spec.timestamps));
    end
    spec.freqs = specslope.freqs; 
    clear specslope
end
%%
spec.winsize = 1;
spec.chanlayers = depthinfo.layer(inCTX);
%% Take Mean Specgram by layer and calculate irasa

LAYERS = depthinfo.lnames;

for ll = 1:length(LAYERS)
    layerchans = strcmp(LAYERS{ll},spec.chanlayers);
    repchan(ll) = round(median(find(layerchans)));
    spec.LayerPSS(:,ll) = nanmean(spec.PSS(:,layerchans),2);
    % Median Normalize Spectrogram
    %spec.Layer(:,:,ll) = NormToInt(spec.Layer(:,:,ll),'median');
end
%spec = rmfield(spec,'data');


%%

spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
depthinfo_allchans = rescaleCx(basePath,'BADOUT',false);
for cc = 1:spikes.numcells
    spikes.layer(cc) = depthinfo_allchans.layer(depthinfo_allchans.channels==spikes.maxWaveformCh(cc));
end


%% Spike Rate
dt = 0.01;
binsize = 0.06;

spkmat = bz_SpktToSpkmat(spikes.times,'dt',dt,'binsize',binsize,...
    'win',sponttimes,'units','rate','bintype','gaussian');
spkmat.NWh = InIntervals(spkmat.timestamps,EMGwhisk.ints.NWh);

spkmat.poprate.All = mean(spkmat.data(:,CellClass.pE),2);
for ll = 1:length(LAYERS)
    layercells = strcmp(LAYERS{ll},spikes.layer & CellClass.pE);
    spkmat.poprate.(LAYERS{ll}) = mean(spkmat.data(:,layercells),2);
end
for ll = 1:length(LAYERS)
    spkmat.PSS.(LAYERS{ll}) = interp1(spec.timestamps,spec.LayerPSS(:,ll),spkmat.timestamps,'nearest');
    %spkmat.PSS.(LAYERS{ll}) = NormToInt(spkmat.PSS.(LAYERS{ll}),'percentile')
end
%%
minX = 100;
for sll = 1:length(LAYERS)
    for pll = 1:length(LAYERS)
[ ConditionalRate.(LAYERS{pll}).(LAYERS{sll})] = ConditionalHist(spkmat.PSS.(LAYERS{pll}),spkmat.poprate.(LAYERS{sll}),...
        'Xbounds',[-1.75 0],'numXbins',40,'Ybounds',[0 20],'numYbins',40,'minX',minX);
    
[ ConditionalRate.log.(LAYERS{pll}).(LAYERS{sll})] = ConditionalHist(spkmat.PSS.(LAYERS{pll}),log10(spkmat.poprate.(LAYERS{sll})),...
        'Xbounds',[-1.75 0],'numXbins',40,'Ybounds',[-0.5 1.5],'numYbins',40,'minX',minX);
    
    
[ ConditionalRate.NWh.(LAYERS{pll}).(LAYERS{sll})] = ConditionalHist(spkmat.PSS.(LAYERS{pll})(spkmat.NWh),...
    log10(spkmat.poprate.(LAYERS{sll})(spkmat.NWh)),...
        'Xbounds',[-1.75 0],'numXbins',40,'Ybounds',[-0.5 1.5],'numYbins',40,'minX',minX);
    %Ints: NWh/Wh
    %Pupil/Whisk onset? (other script)
    end
end
    %%
figure
for sll = 1:length(LAYERS)
    for pll = 1:length(LAYERS)
        subplot(6,6,pll+(sll-1)*6)
    imagesc(ConditionalRate.log.(LAYERS{pll}).(LAYERS{sll}).Xbins,...
        ConditionalRate.log.(LAYERS{pll}).(LAYERS{sll}).Ybins,...
        ConditionalRate.log.(LAYERS{pll}).(LAYERS{sll}).pYX')
    alpha(gca,single(ConditionalRate.log.(LAYERS{pll}).(LAYERS{sll}).pYX'>1e-4))
    axis xy
    clim([0 0.1])
    colorbar
    crameri turku
    LogScale('y',10)
    if sll==6
        xlabel(['PSS, ',(LAYERS{pll})])
    end
    if pll==1
        ylabel({(LAYERS{sll}),'Pop Rate (Hz)'})
    end
    
    end
end
NiceSave('PopRateandPSS',figfolder,baseName)
%%
clear PSSConditionalISI
for pll = 1:length(LAYERS)
    bz_Counter(pll,length(LAYERS),'PSS Layer')
    PSS.timestamps = spec.timestamps;
    PSS.data = spec.LayerPSS(:,pll);
    [PSSConditionalISI.(LAYERS{pll})] = bz_ConditionalISI(spikes.times,PSS,...
        'showfig',false,'ISIDist',true);%,'ints',ints.(states{ss}));%,...
    
    for ll = 1:length(LAYERS)
        layercells = strcmp(LAYERS{ll},spikes.layer) & CellClass.pE;
        layermean.(LAYERS{pll}).(LAYERS{ll}) = nanmean(PSSConditionalISI.(LAYERS{pll}).Dist.pYX(:,:,layercells),3);
    end
    
end
PSSConditionalISI.celllayers = spikes.layer;
PSSConditionalISI.CellClass = CellClass;
%%
figure
for sll = 1:length(LAYERS)
    for pll = 1:length(LAYERS)
        subplot(6,6,pll+(sll-1)*6)
    imagesc(PSSConditionalISI.(LAYERS{pll}).Dist.Xbins(1,:,1),PSSConditionalISI.(LAYERS{pll}).Dist.Ybins(1,:,1),layermean.(LAYERS{pll}).(LAYERS{sll})')
    axis xy
    %crameri bilbao
    if sll==6
        xlabel(['PSS, ',(LAYERS{pll})])
    end
    if pll==1
        ylabel({(LAYERS{sll}),'ISI'})
    end
    
    LogScale('y',10,'exp',10)
    
    end
end
NiceSave('ISIandPSS',figfolder,baseName)
%%
% for cc = 1:length(spec.channels)
%     bz_Counter(cc,length(spec.channels),'channel')
% PSS.timestamps = spec.timestamps;
% PSS.data = spec.PSS(:,cc);
% [PSSConditionalISI(cc)] = bz_ConditionalISI(spikes.times,PSS,...
%     'showfig',false,'ISIDist',true);%,'ints',ints.(states{ss}));%,...
% MI(cc,:) = squeeze(PSSConditionalISI(cc).MutInf);
% end

