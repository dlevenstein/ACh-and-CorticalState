basePath = pwd;
[baseFolder,baseName] = fileparts(basePath);

%% Getting spikes from Kilosort folder and returning to home folder
kilofolder = dir('Kilosort_*');
cd(fullfile(kilofolder.folder,[kilofolder.name]));
spikes = bz_GetSpikes('saveMat',true,'noPrompts',true,'forceReload',true);

cd(basePath);

%%
filterparms.deltafilter = [1 10];%heuristically defined.  room for improvement here.
filterparms.gammafilter = [100 400];
filterparms.gammasmoothwin = 0.08; %window for smoothing gamma power (s)
filterparms.gammanormwin = 20; %window for gamma normalization (s)

%%
load(fullfile(basePath,[baseName,'.sessionInfo.mat']));
dtchan = sessionInfo.slowdetectorchannel;

%%
mergefile = fullfile(basePath,[baseName,'.MergePoints.events.mat']);
load(mergefile,'MergePoints');

spontidx = find(startsWith(MergePoints.foldernames,"Spont"));
sponttimes = [MergePoints.timestamps(spontidx(1),1) MergePoints.timestamps(spontidx(end),2)];

% THIS version w/o Kilosort spikes...
% DetectSlowWaves(basePath,'NREMInts',sponttimes,'CTXChans',dtchan,'filterparms',...
%     filterparms,'noSpikes',true,'showFig',true,'forceReload',true,'noPrompts',true);

DetectSlowWaves(basePath,'NREMInts',sponttimes,'CTXChans',dtchan,'filterparms',...
    filterparms,'spikes',spikes,'showFig',true,'forceReload',true,'noPrompts',true);
