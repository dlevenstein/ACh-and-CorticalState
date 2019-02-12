basePath = pwd;
[baseFolder,baseName] = fileparts(basePath);

%%
filterparms.deltafilter = [1 10];%heuristically defined.  room for improvement here.
%filterparms.deltafilter = [0.5 10];
filterparms.gammafilter = [100 400];
filterparms.gammasmoothwin = 0.08; %window for smoothing gamma power (s)
filterparms.gammanormwin = 20; %window for gamma normalization (s)

%%
load(fullfile(basePath,[baseName,'.sessionInfo.mat']));
dtchan = sessionInfo.slowdetectorchannel;
L5chans = [0 1 30 16 2 29 6 4 36 59 44 34 61 45 33 47 17 63 32 60];

%%
mergefile = fullfile(basePath,[baseName,'.MergePoints.events.mat']);
load(mergefile,'MergePoints');

spontidx = find(startsWith(MergePoints.foldernames,"Spont"));
sponttimes = [MergePoints.timestamps(spontidx(1),1) MergePoints.timestamps(spontidx(end),2)];

DetectSlowWaves(basePath,'NREMInts',sponttimes,'DetectionChannel',2,'filterparms',filterparms,'MUAspikes',true,'forceReload',true,'noPrompts',true);
%DetectSlowWaves(basePath,'NREMInts',sponttimes,'CTXChans',L5chans,'filterparms',filterparms,'MUAspikes',true,'noSpikes',true,'forceReload',true,'noPrompts',true);