basePath = pwd;
[baseFolder,baseName] = fileparts(basePath);

%%
load(fullfile(basePath,[baseName,'.sessionInfo.mat']));
dtchan = sessionInfo.slowdetectorchannel;

%%
mergefile = fullfile(basePath,[baseName,'.MergePoints.events.mat']);
load(mergefile,'MergePoints');

spontidx = find(startsWith(MergePoints.foldernames,"Spont"));
sponttimes = [MergePoints.timestamps(spontidx(1),1) MergePoints.timestamps(spontidx(end),2)];

DetectSlowWaves(basePath,'NREMInts',sponttimes,'DetectionChannel',dtchan,'MUAspikes',true,'noPrompts',true);