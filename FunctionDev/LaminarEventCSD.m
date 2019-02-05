basePath = pwd;
[baseFolder,baseName] = fileparts(basePath);

savfile = fullfile(basePath,[baseName,'.sessionInfo.mat']);
load(savfile,'sessionInfo');
badchannels = sessionInfo.badchannels;
usechannels = sessionInfo.AnatGrps.Channels;
usechannels(ismember(usechannels,badchannels))=[];

lfp = bz_GetLFP('all','noPrompts',true);

savfile = fullfile(basePath,[baseName,'.SlowWaves.events.mat']);
load(savfile,'SlowWaves');

bz_eventCSD(lfp,SlowWaves.ints.UP(:,1),'channels',usechannels,'twin',[0.25 0.25],'spat_sm',8,'saveName','SlowWavesCSD');