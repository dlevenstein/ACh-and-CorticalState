basePath = pwd;
[baseFolder,baseName] = fileparts(basePath);

savfile = fullfile(basePath,[baseName,'.sessionInfo.mat']);
load(savfile,'sessionInfo');
badchannels = sessionInfo.badchannels;
usechannels = sessionInfo.AnatGrps.Channels;
usechannels(ismember(usechannels,badchannels))=[];

mua = [];

savfile = fullfile(basePath,[baseName,'.SlowWaves.events.mat']);
load(savfile,'SlowWaves');
eventMUA(mua,SlowWaves.ints.UP(:,1),'channels',usechannels,'twin',[0.25 0.25],'spat_sm',3,'saveName','SlowWavesMUA');

savfile = fullfile(basePath,[baseName,'.EMGwhisk.behavior.mat']);
load(savfile,'EMGwhisk');
eventMUA(mua,EMGwhisk.ints.Wh(:,1),'channels',usechannels,'twin',[0.75 0.75],'spat_sm',3,'saveName','EMGwhiskMUA');

savfile = fullfile(basePath,[baseName,'.Piezotouch.behavior.mat']);
if isfile(savfile) > 0
load(savfile,'Piezotouch');
eventMUA(mua,Piezotouch.ints.Touch(:,1),'channels',usechannels,'twin',[0.75 0.75],'spat_sm',3,'saveName','PiezotouchMUA');
end