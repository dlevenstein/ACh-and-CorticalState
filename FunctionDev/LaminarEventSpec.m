basePath = pwd;
[baseFolder,baseName] = fileparts(basePath);
%savefile = fullfile(basePath,[baseName,'.',num2str(i),'.LayerID.lfp.mat']);

%%
load(fullfile(basePath,[baseName,'.sessionInfo.mat']),'sessionInfo');
badchannels = sessionInfo.badchannels;
usechannels = sessionInfo.AnatGrps.Channels;
usechannels(ismember(usechannels,badchannels))=[];

savfile = fullfile(basePath,[baseName,'.SlowWaves.events.mat']);
load(savfile,'SlowWaves');
eventSpec(SlowWaves.ints.UP(:,1),'channels',usechannels,'twin',[0.25 0.25],'spat_sm',8,'saveName','SlowWavesCSD');

savfile = fullfile(basePath,[baseName,'.EMGwhisk.behavior.mat']);
load(savfile,'EMGwhisk');
eventSpec(EMGwhisk.ints.Wh(:,1),'channels',usechannels,'twin',[0.75 0.75],'spat_sm',8,'saveName','EMGwhiskCSD');

savfile = fullfile(basePath,[baseName,'.Piezotouch.behavior.mat']);
if isfile(savfile) > 0
    load(savfile,'Piezotouch');
    eventSpec(Piezotouch.ints.Touch(:,1),'channels',usechannels,'twin',[0.75 0.75],'spat_sm',8,'saveName','PiezotouchCSD');
end