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

[B,I] = sort(SlowWaves.SWpeakmag,'descend');

eventCSD(lfp,SlowWaves.ints.UP(I(1:round(0.1*length(B))),1),'channels',usechannels,'twin',[0.25 0.25],'spat_sm',8,'saveName','SlowWavesCSD');

% eventCSD(lfp,SlowWaves.ints.UP(:,1),'channels',usechannels,'twin',[0.25 0.25],'spat_sm',8,'saveName','SlowWavesCSD');

% savfile = fullfile(basePath,[baseName,'.EMGwhisk.behavior.mat']);
% load(savfile,'EMGwhisk');
% eventCSD(lfp,EMGwhisk.ints.Wh(:,1),'channels',usechannels,'twin',[0.75 0.75],'spat_sm',8,'saveName','EMGwhiskCSD');
% 
% savfile = fullfile(basePath,[baseName,'.Piezotouch.behavior.mat']);
% if isfile(savfile) > 0
% load(savfile,'Piezotouch');
% eventCSD(lfp,Piezotouch.ints.Touch(:,1),'channels',usechannels,'twin',[0.75 0.75],'spat_sm',8,'saveName','PiezotouchCSD');
% end