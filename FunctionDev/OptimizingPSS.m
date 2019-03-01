basePath = pwd;
baseName = bz_BasenameFromBasepath(basePath);
sessionInfo = bz_getSessionInfo(basePath,'noPrompts',true);
channels = sessionInfo.channels;
usechannels = sessionInfo.AnatGrps.Channels;
badchannels = sessionInfo.badchannels;
badidx = ismember(usechannels,badchannels);
usechannels(badidx) = [];

%Pending cortical channels

%% Load Pupil
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

%% Load the PSS
%load([basePath,filesep,baseName,'.PowerSpectrumSlope.lfp.mat'])

%% Load Whisks
EMGwhisk = bz_LoadBehavior(basePath,'EMGwhisk');

%Specifying SPONT whisking
load(fullfile(basePath,[baseName,'.MergePoints.events.mat']),'MergePoints');
sidx = find(startsWith(MergePoints.foldernames,"Spont"));
sponttimes = [MergePoints.timestamps(sidx(1),1) MergePoints.timestamps(sidx(end),2)];

spontidx = find(EMGwhisk.ints.Wh(:,2) < sponttimes(2));
EMGwhisk.ints.Wh = EMGwhisk.ints.Wh(spontidx,:);

spontidx = find(EMGwhisk.ints.NWh(:,2) < sponttimes(2));
EMGwhisk.ints.NWh = EMGwhisk.ints.NWh(spontidx,:);

spontidx = find(EMGwhisk.timestamps < sponttimes(2));
EMGwhisk.timestamps = EMGwhisk.timestamps(spontidx);
EMGwhisk.EMGenvelope = EMGwhisk.EMGenvelope(spontidx);
EMGwhisk.EMG = EMGwhisk.EMG(spontidx);
EMGwhisk.EMGsm = EMGwhisk.EMGsm(spontidx);

%%
lowfilter = [0.01 0.1];
pupil4filter = pupildilation;

lowpupildata = bz_Filter(pupil4filter,'passband',lowfilter,'filter' ,'fir1','order',3);
%highpupildata = bz_Filter(pupil4filter,'passband',highfilter,'filter' ,'fir1');

%Get the pupil phase/power of each whisk start
EMGwhisk.phase = interp1(lowpupildata.timestamps,lowpupildata.phase,...
    EMGwhisk.ints.Wh(:,1),'nearest');
EMGwhisk.power = interp1(lowpupildata.timestamps,log10(lowpupildata.amp),...
    EMGwhisk.ints.Wh(:,1),'nearest');
%Get the closest PSS timepoint to each whisk... (dt issue...)

%%
lowerbound = [1:20];
upperbound = [60:120];

WhiskPSScorr.EMG = zeros(size(sessionInfo.AnatGrps.Channels));
WhiskPSScorr.pup = zeros(size(sessionInfo.AnatGrps.Channels));
WhiskPSScorr.dpdt = zeros(size(sessionInfo.AnatGrps.Channels));
WhiskPSScorr.phasecoupling = zeros(size(sessionInfo.AnatGrps.Channels));

%%
lfp = bz_GetLFP('all','basepath',basePath,'noPrompts',true);

%%
dt = 0.5;
winsize = 2;
[PSS] = bz_PowerSpectrumSlope(lfp,winsize,dt,'frange',[lowerbound(f) upperbound(ff)],...
    'channels',usechannels,'Redetect',true);

clear lfp

%%
PSS.EMG = interp1(EMGwhisk.timestamps,EMGwhisk.EMGenvelope,PSS.timestamps);
PSS.pupilsize = interp1(pupildilation.timestamps,pupildilation.data,...
    PSS.timestamps,'nearest');
PSS.dpdt = interp1(pupildilation.timestamps(1:end-1),pupildilation.dpdt,...
    PSS.timestamps,'nearest');
PSS.pupilphase = interp1(lowpupildata.timestamps,lowpupildata.phase,...
    PSS.timestamps,'nearest');

%%
[WhiskPSScorr.EMG,WhiskPSScorr.EMG_p] =...
    corr(log10(PSS.EMG),PSS.data,...
    'type','spearman','rows','complete');
[WhiskPSScorr.pup,WhiskPSScorr.pup_p] =...
    corr(log10(PSS.pupilsize),PSS.data,...
    'type','spearman','rows','complete');
[WhiskPSScorr.dpdt,WhiskPSScorr.dpdt_p] =...
    corr(PSS.dpdt,PSS.data,...
    'type','spearman','rows','complete');
WhiskPSScorr.phasecoupling = ...
    abs(nanmean((PSS.data./nanmean(PSS.data)).*exp(1i.*PSS.pupilphase)));

%%
load(fullfile(basePath,[baseName,'.grandWhiskPSScorr.lfp.mat']),'grandWhiskPSScorr');

grandWhiskPSScorr.EMG(f,ff) = max(WhiskPSScorr.EMG);
grandWhiskPSScorr.pup(f,ff) = max(WhiskPSScorr.pup);
grandWhiskPSScorr.dpdt(f,ff) = max(WhiskPSScorr.dpdt);
grandWhiskPSScorr.phasecoupling(f,ff) = max(WhiskPSScorr.phasecoupling);

save(fullfile(basePath,[baseName,'.grandWhiskPSScorr.lfp.mat']),'grandWhiskPSScorr');