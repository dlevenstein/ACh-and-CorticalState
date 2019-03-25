basePath = pwd;
baseName = bz_BasenameFromBasepath(basePath);
sessionInfo = bz_getSessionInfo(basePath,'noPrompts',true);

% badchannels = sessionInfo.badchannels;
% usechannels = sessionInfo.AnatGrps.Channels;
% usechannels(ismember(usechannels,badchannels))=[];
channels = sessionInfo.channels;

% Obtaining normalized layer info
rescaled = rescaleCx(basePath);
usechannels = rescaled.channels(~isnan(rescaled.ndepth));
normdepth = rescaled.ndepth(~isnan(rescaled.ndepth));

L1idx = find(rescaled.ndepth >= 0 & rescaled.ndepth <= 0.1);
L1idx = rescaled.channels(L1idx);

L23idx = find(rescaled.ndepth >= 0.1 & rescaled.ndepth <= 0.35);
L23idx = rescaled.channels(L23idx);

L4idx = find(rescaled.ndepth >= 0.35 & rescaled.ndepth <= 0.5);
L4idx = rescaled.channels(L4idx);

L5aidx = find(rescaled.ndepth >= 0.5 & rescaled.ndepth <= 0.6);
L5aidx = rescaled.channels(L5aidx);

L56idx = find(rescaled.ndepth >= 0.6 & rescaled.ndepth <= 0.9);
L56idx = rescaled.channels(L56idx);

L6idx = find(rescaled.ndepth >= 0.9 & rescaled.ndepth <= 1);
L6idx = rescaled.channels(L6idx);

% Specifying folders
figfolder = fullfile(basePath,'AnalysisFigures');
savefile = fullfile(basePath,[baseName,'.LaminarSpectralAnalysis.mat']);
savefolder = fullfile(basePath,'WaveSpec');
%savefolder = fullfile(basePath,'WaveSpec2');

%% Loading behavior...
% Pupil diameter
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

% EMG
EMGwhisk = bz_LoadBehavior(basePath,'EMGwhisk');

% Specifying SPONT whisking
load(fullfile(basePath,[baseName,'.MergePoints.events.mat']),'MergePoints');
sidx = find(startsWith(MergePoints.foldernames,"Spont"));
%sponttimes = [MergePoints.timestamps(sidx(1),1) MergePoints.timestamps(sidx(end),2)];
sponttimes = [MergePoints.timestamps(sidx(1),1) MergePoints.timestamps(sidx(1),2)/8];

spontidx = find(EMGwhisk.ints.Wh(:,2) < sponttimes(2));
EMGwhisk.ints.Wh = EMGwhisk.ints.Wh(spontidx,:);

spontidx = find(EMGwhisk.ints.NWh(:,2) < sponttimes(2));
EMGwhisk.ints.NWh = EMGwhisk.ints.NWh(spontidx,:);

spontidx = find(EMGwhisk.timestamps < sponttimes(2));
EMGwhisk.timestamps = EMGwhisk.timestamps(spontidx);
EMGwhisk.EMGenvelope = EMGwhisk.EMGenvelope(spontidx);
EMGwhisk.EMG = EMGwhisk.EMG(spontidx);
EMGwhisk.EMGsm = EMGwhisk.EMGsm(spontidx);

% Tentative for troubleshooting
spontidx = find(pupildilation.timestamps < sponttimes(2));
pupildilation.data = pupildilation.data(spontidx);
pupildilation.dpdt = pupildilation.dpdt(spontidx);
pupildilation.timestamps = pupildilation.timestamps(spontidx);

%% Partition lo/hi Pupil dilation indices
% Getting intervals in spec times
load(fullfile(savefolder,[baseName,'.',num2str(0),'.WaveSpec.lfp.mat']));
%load(fullfile(savefolder,[baseName,'.',num2str(0),'.WaveSpec2.lfp.mat']));

% Pupil phase
% lowfilter = [0.01 0.1];
% pupil4filter = pupildilation;
% lowpupildata = bz_Filter(pupil4filter,'passband',lowfilter,'filter' ,'fir1','order',3);
% x = lowpupildata.timestamps;
% x = repmat(x,[1 5]);
% y = repmat((0:4) * (1/(lowpupildata.samplingRate*5)) ,[length(x) 1]);
% x = x + y;
% y = reshape(x',[length(x)*5 1]);
% lowpupildata.timestamps = y;
% lowpupildata.amp = resample(lowpupildata.amp,5,1);
%
% pupthresh = nanmedian(log10(lowpupildata.amp));
% highpup = log10(lowpupildata.amp)>pupthresh;

pupthresh = nanmedian(pupildilation.data);
highpup = pupildilation.data>pupthresh;

% eventshipupil = interp1(wavespec.timestamps,...
%     wavespec.timestamps,...
%     lowpupildata.timestamps(highpup),'nearest').*wavespec.samplingRate;
%
% eventslopupil = interp1(wavespec.timestamps,...
%     wavespec.timestamps,...
%     lowpupildata.timestamps(~highpup),'nearest').*wavespec.samplingRate;

allidx.hiP = interp1(wavespec.timestamps,...
    wavespec.timestamps,...
    pupildilation.timestamps(highpup),'nearest');
for i = 1:length(allidx.hiP)
    allidx.hiP(i) = find(wavespec.timestamps == allidx.hiP(i)); 
end

allidx.loP = interp1(wavespec.timestamps,...
    wavespec.timestamps,...
    pupildilation.timestamps(~highpup),'nearest');
for i = 1:length(allidx.loP)
    allidx.loP(i) = find(wavespec.timestamps == allidx.loP(i)); 
end

%% PUPIL dilation ONsets/OFFsets
% Set the Pupil dilation thresholds by troughs 
tempPup = pupildilation.dpdt;
tempPup(tempPup<0) = nan;
%Pupz = NormToInt(tempPup,'modZ'); %Modified Z score - robust to outliers
Pupz = (tempPup - nanmean(tempPup))./nanstd(tempPup,0,1);

% find by "gradient descent"(ish) from initial guess (0.5)
tPupbins = linspace(-1.5,2,100);
tPuphist = hist(log10(Pupz),tPupbins);
Pupgrad = smooth(gradient(tPuphist),4);

% Find troughs (gradient crossing from - to +)
troughidx = find(diff(Pupgrad<0)==1);
troughs = 10.^tPupbins(troughidx);

pupthreshold = troughs(find(troughs>0.5,1,'first'));
pup_thresh = Pupz > pupthreshold;
pup_on = round(find(pup_thresh(2:end) > pup_thresh(1:end-1))+1); %Pup onsets (si)
pup_off = round(find(pup_thresh(2:end) < pup_thresh(1:end-1))+1); %Pup offsets (si)

% If data starts/ends in the middle of an epoch, drop first/last trigger
if pup_off(1)<pup_on(1)
    pup_off = pup_off(2:end);
end
if pup_off(end) < pup_on(end)
    pup_on = pup_on(1:end-1);
end

for i = 1:length(pup_on)
    tidx = find(pupildilation.dpdt(pup_on(i):pup_off(i))...
        == max(pupildilation.dpdt(pup_on(i):pup_off(i))),1,'first');
    pup_on(i) = pup_on(i)+tidx;
end

% figure; 
% plot(pupildilation.dpdt,'k'); hold on;
% plot(pup_on,pupildilation.dpdt(pup_on),'ro')

%temppup_on = pup_on;
for i = 1:length(pup_on)
    tidx = find(pupildilation.dpdt(1:pup_on(i)) < 0,1,'last');
    pup_on(i) = tidx;
end

% figure; 
% plot(pupildilation.dpdt,'k'); hold on;
% plot(temppup_on,pupildilation.dpdt(temppup_on),'ro')
% plot(pup_on,pupildilation.dpdt(pup_on),'go')

for i = 1:length(pup_on)
    tidx = find(pupildilation.data(pup_on(i):end-1) < pupildilation.data(pup_on(i)),1,'first');
    if ~isempty(tidx)
        pup_off(i) = pup_on(i)+tidx;
    else
        pup_off(i) = length(pupildilation.data)-1;
    end
end

% Convert to seconds, merge and exclude transients
pup_on = pupildilation.timestamps(pup_on);
pup_off = pupildilation.timestamps(pup_off);
% [ Pupints ] = MergeSeparatedInts( [pup_on,pup_off],0.1 );
Pupints = [pup_on,pup_off];
% [pup_on,pup_off] = MinEpochLength(Pupints(:,1),Pupints(:,2),0.1,1);
% Pupints = [pup_on,pup_off];

eventsidx.Pup = interp1(wavespec.timestamps,wavespec.timestamps,...
    pup_on,'nearest');
for i = 1:length(eventsidx.Pup)
    eventsidx.Pup(i) = find(wavespec.timestamps == eventsidx.Pup(i)); 
end

% Durations and peaks
Pupdur = pup_off-pup_on;
pup_peak = length(pup_on);
for i = 1:length(pup_on)
    tsidx = find(pupildilation.timestamps == pup_on(i));
    teidx = find(pupildilation.timestamps == pup_off(i));
    pup_peak(i) = max(pupildilation.data(tsidx:teidx));
end

%% Partition lo/hi, short/long Whisking epochs/indices
% For NWh
nwhdurs = EMGwhisk.ints.NWh(:,2)-EMGwhisk.ints.NWh(:,1);

eventsidx.NWh = interp1(wavespec.timestamps,wavespec.timestamps,...
    EMGwhisk.ints.NWh,'nearest');
allidx.NWh = [];
for e = 1:size(eventsidx.NWh,1)
    tempidx = find([wavespec.timestamps >= eventsidx.NWh(e,1) & wavespec.timestamps <= eventsidx.NWh(e,2)]);
    allidx.NWh = cat(1,allidx.NWh,tempidx);
end
for i = 1:length(eventsidx.NWh)
    eventsidx.NWh(i) = find(wavespec.timestamps == eventsidx.NWh(i)); 
end

% For >0.5 s NWh
tevents = interp1(wavespec.timestamps,wavespec.timestamps,...
    EMGwhisk.ints.NWh(nwhdurs>=0.5,:),'nearest');
tallidx = [];
for e = 1:size(tevents,1)
    tempidx = find([wavespec.timestamps >= tevents(e,1) & wavespec.timestamps <= tevents(e,2)]);
    tallidx = cat(1,tallidx,tempidx);
end
for i = 1:length(tevents)
    tevents(i) = find(wavespec.timestamps == tevents(i)); 
end
eventsidx.lNWh = tevents;
allidx.lNWh = tallidx;

% For Wh
whpeak = NaN(size(EMGwhisk.ints.Wh,1),1);
for i = 1:size(EMGwhisk.ints.Wh,1)
    tempidx = find(EMGwhisk.timestamps >= EMGwhisk.ints.Wh(i,1)...
        & EMGwhisk.timestamps <= EMGwhisk.ints.Wh(i,2));
    whpeak(i) = nanmean(EMGwhisk.EMGsm(tempidx),1);
end
whthresh = nanmedian(whpeak);
whdurs = EMGwhisk.ints.Wh(:,2)-EMGwhisk.ints.Wh(:,1);

% all Whisking
tevents = interp1(wavespec.timestamps,wavespec.timestamps,...
    EMGwhisk.ints.Wh,'nearest');
tallidx = [];
for e = 1:size(tevents,1)
    tempidx = find([wavespec.timestamps >= tevents(e,1) & wavespec.timestamps <= tevents(e,2)]);
    tallidx = cat(1,tallidx,tempidx);
end
for i = 1:length(tevents)
    tevents(i) = find(wavespec.timestamps == tevents(i)); 
end
eventsidx.Wh = tevents;
allidx.Wh = tallidx;

% <0.5 s Whisking
tevents = interp1(wavespec.timestamps,wavespec.timestamps,...
    EMGwhisk.ints.Wh(whdurs<0.5,:),'nearest');
tallidx = [];
for e = 1:size(tevents,1)
    tempidx = find([wavespec.timestamps >= tevents(e,1) & wavespec.timestamps <= tevents(e,2)]);
    tallidx = cat(1,tallidx,tempidx);
end
for i = 1:length(tevents)
    tevents(i) = find(wavespec.timestamps == tevents(i)); 
end
eventsidx.sWh = tevents;
allidx.sWh = tallidx;

% >0.5 s Whisking
tevents = interp1(wavespec.timestamps,wavespec.timestamps,...
    EMGwhisk.ints.Wh(whdurs>=0.5,:),'nearest');
tallidx = [];
for e = 1:size(tevents,1)
    tempidx = find([wavespec.timestamps >= tevents(e,1) & wavespec.timestamps <= tevents(e,2)]);
    tallidx = cat(1,tallidx,tempidx);
end
for i = 1:length(tevents)
    tevents(i) = find(wavespec.timestamps == tevents(i)); 
end
eventsidx.lWh = tevents;
allidx.lWh = tallidx;

% <median amp Whisking
tevents = interp1(wavespec.timestamps,wavespec.timestamps,...
    EMGwhisk.ints.Wh(whpeak<whthresh,:),'nearest');
tallidx = [];
for e = 1:size(tevents,1)
    tempidx = find([wavespec.timestamps >= tevents(e,1) & wavespec.timestamps <= tevents(e,2)]);
    tallidx = cat(1,tallidx,tempidx);
end
for i = 1:length(tevents)
    tevents(i) = find(wavespec.timestamps == tevents(i)); 
end
eventsidx.loWh = tevents;
allidx.loWh = tallidx;

% >median amp Whisking
tevents = interp1(wavespec.timestamps,wavespec.timestamps,...
    EMGwhisk.ints.Wh(whpeak>=whthresh,:),'nearest');
tallidx = [];
for e = 1:size(tevents,1)
    tempidx = find([wavespec.timestamps >= tevents(e,1) & wavespec.timestamps <= tevents(e,2)]);
    tallidx = cat(1,tallidx,tempidx);
end
for i = 1:length(tevents)
    tevents(i) = find(wavespec.timestamps == tevents(i)); 
end
eventsidx.hiWh = tevents;
allidx.hiWh = tallidx;

% for clear NWh baseline prior to Wh onset
tempidx = [];
for i = 1:size(EMGwhisk.ints.Wh,1)
    idx = find(EMGwhisk.ints.NWh(:,1) > EMGwhisk.ints.Wh(i,1)-5 & EMGwhisk.ints.NWh(:,2) < EMGwhisk.ints.Wh(i,1));
    if ~isempty(idx)
        tempidx = cat(1,tempidx,i);
    else
    end
end
tevents = interp1(wavespec.timestamps,wavespec.timestamps,...
    EMGwhisk.ints.Wh(tempidx,:),'nearest');
for i = 1:length(tevents)
    tevents(i) = find(wavespec.timestamps == tevents(i)); 
end
eventsidx.clearWh = tevents;

%% Loading binary data...
% sessionInfo = bz_getSessionInfo(basePath, 'noPrompts', true);
% datSampleRate = sessionInfo.rates.wideband;
% datfilename = fullfile(basePath,[baseName,'.dat']);
% channels = sessionInfo.channels;
% 
% downfactor = 25;
% datlfp.data = bz_LoadBinary(datfilename,...
%               'frequency',datSampleRate,'nchannels',sessionInfo.nChannels,...
%               'channels',channels+1,'downsample',downfactor);
% 
% datlfp.samplingRate = datSampleRate./downfactor;
% datlfp.timestamps = [0:(length(datlfp.data)-1)]'/datlfp.samplingRate;  %To be overwritten later...
% datlfp.channels = channels;

%% Columnar Power spectra and oscillospecs by state
maxRescaleFactor = 2.9; 
numberRescalesfreq = maxRescaleFactor*wavespec.freqs(1);
numberRescalesidx = find(wavespec.freqs >= numberRescalesfreq);
numberRescales = numberRescalesidx(1)-1;
validFreqInds = numberRescales + 1:wavespec.nfreqs - numberRescales - 1;

% Allocating...
columnspec_all.db = NaN(size(wavespec.data,2),length(channels));
columnspec_all.mednorm = NaN(size(wavespec.data,2),length(channels));
columnspec_all.modz = NaN(size(wavespec.data,2),length(channels));
columnspec_all.frac = NaN(length(validFreqInds),length(channels));
columnspec_all.osci = NaN(length(validFreqInds),length(channels));

columnspec_NWh.db = NaN(size(wavespec.data,2),length(channels));
columnspec_NWh.mednorm = NaN(size(wavespec.data,2),length(channels));
columnspec_NWh.modz = NaN(size(wavespec.data,2),length(channels));
columnspec_NWh.frac = NaN(length(validFreqInds),length(channels));
columnspec_NWh.osci = NaN(length(validFreqInds),length(channels));

columnspec_lNWh.db = NaN(size(wavespec.data,2),length(channels));
columnspec_lNWh.mednorm = NaN(size(wavespec.data,2),length(channels));
columnspec_lNWh.modz = NaN(size(wavespec.data,2),length(channels));
columnspec_lNWh.frac = NaN(length(validFreqInds),length(channels));
columnspec_lNWh.osci = NaN(length(validFreqInds),length(channels));

columnspec_Wh.db = NaN(size(wavespec.data,2),length(channels));
columnspec_Wh.mednorm = NaN(size(wavespec.data,2),length(channels));
columnspec_Wh.modz = NaN(size(wavespec.data,2),length(channels));
columnspec_Wh.frac = NaN(length(validFreqInds),length(channels));
columnspec_Wh.osci = NaN(length(validFreqInds),length(channels));

columnspec_sWh.db = NaN(size(wavespec.data,2),length(channels));
columnspec_sWh.mednorm = NaN(size(wavespec.data,2),length(channels));
columnspec_sWh.modz = NaN(size(wavespec.data,2),length(channels));
columnspec_sWh.frac = NaN(length(validFreqInds),length(channels));
columnspec_sWh.osci = NaN(length(validFreqInds),length(channels));

columnspec_lWh.db = NaN(size(wavespec.data,2),length(channels));
columnspec_lWh.mednorm = NaN(size(wavespec.data,2),length(channels));
columnspec_lWh.modz = NaN(size(wavespec.data,2),length(channels));
columnspec_lWh.frac = NaN(length(validFreqInds),length(channels));
columnspec_lWh.osci = NaN(length(validFreqInds),length(channels));

columnspec_loWh.db = NaN(size(wavespec.data,2),length(channels));
columnspec_loWh.mednorm = NaN(size(wavespec.data,2),length(channels));
columnspec_loWh.modz = NaN(size(wavespec.data,2),length(channels));
columnspec_loWh.frac = NaN(length(validFreqInds),length(channels));
columnspec_loWh.osci = NaN(length(validFreqInds),length(channels));

columnspec_hiWh.db = NaN(size(wavespec.data,2),length(channels));
columnspec_hiWh.mednorm = NaN(size(wavespec.data,2),length(channels));
columnspec_hiWh.modz = NaN(size(wavespec.data,2),length(channels));
columnspec_hiWh.frac = NaN(length(validFreqInds),length(channels));
columnspec_hiWh.osci = NaN(length(validFreqInds),length(channels));

columnspec_loP.db = NaN(size(wavespec.data,2),length(channels));
columnspec_loP.mednorm = NaN(size(wavespec.data,2),length(channels));
columnspec_loP.modz = NaN(size(wavespec.data,2),length(channels));
columnspec_loP.frac = NaN(length(validFreqInds),length(channels));
columnspec_loP.osci = NaN(length(validFreqInds),length(channels));

columnspec_hiP.db = NaN(size(wavespec.data,2),length(channels));
columnspec_hiP.mednorm = NaN(size(wavespec.data,2),length(channels));
columnspec_hiP.modz = NaN(size(wavespec.data,2),length(channels));
columnspec_hiP.frac = NaN(length(validFreqInds),length(channels));
columnspec_hiP.osci = NaN(length(validFreqInds),length(channels));

for i = 1:length(channels)
    i
    % Loading spectrograms
    load(fullfile(savefolder,[baseName,'.',num2str(channels(i)),'.WaveSpec.lfp.mat']));
    %load(fullfile(savefolder,[baseName,'.',num2str(channels(i)),'.WaveSpec2.lfp.mat']));
    
    wavespec.dataz = NormToInt(log10(abs(wavespec.data)),'modZ');
    wavespec.datan = log10(abs(wavespec.data))./nanmedian(log10(abs(wavespec.data)),1);
    
    % WaveIRASA wavelet spec
    [wavespec.frac,wavespec.osci,wavespec.validfreq] = WaveIRASA(wavespec);
    
    % Averaging...
    columnspec_all.db(:,i) = nanmedian(log10(abs(wavespec.data)),1);
    columnspec_all.mednorm(:,i) = nanmedian(wavespec.datan,1);
    columnspec_all.modz(:,i) = nanmedian(wavespec.dataz,1);
    columnspec_all.frac(:,i) = nanmedian(log10(wavespec.frac),1);
    columnspec_all.osci(:,i) = nanmedian(wavespec.osci,1);
    %columnspec_all.osci(columnspec_all.osci(:,i)<0,i) = 0;
    
    columnspec_NWh.db(:,i) = nanmedian(log10(abs(wavespec.data(allidx.NWh,:))),1);
    columnspec_NWh.mednorm(:,i) = nanmedian(wavespec.datan(allidx.NWh,:),1);
    columnspec_NWh.modz(:,i) = nanmedian(wavespec.dataz(allidx.NWh,:),1);
    columnspec_NWh.frac(:,i) = nanmedian(log10(wavespec.frac(allidx.NWh,:)),1);
    columnspec_NWh.osci(:,i) = nanmedian(wavespec.osci(allidx.NWh,:),1);
    %columnspec_NWh.osci(columnspec_NWh.osci(:,i)<0,i) = 0;
    
    columnspec_lNWh.db(:,i) = nanmedian(log10(abs(wavespec.data(allidx.lNWh,:))),1);
    columnspec_lNWh.mednorm(:,i) = nanmedian(wavespec.datan(allidx.lNWh,:),1);
    columnspec_lNWh.modz(:,i) = nanmedian(wavespec.dataz(allidx.lNWh,:),1);
    columnspec_lNWh.frac(:,i) = nanmedian(log10(wavespec.frac(allidx.lNWh,:)),1);
    columnspec_lNWh.osci(:,i) = nanmedian(wavespec.osci(allidx.lNWh,:),1);
    
    columnspec_Wh.db(:,i) = nanmedian(log10(abs(wavespec.data(allidx.Wh,:))),1);
    columnspec_Wh.mednorm(:,i) = nanmedian(wavespec.datan(allidx.Wh,:),1);
    columnspec_Wh.modz(:,i) = nanmedian(wavespec.dataz(allidx.Wh,:),1);
    columnspec_Wh.frac(:,i) = nanmedian(log10(wavespec.frac(allidx.Wh,:)),1);
    columnspec_Wh.osci(:,i) = nanmedian(wavespec.osci(allidx.Wh,:),1);
    %columnspec_Wh.osci(columnspec_Wh.osci(:,i)<0,i) = 0;
    
    columnspec_sWh.db(:,i) = nanmedian(log10(abs(wavespec.data(allidx.sWh,:))),1);
    columnspec_sWh.mednorm(:,i) = nanmedian(wavespec.datan(allidx.sWh,:),1);
    columnspec_sWh.modz(:,i) = nanmedian(wavespec.dataz(allidx.sWh,:),1);
    columnspec_sWh.frac(:,i) = nanmedian(log10(wavespec.frac(allidx.sWh,:)),1);
    columnspec_sWh.osci(:,i) = nanmedian(wavespec.osci(allidx.sWh,:),1);
    
    columnspec_lWh.db(:,i) = nanmedian(log10(abs(wavespec.data(allidx.lWh,:))),1);
    columnspec_lWh.mednorm(:,i) = nanmedian(wavespec.datan(allidx.lWh,:),1);
    columnspec_lWh.modz(:,i) = nanmedian(wavespec.dataz(allidx.lWh,:),1);
    columnspec_lWh.frac(:,i) = nanmedian(log10(wavespec.frac(allidx.lWh,:)),1);
    columnspec_lWh.osci(:,i) = nanmedian(wavespec.osci(allidx.lWh,:),1);
    
    columnspec_loWh.db(:,i) = nanmedian(log10(abs(wavespec.data(allidx.loWh,:))),1);
    columnspec_loWh.mednorm(:,i) = nanmedian(wavespec.datan(allidx.loWh,:),1);
    columnspec_loWh.modz(:,i) = nanmedian(wavespec.dataz(allidx.loWh,:),1);
    columnspec_loWh.frac(:,i) = nanmedian(log10(wavespec.frac(allidx.loWh,:)),1);
    columnspec_loWh.osci(:,i) = nanmedian(wavespec.osci(allidx.loWh,:),1);
    
    columnspec_hiWh.db(:,i) = nanmedian(log10(abs(wavespec.data(allidx.hiWh,:))),1);
    columnspec_hiWh.mednorm(:,i) = nanmedian(wavespec.datan(allidx.hiWh,:),1);
    columnspec_hiWh.modz(:,i) = nanmedian(wavespec.dataz(allidx.hiWh,:),1);
    columnspec_hiWh.frac(:,i) = nanmedian(log10(wavespec.frac(allidx.hiWh,:)),1);
    columnspec_hiWh.osci(:,i) = nanmedian(wavespec.osci(allidx.hiWh,:),1);
    
    columnspec_loP.db(:,i) = nanmedian(log10(abs(wavespec.data(allidx.loP,:))),1);
    columnspec_loP.mednorm(:,i) = nanmedian(wavespec.datan(allidx.loP,:),1);
    columnspec_loP.modz(:,i) = nanmedian(wavespec.dataz(allidx.loP,:),1);
    columnspec_loP.frac(:,i) = nanmedian(log10(wavespec.frac(allidx.loP,:)),1);
    columnspec_loP.osci(:,i) = nanmedian(wavespec.osci(allidx.loP,:),1);
    %columnspec_loP.osci(columnspec_loP.osci(:,i)<0,i) = 0;
    
    columnspec_hiP.db(:,i) = nanmedian(log10(abs(wavespec.data(allidx.hiP,:))),1);
    columnspec_hiP.mednorm(:,i) = nanmedian(wavespec.datan(allidx.hiP,:),1);
    columnspec_hiP.modz(:,i) = nanmedian(wavespec.dataz(allidx.hiP,:),1);
    columnspec_hiP.frac(:,i) = nanmedian(log10(wavespec.frac(allidx.hiP,:)),1);
    columnspec_hiP.osci(:,i) = nanmedian(wavespec.osci(allidx.hiP,:),1);
    %columnspec_hiP.osci(columnspec_hiP.osci(:,i)<0,i) = 0;
    
end

% Saving to struct
LayerSpectral.freqs = wavespec.freqs;
LayerSpectral.validfreqs = wavespec.validfreq;
LayerSpectral.columnspec_all = columnspec_all;
LayerSpectral.columnspec_NWh = columnspec_NWh;
LayerSpectral.columnspec_lNWh = columnspec_lNWh;
LayerSpectral.columnspec_Wh = columnspec_Wh;
LayerSpectral.columnspec_sWh = columnspec_sWh;
LayerSpectral.columnspec_lWh = columnspec_lWh;
LayerSpectral.columnspec_loWh = columnspec_loWh;
LayerSpectral.columnspec_hiWh = columnspec_hiWh;
LayerSpectral.columnspec_loP = columnspec_loP;
LayerSpectral.columnspec_hiP = columnspec_hiP;

%% FIGURE: NWh vs Wh Columnar Specs
cmax = max(max(columnspec_Wh.db(:,usechannels+1)-columnspec_NWh.db(:,usechannels+1)));

figure;
subplot(4,3,1);
imagesc(log10(wavespec.freqs),normdepth,...
    (columnspec_Wh.db(:,usechannels+1)-columnspec_NWh.db(:,usechannels+1))');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmax*-1 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
% xtickangle(45);
% xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('Pspec Wh-NWh');

cmin = min([min(min(columnspec_Wh.db(:,usechannels+1)))...
    min(min(columnspec_NWh.db(:,usechannels+1)))]);
cmax = max([max(max(columnspec_Wh.db(:,usechannels+1)))...
    max(max(columnspec_NWh.db(:,usechannels+1)))]);

subplot(4,3,2);
imagesc(log10(wavespec.freqs),normdepth,columnspec_NWh.db(:,usechannels+1)');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
% xtickangle(45);
% xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('PSpec NWh');

subplot(4,3,3);
imagesc(log10(wavespec.freqs),normdepth,columnspec_Wh.db(:,usechannels+1)');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
% xtickangle(45);
% xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('PSpec Wh');

% modZ
cmax = max(max(columnspec_Wh.modz(:,usechannels+1)-columnspec_NWh.modz(:,usechannels+1)));

subplot(4,3,4);
imagesc(log10(wavespec.freqs),normdepth,...
    (columnspec_Wh.modz(:,usechannels+1)-columnspec_NWh.modz(:,usechannels+1))');
axis tight
LogScale('x',10)
colormap('jet');
caxis([cmax*-1 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
% xtickangle(45);
% xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('modZ-Pspec Wh-NWh');

cmin = min([min(min(columnspec_Wh.modz(:,usechannels+1)))...
    min(min(columnspec_NWh.modz(:,usechannels+1)))]);
cmax = max([max(max(columnspec_Wh.modz(:,usechannels+1)))...
    max(max(columnspec_NWh.modz(:,usechannels+1)))]);

subplot(4,3,5);
imagesc(log10(wavespec.freqs),normdepth,columnspec_NWh.modz(:,usechannels+1)');
axis tight
LogScale('x',10)
colormap('jet');
caxis([cmax*-1 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
% xtickangle(45);
% xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('modZ-Pspec NWh');

subplot(4,3,6);
imagesc(log10(wavespec.freqs),normdepth,columnspec_Wh.modz(:,usechannels+1)');
axis tight
LogScale('x',10)
colormap('jet');
caxis([cmax*-1 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
% xtickangle(45);
% xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('modZ-Pspec Wh');

% frac
cmax = max(max(columnspec_Wh.frac(:,usechannels+1)-columnspec_NWh.frac(:,usechannels+1)));

subplot(4,3,7);
imagesc(log10(wavespec.validfreq),normdepth,...
    (columnspec_Wh.frac(:,usechannels+1)-columnspec_NWh.frac(:,usechannels+1))');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmax*-1 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
% xtickangle(45);
% xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('Fractal Wh-NWh');

cmin = min([min(min(columnspec_Wh.frac(:,usechannels+1)))...
    min(min(columnspec_NWh.frac(:,usechannels+1)))]);
cmax = max([max(max(columnspec_Wh.frac(:,usechannels+1)))...
    max(max(columnspec_NWh.frac(:,usechannels+1)))]);

subplot(4,3,8);
imagesc(log10(wavespec.validfreq),normdepth,columnspec_NWh.frac(:,usechannels+1)');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
% xtickangle(45);
% xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('Fractal NWh');

subplot(4,3,9);
imagesc(log10(wavespec.validfreq),normdepth,columnspec_Wh.frac(:,usechannels+1)');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
% xtickangle(45);
% xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('Fractal Wh');

% osci
cmax = max(max(columnspec_Wh.osci(:,usechannels+1)-columnspec_NWh.osci(:,usechannels+1)));

subplot(4,3,10);
imagesc(log10(wavespec.validfreq),normdepth,...
    (columnspec_Wh.osci(:,usechannels+1)-columnspec_NWh.osci(:,usechannels+1))');
axis tight
LogScale('x',10)
colormap('jet');
caxis([cmax*-1 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
% xtickangle(45);
% xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
xlabel('f (Hz)');
title('Osci Wh-NWh');

cmin = min([min(min(columnspec_Wh.osci(:,usechannels+1)))...
    min(min(columnspec_NWh.osci(:,usechannels+1)))]);
cmax = max([max(max(columnspec_Wh.osci(:,usechannels+1)))...
    max(max(columnspec_NWh.osci(:,usechannels+1)))]);

subplot(4,3,11);
imagesc(log10(wavespec.validfreq),normdepth,columnspec_NWh.osci(:,usechannels+1)');
axis tight
LogScale('x',10)
colormap('jet');
caxis([cmin cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
% xtickangle(45);
% xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
xlabel('f (Hz)');
title('Osci NWh');

subplot(4,3,12);
imagesc(log10(wavespec.validfreq),normdepth,columnspec_Wh.osci(:,usechannels+1)');
axis tight
LogScale('x',10)
colormap('jet');
caxis([cmin cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
% xtickangle(45);
% xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
xlabel('f (Hz)');
title('Osci Wh');

NiceSave('ColumnarSpecs_Wh_NWh',figfolder,baseName)

%% FIGURE 2: lo-hi Pupil
cmax = max(max(columnspec_hiP.db(:,usechannels+1)-columnspec_loP.db(:,usechannels+1)));

figure;
subplot(4,3,1);
imagesc(log10(wavespec.freqs),normdepth,...
    (columnspec_hiP.db(:,usechannels+1)-columnspec_loP.db(:,usechannels+1))');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmax*-1 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
% xtickangle(45);
% xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('Pspec hiP-loP');

cmin = min([min(min(columnspec_hiP.db(:,usechannels+1)))...
    min(min(columnspec_loP.db(:,usechannels+1)))]);
cmax = max([max(max(columnspec_hiP.db(:,usechannels+1)))...
    max(max(columnspec_loP.db(:,usechannels+1)))]);

subplot(4,3,2);
imagesc(log10(wavespec.freqs),normdepth,columnspec_loP.db(:,usechannels+1)');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
% xtickangle(45);
% xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('PSpec loP');

subplot(4,3,3);
imagesc(log10(wavespec.freqs),normdepth,columnspec_hiP.db(:,usechannels+1)');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
% xtickangle(45);
% xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('PSpec hiP');

% modZ
cmax = max(max(columnspec_hiP.modz(:,usechannels+1)-columnspec_loP.modz(:,usechannels+1)));

subplot(4,3,4);
imagesc(log10(wavespec.freqs),normdepth,...
    (columnspec_hiP.modz(:,usechannels+1)-columnspec_loP.modz(:,usechannels+1))');
axis tight
LogScale('x',10)
colormap('jet');
caxis([cmax*-1 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
% xtickangle(45);
% xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('modZ-Pspec hiP-loP');

cmin = min([min(min(columnspec_hiP.modz(:,usechannels+1)))...
    min(min(columnspec_loP.modz(:,usechannels+1)))]);
cmax = max([max(max(columnspec_hiP.modz(:,usechannels+1)))...
    max(max(columnspec_loP.modz(:,usechannels+1)))]);

subplot(4,3,5);
imagesc(log10(wavespec.freqs),normdepth,columnspec_loP.modz(:,usechannels+1)');
axis tight
LogScale('x',10)
colormap('jet');
caxis([cmax*-1 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
% xtickangle(45);
% xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('modZ-Pspec loP');

subplot(4,3,6);
imagesc(log10(wavespec.freqs),normdepth,columnspec_hiP.modz(:,usechannels+1)');
axis tight
LogScale('x',10)
colormap('jet');
caxis([cmax*-1 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
% xtickangle(45);
% xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('modZ-Pspec hiP');

% frac
cmax = max(max(columnspec_hiP.frac(:,usechannels+1)-columnspec_loP.frac(:,usechannels+1)));

subplot(4,3,7);
imagesc(log10(wavespec.validfreq),normdepth,...
    (columnspec_hiP.frac(:,usechannels+1)-columnspec_loP.frac(:,usechannels+1))');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmax*-1 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
% xtickangle(45);
% xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('Fractal hiP-loP');

cmin = min([min(min(columnspec_hiP.frac(:,usechannels+1)))...
    min(min(columnspec_loP.frac(:,usechannels+1)))]);
cmax = max([max(max(columnspec_hiP.frac(:,usechannels+1)))...
    max(max(columnspec_loP.frac(:,usechannels+1)))]);

subplot(4,3,8);
imagesc(log10(wavespec.validfreq),normdepth,columnspec_loP.frac(:,usechannels+1)');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
% xtickangle(45);
% xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('Fractal loP');

subplot(4,3,9);
imagesc(log10(wavespec.validfreq),normdepth,columnspec_hiP.frac(:,usechannels+1)');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
% xtickangle(45);
% xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('Fractal hiP');

% osci
cmax = max(max(columnspec_hiP.osci(:,usechannels+1)-columnspec_loP.osci(:,usechannels+1)));

subplot(4,3,10);
imagesc(log10(wavespec.validfreq),normdepth,...
    (columnspec_hiP.osci(:,usechannels+1)-columnspec_loP.osci(:,usechannels+1))');
axis tight
LogScale('x',10)
colormap('jet');
caxis([cmax*-1 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
% xtickangle(45);
% xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
xlabel('f (Hz)');
title('Osci hiP-loP');

cmin = min([min(min(columnspec_hiP.osci(:,usechannels+1)))...
    min(min(columnspec_loP.osci(:,usechannels+1)))]);
cmax = max([max(max(columnspec_hiP.osci(:,usechannels+1)))...
    max(max(columnspec_loP.osci(:,usechannels+1)))]);

subplot(4,3,11);
imagesc(log10(wavespec.validfreq),normdepth,columnspec_loP.osci(:,usechannels+1)');
axis tight
LogScale('x',10)
colormap('jet');
caxis([cmin cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
% xtickangle(45);
% xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
xlabel('f (Hz)');
title('Osci loP');

subplot(4,3,12);
imagesc(log10(wavespec.validfreq),normdepth,columnspec_hiP.osci(:,usechannels+1)');
axis tight
LogScale('x',10)
colormap('jet');
caxis([cmin cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
% xtickangle(45);
% xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
xlabel('f (Hz)');
title('Osci hiP');

NiceSave('ColumnarSpecs_hiP_loP',figfolder,baseName)

%% Laminar eventSpec by state
twin = [0.75 0.75].*wavespec.samplingRate;
laminarspec = NaN(size(wavespec.data,1),size(wavespec.data,2),6);
laminarfrac = NaN(size(wavespec.data,1),length(validFreqInds),6);
laminarosci = NaN(size(wavespec.data,1),length(validFreqInds),6);

% Layer 1
lidx = 1;
chidx = L1idx;
eventspec.Wh.spec = []; eventspec.Wh.frac = []; eventspec.Wh.pss = []; eventspec.Wh.osci = [];
eventspec.clearWh.spec = []; eventspec.clearWh.frac = []; eventspec.clearWh.pss = []; eventspec.clearWh.osci = [];
eventspec.sWh.spec = []; eventspec.sWh.frac = []; eventspec.sWh.pss = []; eventspec.sWh.osci = [];
eventspec.lWh.spec = []; eventspec.lWh.frac = []; eventspec.lWh.pss = []; eventspec.lWh.osci = [];
eventspec.loWh.spec = []; eventspec.loWh.frac = []; eventspec.loWh.pss = []; eventspec.loWh.osci = [];
eventspec.hiWh.spec = []; eventspec.hiWh.frac = []; eventspec.hiWh.pss = []; eventspec.hiWh.osci = [];
eventspec.Pup.spec = []; eventspec.Pup.frac = []; eventspec.Pup.pss = []; eventspec.Pup.osci = [];
for x = 1:length(chidx)
    x
    % Loading spectrograms
    load(fullfile(savefolder,[baseName,'.',num2str(chidx(x)),'.WaveSpec.lfp.mat']));
    %load(fullfile(savefolder,[baseName,'.',num2str(chidx(x)),'.WaveSpec2.lfp.mat']));
    wavespec.dataz = NormToInt(log10(abs(wavespec.data)),'modZ');
    
    % WaveIRASA wavelet spec
    [wavespec.frac,wavespec.osci,wavespec.validfreq] = WaveIRASA(wavespec);
    %wavespec.osci(wavespec.osci<0) = 0;
    wavespec = WaveIRASA_plawfit(wavespec, [2.5 100]);

    % Averaging...
    tempSpec = cat(3,laminarspec(:,:,lidx),log10(abs(wavespec.data)));
    laminarspec(:,:,lidx) = nansum(tempSpec,3);
    
    tempFrac = cat(3,laminarfrac(:,:,lidx),log10(wavespec.frac));
    laminarfrac(:,:,lidx) = nansum(tempFrac,3);
    
    tempOsci = cat(3,laminarosci(:,:,lidx),wavespec.osci);
    laminarosci(:,:,lidx) = nansum(tempOsci,3);
   
    % eventSpec
    teventSpec = eventSpec(wavespec,eventsidx.Wh(:,1));
    eventspec.Wh = bz_CollapseStruct([eventspec.Wh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.clearWh(:,1));
    eventspec.clearWh = bz_CollapseStruct([eventspec.clearWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.sWh(:,1));
    eventspec.sWh = bz_CollapseStruct([eventspec.sWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.lWh(:,1));
    eventspec.lWh = bz_CollapseStruct([eventspec.lWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.loWh(:,1));
    eventspec.loWh = bz_CollapseStruct([eventspec.loWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.hiWh(:,1));
    eventspec.hiWh = bz_CollapseStruct([eventspec.hiWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.Pup(:,1));
    eventspec.Pup = bz_CollapseStruct([eventspec.Pup teventSpec],3,'justcat',true);
end
eventspec.Wh.spec = squeeze(nanmean(eventspec.Wh.spec,3)); 
eventspec.Wh.frac = squeeze(nanmean(eventspec.Wh.frac,3));  
eventspec.Wh.pss = squeeze(nanmean(eventspec.Wh.pss,3));  
eventspec.Wh.osci = squeeze(nanmean(eventspec.Wh.osci,3)); 
eventspec.clearWh.spec = squeeze(nanmean(eventspec.clearWh.spec,3)); 
eventspec.clearWh.frac = squeeze(nanmean(eventspec.clearWh.frac,3));
eventspec.clearWh.pss = squeeze(nanmean(eventspec.clearWh.pss,3));  
eventspec.clearWh.osci = squeeze(nanmean(eventspec.clearWh.osci,3));
eventspec.sWh.spec = squeeze(nanmean(eventspec.sWh.spec,3)); 
eventspec.sWh.frac = squeeze(nanmean(eventspec.sWh.frac,3));  
eventspec.sWh.pss = squeeze(nanmean(eventspec.sWh.pss,3));  
eventspec.sWh.osci = squeeze(nanmean(eventspec.sWh.osci,3));
eventspec.lWh.spec = squeeze(nanmean(eventspec.lWh.spec,3)); 
eventspec.lWh.frac = squeeze(nanmean(eventspec.lWh.frac,3));  
eventspec.lWh.pss = squeeze(nanmean(eventspec.lWh.pss,3));  
eventspec.lWh.osci = squeeze(nanmean(eventspec.lWh.osci,3));
eventspec.loWh.spec = squeeze(nanmean(eventspec.loWh.spec,3)); 
eventspec.loWh.frac = squeeze(nanmean(eventspec.loWh.frac,3));  
eventspec.loWh.pss = squeeze(nanmean(eventspec.loWh.pss,3));  
eventspec.loWh.osci = squeeze(nanmean(eventspec.loWh.osci,3));
eventspec.hiWh.spec = squeeze(nanmean(eventspec.hiWh.spec,3)); 
eventspec.hiWh.frac = squeeze(nanmean(eventspec.hiWh.frac,3));
eventspec.hiWh.pss = squeeze(nanmean(eventspec.hiWh.pss,3));  
eventspec.hiWh.osci = squeeze(nanmean(eventspec.hiWh.osci,3));
eventspec.Pup.spec = squeeze(nanmean(eventspec.Pup.spec,3)); 
eventspec.Pup.frac = squeeze(nanmean(eventspec.Pup.frac,3)); 
eventspec.Pup.pss = squeeze(nanmean(eventspec.Pup.pss,3)); 
eventspec.Pup.osci = squeeze(nanmean(eventspec.Pup.osci,3)); 
laminarspec(:,:,lidx) = laminarspec(:,:,lidx)./length(chidx);
laminarfrac(:,:,lidx) = laminarfrac(:,:,lidx)./length(chidx);
laminarosci(:,:,lidx) = laminarosci(:,:,lidx)./length(chidx);
L1eventSpec = eventspec;

% Layer 2/3
lidx = 2;
chidx = L23idx;
eventspec.Wh.spec = []; eventspec.Wh.frac = []; eventspec.Wh.pss = []; eventspec.Wh.osci = [];
eventspec.clearWh.spec = []; eventspec.clearWh.frac = []; eventspec.clearWh.pss = []; eventspec.clearWh.osci = [];
eventspec.sWh.spec = []; eventspec.sWh.frac = []; eventspec.sWh.pss = []; eventspec.sWh.osci = [];
eventspec.lWh.spec = []; eventspec.lWh.frac = []; eventspec.lWh.pss = []; eventspec.lWh.osci = [];
eventspec.loWh.spec = []; eventspec.loWh.frac = []; eventspec.loWh.pss = []; eventspec.loWh.osci = [];
eventspec.hiWh.spec = []; eventspec.hiWh.frac = []; eventspec.hiWh.pss = []; eventspec.hiWh.osci = [];
eventspec.Pup.spec = []; eventspec.Pup.frac = []; eventspec.Pup.pss = []; eventspec.Pup.osci = [];
for x = 1:length(chidx)
    x
    % Loading spectrograms
    load(fullfile(savefolder,[baseName,'.',num2str(chidx(x)),'.WaveSpec.lfp.mat']));
    %load(fullfile(savefolder,[baseName,'.',num2str(chidx(x)),'.WaveSpec2.lfp.mat']));
    wavespec.dataz = NormToInt(log10(abs(wavespec.data)),'modZ');
    
    % WaveIRASA wavelet spec
    [wavespec.frac,wavespec.osci,wavespec.validfreq] = WaveIRASA(wavespec);
    %wavespec.osci(wavespec.osci<0) = 0;
    wavespec = WaveIRASA_plawfit(wavespec, [2.5 100]);

    % Averaging...
    tempSpec = cat(3,laminarspec(:,:,lidx),log10(abs(wavespec.data)));
    laminarspec(:,:,lidx) = nansum(tempSpec,3);
    
    tempFrac = cat(3,laminarfrac(:,:,lidx),log10(wavespec.frac));
    laminarfrac(:,:,lidx) = nansum(tempFrac,3);
    
    tempOsci = cat(3,laminarosci(:,:,lidx),wavespec.osci);
    laminarosci(:,:,lidx) = nansum(tempOsci,3);
   
    % eventSpec
    teventSpec = eventSpec(wavespec,eventsidx.Wh(:,1));
    eventspec.Wh = bz_CollapseStruct([eventspec.Wh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.clearWh(:,1));
    eventspec.clearWh = bz_CollapseStruct([eventspec.clearWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.sWh(:,1));
    eventspec.sWh = bz_CollapseStruct([eventspec.sWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.lWh(:,1));
    eventspec.lWh = bz_CollapseStruct([eventspec.lWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.loWh(:,1));
    eventspec.loWh = bz_CollapseStruct([eventspec.loWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.hiWh(:,1));
    eventspec.hiWh = bz_CollapseStruct([eventspec.hiWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.Pup(:,1));
    eventspec.Pup = bz_CollapseStruct([eventspec.Pup teventSpec],3,'justcat',true);
end
eventspec.Wh.spec = squeeze(nanmean(eventspec.Wh.spec,3)); 
eventspec.Wh.frac = squeeze(nanmean(eventspec.Wh.frac,3));  
eventspec.Wh.pss = squeeze(nanmean(eventspec.Wh.pss,3));  
eventspec.Wh.osci = squeeze(nanmean(eventspec.Wh.osci,3)); 
eventspec.clearWh.spec = squeeze(nanmean(eventspec.clearWh.spec,3)); 
eventspec.clearWh.frac = squeeze(nanmean(eventspec.clearWh.frac,3));
eventspec.clearWh.pss = squeeze(nanmean(eventspec.clearWh.pss,3));  
eventspec.clearWh.osci = squeeze(nanmean(eventspec.clearWh.osci,3));
eventspec.sWh.spec = squeeze(nanmean(eventspec.sWh.spec,3)); 
eventspec.sWh.frac = squeeze(nanmean(eventspec.sWh.frac,3));  
eventspec.sWh.pss = squeeze(nanmean(eventspec.sWh.pss,3));  
eventspec.sWh.osci = squeeze(nanmean(eventspec.sWh.osci,3));
eventspec.lWh.spec = squeeze(nanmean(eventspec.lWh.spec,3)); 
eventspec.lWh.frac = squeeze(nanmean(eventspec.lWh.frac,3));  
eventspec.lWh.pss = squeeze(nanmean(eventspec.lWh.pss,3));  
eventspec.lWh.osci = squeeze(nanmean(eventspec.lWh.osci,3));
eventspec.loWh.spec = squeeze(nanmean(eventspec.loWh.spec,3)); 
eventspec.loWh.frac = squeeze(nanmean(eventspec.loWh.frac,3));  
eventspec.loWh.pss = squeeze(nanmean(eventspec.loWh.pss,3));  
eventspec.loWh.osci = squeeze(nanmean(eventspec.loWh.osci,3));
eventspec.hiWh.spec = squeeze(nanmean(eventspec.hiWh.spec,3)); 
eventspec.hiWh.frac = squeeze(nanmean(eventspec.hiWh.frac,3));
eventspec.hiWh.pss = squeeze(nanmean(eventspec.hiWh.pss,3));  
eventspec.hiWh.osci = squeeze(nanmean(eventspec.hiWh.osci,3));
eventspec.Pup.spec = squeeze(nanmean(eventspec.Pup.spec,3)); 
eventspec.Pup.frac = squeeze(nanmean(eventspec.Pup.frac,3)); 
eventspec.Pup.pss = squeeze(nanmean(eventspec.Pup.pss,3)); 
eventspec.Pup.osci = squeeze(nanmean(eventspec.Pup.osci,3)); 
laminarspec(:,:,lidx) = laminarspec(:,:,lidx)./length(chidx);
laminarfrac(:,:,lidx) = laminarfrac(:,:,lidx)./length(chidx);
laminarosci(:,:,lidx) = laminarosci(:,:,lidx)./length(chidx);
L23eventSpec = eventspec;

% Layer 4
lidx = 3;
chidx = L4idx;
eventspec.Wh.spec = []; eventspec.Wh.frac = []; eventspec.Wh.pss = []; eventspec.Wh.osci = [];
eventspec.clearWh.spec = []; eventspec.clearWh.frac = []; eventspec.clearWh.pss = []; eventspec.clearWh.osci = [];
eventspec.sWh.spec = []; eventspec.sWh.frac = []; eventspec.sWh.pss = []; eventspec.sWh.osci = [];
eventspec.lWh.spec = []; eventspec.lWh.frac = []; eventspec.lWh.pss = []; eventspec.lWh.osci = [];
eventspec.loWh.spec = []; eventspec.loWh.frac = []; eventspec.loWh.pss = []; eventspec.loWh.osci = [];
eventspec.hiWh.spec = []; eventspec.hiWh.frac = []; eventspec.hiWh.pss = []; eventspec.hiWh.osci = [];
eventspec.Pup.spec = []; eventspec.Pup.frac = []; eventspec.Pup.pss = []; eventspec.Pup.osci = [];
for x = 1:length(chidx)
    x
    % Loading spectrograms
    load(fullfile(savefolder,[baseName,'.',num2str(chidx(x)),'.WaveSpec.lfp.mat']));
    %load(fullfile(savefolder,[baseName,'.',num2str(chidx(x)),'.WaveSpec2.lfp.mat']));
    wavespec.dataz = NormToInt(log10(abs(wavespec.data)),'modZ');
    
    % WaveIRASA wavelet spec
    [wavespec.frac,wavespec.osci,wavespec.validfreq] = WaveIRASA(wavespec);
    %wavespec.osci(wavespec.osci<0) = 0;
    wavespec = WaveIRASA_plawfit(wavespec, [2.5 100]);

    % Averaging...
    tempSpec = cat(3,laminarspec(:,:,lidx),log10(abs(wavespec.data)));
    laminarspec(:,:,lidx) = nansum(tempSpec,3);
    
    tempFrac = cat(3,laminarfrac(:,:,lidx),log10(wavespec.frac));
    laminarfrac(:,:,lidx) = nansum(tempFrac,3);
    
    tempOsci = cat(3,laminarosci(:,:,lidx),wavespec.osci);
    laminarosci(:,:,lidx) = nansum(tempOsci,3);
   
    % eventSpec
    teventSpec = eventSpec(wavespec,eventsidx.Wh(:,1));
    eventspec.Wh = bz_CollapseStruct([eventspec.Wh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.clearWh(:,1));
    eventspec.clearWh = bz_CollapseStruct([eventspec.clearWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.sWh(:,1));
    eventspec.sWh = bz_CollapseStruct([eventspec.sWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.lWh(:,1));
    eventspec.lWh = bz_CollapseStruct([eventspec.lWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.loWh(:,1));
    eventspec.loWh = bz_CollapseStruct([eventspec.loWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.hiWh(:,1));
    eventspec.hiWh = bz_CollapseStruct([eventspec.hiWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.Pup(:,1));
    eventspec.Pup = bz_CollapseStruct([eventspec.Pup teventSpec],3,'justcat',true);
end
eventspec.Wh.spec = squeeze(nanmean(eventspec.Wh.spec,3)); 
eventspec.Wh.frac = squeeze(nanmean(eventspec.Wh.frac,3));  
eventspec.Wh.pss = squeeze(nanmean(eventspec.Wh.pss,3));  
eventspec.Wh.osci = squeeze(nanmean(eventspec.Wh.osci,3)); 
eventspec.clearWh.spec = squeeze(nanmean(eventspec.clearWh.spec,3)); 
eventspec.clearWh.frac = squeeze(nanmean(eventspec.clearWh.frac,3));
eventspec.clearWh.pss = squeeze(nanmean(eventspec.clearWh.pss,3));  
eventspec.clearWh.osci = squeeze(nanmean(eventspec.clearWh.osci,3));
eventspec.sWh.spec = squeeze(nanmean(eventspec.sWh.spec,3)); 
eventspec.sWh.frac = squeeze(nanmean(eventspec.sWh.frac,3));  
eventspec.sWh.pss = squeeze(nanmean(eventspec.sWh.pss,3));  
eventspec.sWh.osci = squeeze(nanmean(eventspec.sWh.osci,3));
eventspec.lWh.spec = squeeze(nanmean(eventspec.lWh.spec,3)); 
eventspec.lWh.frac = squeeze(nanmean(eventspec.lWh.frac,3));  
eventspec.lWh.pss = squeeze(nanmean(eventspec.lWh.pss,3));  
eventspec.lWh.osci = squeeze(nanmean(eventspec.lWh.osci,3));
eventspec.loWh.spec = squeeze(nanmean(eventspec.loWh.spec,3)); 
eventspec.loWh.frac = squeeze(nanmean(eventspec.loWh.frac,3));  
eventspec.loWh.pss = squeeze(nanmean(eventspec.loWh.pss,3));  
eventspec.loWh.osci = squeeze(nanmean(eventspec.loWh.osci,3));
eventspec.hiWh.spec = squeeze(nanmean(eventspec.hiWh.spec,3)); 
eventspec.hiWh.frac = squeeze(nanmean(eventspec.hiWh.frac,3));
eventspec.hiWh.pss = squeeze(nanmean(eventspec.hiWh.pss,3));  
eventspec.hiWh.osci = squeeze(nanmean(eventspec.hiWh.osci,3));
eventspec.Pup.spec = squeeze(nanmean(eventspec.Pup.spec,3)); 
eventspec.Pup.frac = squeeze(nanmean(eventspec.Pup.frac,3)); 
eventspec.Pup.pss = squeeze(nanmean(eventspec.Pup.pss,3)); 
eventspec.Pup.osci = squeeze(nanmean(eventspec.Pup.osci,3)); 
laminarspec(:,:,lidx) = laminarspec(:,:,lidx)./length(chidx);
laminarfrac(:,:,lidx) = laminarfrac(:,:,lidx)./length(chidx);
laminarosci(:,:,lidx) = laminarosci(:,:,lidx)./length(chidx);
L4eventSpec = eventspec;

% Layer 5a
lidx = 4;
chidx = L5aidx;
eventspec.Wh.spec = []; eventspec.Wh.frac = []; eventspec.Wh.pss = []; eventspec.Wh.osci = [];
eventspec.clearWh.spec = []; eventspec.clearWh.frac = []; eventspec.clearWh.pss = []; eventspec.clearWh.osci = [];
eventspec.sWh.spec = []; eventspec.sWh.frac = []; eventspec.sWh.pss = []; eventspec.sWh.osci = [];
eventspec.lWh.spec = []; eventspec.lWh.frac = []; eventspec.lWh.pss = []; eventspec.lWh.osci = [];
eventspec.loWh.spec = []; eventspec.loWh.frac = []; eventspec.loWh.pss = []; eventspec.loWh.osci = [];
eventspec.hiWh.spec = []; eventspec.hiWh.frac = []; eventspec.hiWh.pss = []; eventspec.hiWh.osci = [];
eventspec.Pup.spec = []; eventspec.Pup.frac = []; eventspec.Pup.pss = []; eventspec.Pup.osci = [];
for x = 1:length(chidx)
    x
    % Loading spectrograms
    load(fullfile(savefolder,[baseName,'.',num2str(chidx(x)),'.WaveSpec.lfp.mat']));
    %load(fullfile(savefolder,[baseName,'.',num2str(chidx(x)),'.WaveSpec2.lfp.mat']));
    wavespec.dataz = NormToInt(log10(abs(wavespec.data)),'modZ');
    
    % WaveIRASA wavelet spec
    [wavespec.frac,wavespec.osci,wavespec.validfreq] = WaveIRASA(wavespec);
    %wavespec.osci(wavespec.osci<0) = 0;
    wavespec = WaveIRASA_plawfit(wavespec, [2.5 100]);

    % Averaging...
    tempSpec = cat(3,laminarspec(:,:,lidx),log10(abs(wavespec.data)));
    laminarspec(:,:,lidx) = nansum(tempSpec,3);
    
    tempFrac = cat(3,laminarfrac(:,:,lidx),log10(wavespec.frac));
    laminarfrac(:,:,lidx) = nansum(tempFrac,3);
    
    tempOsci = cat(3,laminarosci(:,:,lidx),wavespec.osci);
    laminarosci(:,:,lidx) = nansum(tempOsci,3);
   
    % eventSpec
    teventSpec = eventSpec(wavespec,eventsidx.Wh(:,1));
    eventspec.Wh = bz_CollapseStruct([eventspec.Wh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.clearWh(:,1));
    eventspec.clearWh = bz_CollapseStruct([eventspec.clearWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.sWh(:,1));
    eventspec.sWh = bz_CollapseStruct([eventspec.sWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.lWh(:,1));
    eventspec.lWh = bz_CollapseStruct([eventspec.lWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.loWh(:,1));
    eventspec.loWh = bz_CollapseStruct([eventspec.loWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.hiWh(:,1));
    eventspec.hiWh = bz_CollapseStruct([eventspec.hiWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.Pup(:,1));
    eventspec.Pup = bz_CollapseStruct([eventspec.Pup teventSpec],3,'justcat',true);
end
eventspec.Wh.spec = squeeze(nanmean(eventspec.Wh.spec,3)); 
eventspec.Wh.frac = squeeze(nanmean(eventspec.Wh.frac,3));  
eventspec.Wh.pss = squeeze(nanmean(eventspec.Wh.pss,3));  
eventspec.Wh.osci = squeeze(nanmean(eventspec.Wh.osci,3)); 
eventspec.clearWh.spec = squeeze(nanmean(eventspec.clearWh.spec,3)); 
eventspec.clearWh.frac = squeeze(nanmean(eventspec.clearWh.frac,3));
eventspec.clearWh.pss = squeeze(nanmean(eventspec.clearWh.pss,3));  
eventspec.clearWh.osci = squeeze(nanmean(eventspec.clearWh.osci,3));
eventspec.sWh.spec = squeeze(nanmean(eventspec.sWh.spec,3)); 
eventspec.sWh.frac = squeeze(nanmean(eventspec.sWh.frac,3));  
eventspec.sWh.pss = squeeze(nanmean(eventspec.sWh.pss,3));  
eventspec.sWh.osci = squeeze(nanmean(eventspec.sWh.osci,3));
eventspec.lWh.spec = squeeze(nanmean(eventspec.lWh.spec,3)); 
eventspec.lWh.frac = squeeze(nanmean(eventspec.lWh.frac,3));  
eventspec.lWh.pss = squeeze(nanmean(eventspec.lWh.pss,3));  
eventspec.lWh.osci = squeeze(nanmean(eventspec.lWh.osci,3));
eventspec.loWh.spec = squeeze(nanmean(eventspec.loWh.spec,3)); 
eventspec.loWh.frac = squeeze(nanmean(eventspec.loWh.frac,3));  
eventspec.loWh.pss = squeeze(nanmean(eventspec.loWh.pss,3));  
eventspec.loWh.osci = squeeze(nanmean(eventspec.loWh.osci,3));
eventspec.hiWh.spec = squeeze(nanmean(eventspec.hiWh.spec,3)); 
eventspec.hiWh.frac = squeeze(nanmean(eventspec.hiWh.frac,3));
eventspec.hiWh.pss = squeeze(nanmean(eventspec.hiWh.pss,3));  
eventspec.hiWh.osci = squeeze(nanmean(eventspec.hiWh.osci,3));
eventspec.Pup.spec = squeeze(nanmean(eventspec.Pup.spec,3)); 
eventspec.Pup.frac = squeeze(nanmean(eventspec.Pup.frac,3)); 
eventspec.Pup.pss = squeeze(nanmean(eventspec.Pup.pss,3)); 
eventspec.Pup.osci = squeeze(nanmean(eventspec.Pup.osci,3)); 
laminarspec(:,:,lidx) = laminarspec(:,:,lidx)./length(chidx);
laminarfrac(:,:,lidx) = laminarfrac(:,:,lidx)./length(chidx);
laminarosci(:,:,lidx) = laminarosci(:,:,lidx)./length(chidx);
L5aeventSpec = eventspec;

% Layer 5b/6
lidx = 5;
chidx = L56idx;
eventspec.Wh.spec = []; eventspec.Wh.frac = []; eventspec.Wh.pss = []; eventspec.Wh.osci = [];
eventspec.clearWh.spec = []; eventspec.clearWh.frac = []; eventspec.clearWh.pss = []; eventspec.clearWh.osci = [];
eventspec.sWh.spec = []; eventspec.sWh.frac = []; eventspec.sWh.pss = []; eventspec.sWh.osci = [];
eventspec.lWh.spec = []; eventspec.lWh.frac = []; eventspec.lWh.pss = []; eventspec.lWh.osci = [];
eventspec.loWh.spec = []; eventspec.loWh.frac = []; eventspec.loWh.pss = []; eventspec.loWh.osci = [];
eventspec.hiWh.spec = []; eventspec.hiWh.frac = []; eventspec.hiWh.pss = []; eventspec.hiWh.osci = [];
eventspec.Pup.spec = []; eventspec.Pup.frac = []; eventspec.Pup.pss = []; eventspec.Pup.osci = [];
for x = 1:length(chidx)
    x
    % Loading spectrograms
    load(fullfile(savefolder,[baseName,'.',num2str(chidx(x)),'.WaveSpec.lfp.mat']));
    %load(fullfile(savefolder,[baseName,'.',num2str(chidx(x)),'.WaveSpec2.lfp.mat']));
    wavespec.dataz = NormToInt(log10(abs(wavespec.data)),'modZ');
    
    % WaveIRASA wavelet spec
    [wavespec.frac,wavespec.osci,wavespec.validfreq] = WaveIRASA(wavespec);
    %wavespec.osci(wavespec.osci<0) = 0;
    wavespec = WaveIRASA_plawfit(wavespec, [2.5 100]);

    % Averaging...
    tempSpec = cat(3,laminarspec(:,:,lidx),log10(abs(wavespec.data)));
    laminarspec(:,:,lidx) = nansum(tempSpec,3);
    
    tempFrac = cat(3,laminarfrac(:,:,lidx),log10(wavespec.frac));
    laminarfrac(:,:,lidx) = nansum(tempFrac,3);
    
    tempOsci = cat(3,laminarosci(:,:,lidx),wavespec.osci);
    laminarosci(:,:,lidx) = nansum(tempOsci,3);
   
    % eventSpec
    teventSpec = eventSpec(wavespec,eventsidx.Wh(:,1));
    eventspec.Wh = bz_CollapseStruct([eventspec.Wh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.clearWh(:,1));
    eventspec.clearWh = bz_CollapseStruct([eventspec.clearWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.sWh(:,1));
    eventspec.sWh = bz_CollapseStruct([eventspec.sWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.lWh(:,1));
    eventspec.lWh = bz_CollapseStruct([eventspec.lWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.loWh(:,1));
    eventspec.loWh = bz_CollapseStruct([eventspec.loWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.hiWh(:,1));
    eventspec.hiWh = bz_CollapseStruct([eventspec.hiWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.Pup(:,1));
    eventspec.Pup = bz_CollapseStruct([eventspec.Pup teventSpec],3,'justcat',true);
end
eventspec.Wh.spec = squeeze(nanmean(eventspec.Wh.spec,3)); 
eventspec.Wh.frac = squeeze(nanmean(eventspec.Wh.frac,3));  
eventspec.Wh.pss = squeeze(nanmean(eventspec.Wh.pss,3));  
eventspec.Wh.osci = squeeze(nanmean(eventspec.Wh.osci,3)); 
eventspec.clearWh.spec = squeeze(nanmean(eventspec.clearWh.spec,3)); 
eventspec.clearWh.frac = squeeze(nanmean(eventspec.clearWh.frac,3));
eventspec.clearWh.pss = squeeze(nanmean(eventspec.clearWh.pss,3));  
eventspec.clearWh.osci = squeeze(nanmean(eventspec.clearWh.osci,3));
eventspec.sWh.spec = squeeze(nanmean(eventspec.sWh.spec,3)); 
eventspec.sWh.frac = squeeze(nanmean(eventspec.sWh.frac,3));  
eventspec.sWh.pss = squeeze(nanmean(eventspec.sWh.pss,3));  
eventspec.sWh.osci = squeeze(nanmean(eventspec.sWh.osci,3));
eventspec.lWh.spec = squeeze(nanmean(eventspec.lWh.spec,3)); 
eventspec.lWh.frac = squeeze(nanmean(eventspec.lWh.frac,3));  
eventspec.lWh.pss = squeeze(nanmean(eventspec.lWh.pss,3));  
eventspec.lWh.osci = squeeze(nanmean(eventspec.lWh.osci,3));
eventspec.loWh.spec = squeeze(nanmean(eventspec.loWh.spec,3)); 
eventspec.loWh.frac = squeeze(nanmean(eventspec.loWh.frac,3));  
eventspec.loWh.pss = squeeze(nanmean(eventspec.loWh.pss,3));  
eventspec.loWh.osci = squeeze(nanmean(eventspec.loWh.osci,3));
eventspec.hiWh.spec = squeeze(nanmean(eventspec.hiWh.spec,3)); 
eventspec.hiWh.frac = squeeze(nanmean(eventspec.hiWh.frac,3));
eventspec.hiWh.pss = squeeze(nanmean(eventspec.hiWh.pss,3));  
eventspec.hiWh.osci = squeeze(nanmean(eventspec.hiWh.osci,3));
eventspec.Pup.spec = squeeze(nanmean(eventspec.Pup.spec,3)); 
eventspec.Pup.frac = squeeze(nanmean(eventspec.Pup.frac,3)); 
eventspec.Pup.pss = squeeze(nanmean(eventspec.Pup.pss,3)); 
eventspec.Pup.osci = squeeze(nanmean(eventspec.Pup.osci,3)); 
laminarspec(:,:,lidx) = laminarspec(:,:,lidx)./length(chidx);
laminarfrac(:,:,lidx) = laminarfrac(:,:,lidx)./length(chidx);
laminarosci(:,:,lidx) = laminarosci(:,:,lidx)./length(chidx);
L56eventSpec = eventspec;

% Layer 6
lidx = 6;
chidx = L6idx;
eventspec.Wh.spec = []; eventspec.Wh.frac = []; eventspec.Wh.pss = []; eventspec.Wh.osci = [];
eventspec.clearWh.spec = []; eventspec.clearWh.frac = []; eventspec.clearWh.pss = []; eventspec.clearWh.osci = [];
eventspec.sWh.spec = []; eventspec.sWh.frac = []; eventspec.sWh.pss = []; eventspec.sWh.osci = [];
eventspec.lWh.spec = []; eventspec.lWh.frac = []; eventspec.lWh.pss = []; eventspec.lWh.osci = [];
eventspec.loWh.spec = []; eventspec.loWh.frac = []; eventspec.loWh.pss = []; eventspec.loWh.osci = [];
eventspec.hiWh.spec = []; eventspec.hiWh.frac = []; eventspec.hiWh.pss = []; eventspec.hiWh.osci = [];
eventspec.Pup.spec = []; eventspec.Pup.frac = []; eventspec.Pup.pss = []; eventspec.Pup.osci = [];
for x = 1:length(chidx)
    x
    % Loading spectrograms
    load(fullfile(savefolder,[baseName,'.',num2str(chidx(x)),'.WaveSpec.lfp.mat']));
    %load(fullfile(savefolder,[baseName,'.',num2str(chidx(x)),'.WaveSpec2.lfp.mat']));
    wavespec.dataz = NormToInt(log10(abs(wavespec.data)),'modZ');
    
    % WaveIRASA wavelet spec
    [wavespec.frac,wavespec.osci,wavespec.validfreq] = WaveIRASA(wavespec);
    %wavespec.osci(wavespec.osci<0) = 0;
    wavespec = WaveIRASA_plawfit(wavespec, [2.5 100]);

    % Averaging...
    tempSpec = cat(3,laminarspec(:,:,lidx),log10(abs(wavespec.data)));
    laminarspec(:,:,lidx) = nansum(tempSpec,3);
    
    tempFrac = cat(3,laminarfrac(:,:,lidx),log10(wavespec.frac));
    laminarfrac(:,:,lidx) = nansum(tempFrac,3);
    
    tempOsci = cat(3,laminarosci(:,:,lidx),wavespec.osci);
    laminarosci(:,:,lidx) = nansum(tempOsci,3);
   
    % eventSpec
    teventSpec = eventSpec(wavespec,eventsidx.Wh(:,1));
    eventspec.Wh = bz_CollapseStruct([eventspec.Wh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.clearWh(:,1));
    eventspec.clearWh = bz_CollapseStruct([eventspec.clearWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.sWh(:,1));
    eventspec.sWh = bz_CollapseStruct([eventspec.sWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.lWh(:,1));
    eventspec.lWh = bz_CollapseStruct([eventspec.lWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.loWh(:,1));
    eventspec.loWh = bz_CollapseStruct([eventspec.loWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.hiWh(:,1));
    eventspec.hiWh = bz_CollapseStruct([eventspec.hiWh teventSpec],3,'justcat',true);
    
    teventSpec = eventSpec(wavespec,eventsidx.Pup(:,1));
    eventspec.Pup = bz_CollapseStruct([eventspec.Pup teventSpec],3,'justcat',true);
end
eventspec.Wh.spec = squeeze(nanmean(eventspec.Wh.spec,3)); 
eventspec.Wh.frac = squeeze(nanmean(eventspec.Wh.frac,3));  
eventspec.Wh.pss = squeeze(nanmean(eventspec.Wh.pss,3));  
eventspec.Wh.osci = squeeze(nanmean(eventspec.Wh.osci,3)); 
eventspec.clearWh.spec = squeeze(nanmean(eventspec.clearWh.spec,3)); 
eventspec.clearWh.frac = squeeze(nanmean(eventspec.clearWh.frac,3));
eventspec.clearWh.pss = squeeze(nanmean(eventspec.clearWh.pss,3));  
eventspec.clearWh.osci = squeeze(nanmean(eventspec.clearWh.osci,3));
eventspec.sWh.spec = squeeze(nanmean(eventspec.sWh.spec,3)); 
eventspec.sWh.frac = squeeze(nanmean(eventspec.sWh.frac,3));  
eventspec.sWh.pss = squeeze(nanmean(eventspec.sWh.pss,3));  
eventspec.sWh.osci = squeeze(nanmean(eventspec.sWh.osci,3));
eventspec.lWh.spec = squeeze(nanmean(eventspec.lWh.spec,3)); 
eventspec.lWh.frac = squeeze(nanmean(eventspec.lWh.frac,3));  
eventspec.lWh.pss = squeeze(nanmean(eventspec.lWh.pss,3));  
eventspec.lWh.osci = squeeze(nanmean(eventspec.lWh.osci,3));
eventspec.loWh.spec = squeeze(nanmean(eventspec.loWh.spec,3)); 
eventspec.loWh.frac = squeeze(nanmean(eventspec.loWh.frac,3));  
eventspec.loWh.pss = squeeze(nanmean(eventspec.loWh.pss,3));  
eventspec.loWh.osci = squeeze(nanmean(eventspec.loWh.osci,3));
eventspec.hiWh.spec = squeeze(nanmean(eventspec.hiWh.spec,3)); 
eventspec.hiWh.frac = squeeze(nanmean(eventspec.hiWh.frac,3));
eventspec.hiWh.pss = squeeze(nanmean(eventspec.hiWh.pss,3));  
eventspec.hiWh.osci = squeeze(nanmean(eventspec.hiWh.osci,3));
eventspec.Pup.spec = squeeze(nanmean(eventspec.Pup.spec,3)); 
eventspec.Pup.frac = squeeze(nanmean(eventspec.Pup.frac,3)); 
eventspec.Pup.pss = squeeze(nanmean(eventspec.Pup.pss,3)); 
eventspec.Pup.osci = squeeze(nanmean(eventspec.Pup.osci,3)); 
laminarspec(:,:,lidx) = laminarspec(:,:,lidx)./length(chidx);
laminarfrac(:,:,lidx) = laminarfrac(:,:,lidx)./length(chidx);
laminarosci(:,:,lidx) = laminarosci(:,:,lidx)./length(chidx);
L6eventSpec = eventspec;

% Saving to struct
LayerSpectral.laminarspec = laminarspec;
LayerSpectral.laminarfrac = laminarfrac;
LayerSpectral.laminarosci = laminarosci;
LayerSpectral.L1eventSpec = L1eventSpec; 
LayerSpectral.L23eventSpec = L23eventSpec; 
LayerSpectral.L4eventSpec = L4eventSpec; 
LayerSpectral.L5aeventSpec = L5aeventSpec; 
LayerSpectral.L56eventSpec = L56eventSpec; 
LayerSpectral.L6eventSpec = L6eventSpec; 

%% FIGURE 3: Whisking-aligned Spec/Frac/Osci
twin = round([0.75 0.75].*wavespec.samplingRate);
taxis = (-(twin(1)/wavespec.samplingRate):(1/wavespec.samplingRate):(twin(2)/wavespec.samplingRate))*1e3;

cmin = min([min(min(L1eventSpec.Wh.spec)) min(min(L23eventSpec.Wh.spec))...
    min(min(L4eventSpec.Wh.spec)) min(min(L5aeventSpec.Wh.spec))...
    min(min(L56eventSpec.Wh.spec)) min(min(L6eventSpec.Wh.spec))]);
cmax = max([max(max(L1eventSpec.Wh.spec)) max(max(L23eventSpec.Wh.spec))...
    max(max(L4eventSpec.Wh.spec)) max(max(L5aeventSpec.Wh.spec))...
    max(max(L56eventSpec.Wh.spec)) max(max(L6eventSpec.Wh.spec))]);

figure;
subplot(6,3,1);
imagesc(taxis,log10(wavespec.freqs),L1eventSpec.Wh.spec');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmax*-1 cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.freqs(1)) log10(wavespec.freqs(end))],'--k');hold on;
title('L1 modZ-PSpec');

subplot(6,3,4);
imagesc(taxis,log10(wavespec.freqs),L23eventSpec.Wh.spec');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmax*-1 cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.freqs(1)) log10(wavespec.freqs(end))],'--k');hold on;
title('L2/3');

subplot(6,3,7);
imagesc(taxis,log10(wavespec.freqs),L4eventSpec.Wh.spec');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmax*-1 cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.freqs(1)) log10(wavespec.freqs(end))],'--k');hold on;
title('L4');

subplot(6,3,10);
imagesc(taxis,log10(wavespec.freqs),L5aeventSpec.Wh.spec');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmax*-1 cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.freqs(1)) log10(wavespec.freqs(end))],'--k');hold on;
title('L5a');

subplot(6,3,13);
imagesc(taxis,log10(wavespec.freqs),L56eventSpec.Wh.spec');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmax*-1 cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.freqs(1)) log10(wavespec.freqs(end))],'--k');hold on;
title('L5b');

subplot(6,3,16);
imagesc(taxis,log10(wavespec.freqs),L6eventSpec.Wh.spec');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmax*-1 cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.freqs(1)) log10(wavespec.freqs(end))],'--k');hold on;
title('L6');

% For Fractal
cmin = min([min(min(L1eventSpec.Wh.frac)) min(min(L23eventSpec.Wh.frac))...
    min(min(L4eventSpec.Wh.frac)) min(min(L5aeventSpec.Wh.frac))...
    min(min(L56eventSpec.Wh.frac)) min(min(L6eventSpec.Wh.frac))]);
cmax = max([max(max(L1eventSpec.Wh.frac)) max(max(L23eventSpec.Wh.frac))...
    max(max(L4eventSpec.Wh.frac)) max(max(L5aeventSpec.Wh.frac))...
    max(max(L56eventSpec.Wh.frac)) max(max(L6eventSpec.Wh.frac))]);

subplot(6,3,2);
imagesc(taxis,log10(wavespec.validfreq),L1eventSpec.Wh.frac');hold on;
colormap jet;
LogScale('y',10);
LogScale('c',10);
caxis([cmax/4*-1 cmax/4]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L1 Fractal');

subplot(6,3,5);
imagesc(taxis,log10(wavespec.validfreq),L23eventSpec.Wh.frac');hold on;
colormap jet;
LogScale('y',10);
LogScale('c',10);
caxis([cmax/4*-1 cmax/4]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L2/3');

subplot(6,3,8);
imagesc(taxis,log10(wavespec.validfreq),L4eventSpec.Wh.frac');hold on;
colormap jet;
LogScale('y',10);
LogScale('c',10);
caxis([cmax/4*-1 cmax/4]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L4');

subplot(6,3,11);
imagesc(taxis,log10(wavespec.validfreq),L5aeventSpec.Wh.frac');hold on;
colormap jet;
LogScale('y',10);
LogScale('c',10);
caxis([cmax/4*-1 cmax/4]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L5a');

subplot(6,3,14);
imagesc(taxis,log10(wavespec.validfreq),L56eventSpec.Wh.frac');hold on;
colormap jet;
LogScale('y',10);
LogScale('c',10);
caxis([cmax/4*-1 cmax/4]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L5b');

subplot(6,3,17);
imagesc(taxis,log10(wavespec.validfreq),L6eventSpec.Wh.frac');hold on;
colormap jet;
LogScale('y',10);
LogScale('c',10);
caxis([cmax/4*-1 cmax/4]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L6');

% for Osci
cmin = min([min(min(L1eventSpec.Wh.osci)) min(min(L23eventSpec.Wh.osci))...
    min(min(L4eventSpec.Wh.osci)) min(min(L5aeventSpec.Wh.osci))...
    min(min(L56eventSpec.Wh.osci)) min(min(L6eventSpec.Wh.osci))]);
cmax = max([max(max(L1eventSpec.Wh.osci)) max(max(L23eventSpec.Wh.osci))...
    max(max(L4eventSpec.Wh.osci)) max(max(L5aeventSpec.Wh.osci))...
    max(max(L56eventSpec.Wh.osci)) max(max(L6eventSpec.Wh.osci))]);

subplot(6,3,3);
imagesc(taxis,log10(wavespec.validfreq),L1eventSpec.Wh.osci');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmin cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L1 Osci');

subplot(6,3,6);
imagesc(taxis,log10(wavespec.validfreq),L23eventSpec.Wh.osci');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmin cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L2/3');

subplot(6,3,9);
imagesc(taxis,log10(wavespec.validfreq),L4eventSpec.Wh.osci');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmin cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L4');

subplot(6,3,12);
imagesc(taxis,log10(wavespec.validfreq),L5aeventSpec.Wh.osci');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmin cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L5a');

subplot(6,3,15);
imagesc(taxis,log10(wavespec.validfreq),L56eventSpec.Wh.osci');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmin cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L5b');

subplot(6,3,18);
imagesc(taxis,log10(wavespec.validfreq),L6eventSpec.Wh.osci');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmin cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L6');

NiceSave('Laminar_WhSpec',figfolder,baseName)

%% FIGURE 4: Pupil dilation-aligned Specs
cmin = min([min(min(L1eventSpec.Pup.spec)) min(min(L23eventSpec.Pup.spec))...
    min(min(L4eventSpec.Pup.spec)) min(min(L5aeventSpec.Pup.spec))...
    min(min(L56eventSpec.Pup.spec)) min(min(L6eventSpec.Pup.spec))]);
cmax = max([max(max(L1eventSpec.Pup.spec)) max(max(L23eventSpec.Pup.spec))...
    max(max(L4eventSpec.Pup.spec)) max(max(L5aeventSpec.Pup.spec))...
    max(max(L56eventSpec.Pup.spec)) max(max(L6eventSpec.Pup.spec))]);

figure;
subplot(6,3,1);
imagesc(taxis,log10(wavespec.freqs),L1eventSpec.Pup.spec');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmax*-1 cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.freqs(1)) log10(wavespec.freqs(end))],'--k');hold on;
title('L1 modZ-PSpec');

subplot(6,3,4);
imagesc(taxis,log10(wavespec.freqs),L23eventSpec.Pup.spec');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmax*-1 cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.freqs(1)) log10(wavespec.freqs(end))],'--k');hold on;
title('L2/3');

subplot(6,3,7);
imagesc(taxis,log10(wavespec.freqs),L4eventSpec.Pup.spec');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmax*-1 cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.freqs(1)) log10(wavespec.freqs(end))],'--k');hold on;
title('L4');

subplot(6,3,10);
imagesc(taxis,log10(wavespec.freqs),L5aeventSpec.Pup.spec');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmax*-1 cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.freqs(1)) log10(wavespec.freqs(end))],'--k');hold on;
title('L5a');

subplot(6,3,13);
imagesc(taxis,log10(wavespec.freqs),L56eventSpec.Pup.spec');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmax*-1 cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.freqs(1)) log10(wavespec.freqs(end))],'--k');hold on;
title('L5b');

subplot(6,3,16);
imagesc(taxis,log10(wavespec.freqs),L6eventSpec.Pup.spec');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmax*-1 cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.freqs(1)) log10(wavespec.freqs(end))],'--k');hold on;
title('L6');

% For Fractal
cmin = min([min(min(L1eventSpec.Pup.frac)) min(min(L23eventSpec.Pup.frac))...
    min(min(L4eventSpec.Pup.frac)) min(min(L5aeventSpec.Pup.frac))...
    min(min(L56eventSpec.Pup.frac)) min(min(L6eventSpec.Pup.frac))]);
cmax = max([max(max(L1eventSpec.Pup.frac)) max(max(L23eventSpec.Pup.frac))...
    max(max(L4eventSpec.Pup.frac)) max(max(L5aeventSpec.Pup.frac))...
    max(max(L56eventSpec.Pup.frac)) max(max(L6eventSpec.Pup.frac))]);

subplot(6,3,2);
imagesc(taxis,log10(wavespec.validfreq),L1eventSpec.Pup.frac');hold on;
colormap jet;
LogScale('y',10);
LogScale('c',10);
caxis([cmax*-1 cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L1 Fractal');

subplot(6,3,5);
imagesc(taxis,log10(wavespec.validfreq),L23eventSpec.Pup.frac');hold on;
colormap jet;
LogScale('y',10);
LogScale('c',10);
caxis([cmax*-1 cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L2/3');

subplot(6,3,8);
imagesc(taxis,log10(wavespec.validfreq),L4eventSpec.Pup.frac');hold on;
colormap jet;
LogScale('y',10);
LogScale('c',10);
caxis([cmax*-1 cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L4');

subplot(6,3,11);
imagesc(taxis,log10(wavespec.validfreq),L5aeventSpec.Pup.frac');hold on;
colormap jet;
LogScale('y',10);
LogScale('c',10);
caxis([cmax*-1 cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L5a');

subplot(6,3,14);
imagesc(taxis,log10(wavespec.validfreq),L56eventSpec.Pup.frac');hold on;
colormap jet;
LogScale('y',10);
LogScale('c',10);
caxis([cmax*-1 cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L5b');

subplot(6,3,17);
imagesc(taxis,log10(wavespec.validfreq),L6eventSpec.Pup.frac');hold on;
colormap jet;
LogScale('y',10);
LogScale('c',10);
caxis([cmax*-1 cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L6');

% for Osci
cmin = min([min(min(L1eventSpec.Pup.osci)) min(min(L23eventSpec.Pup.osci))...
    min(min(L4eventSpec.Pup.osci)) min(min(L5aeventSpec.Pup.osci))...
    min(min(L56eventSpec.Pup.osci)) min(min(L6eventSpec.Pup.osci))]);
cmax = max([max(max(L1eventSpec.Pup.osci)) max(max(L23eventSpec.Pup.osci))...
    max(max(L4eventSpec.Pup.osci)) max(max(L5aeventSpec.Pup.osci))...
    max(max(L56eventSpec.Pup.osci)) max(max(L6eventSpec.Pup.osci))]);

subplot(6,3,3);
imagesc(taxis,log10(wavespec.validfreq),L1eventSpec.Pup.osci');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmin cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L1 Osci');

subplot(6,3,6);
imagesc(taxis,log10(wavespec.validfreq),L23eventSpec.Pup.osci');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmin cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L2/3');

subplot(6,3,9);
imagesc(taxis,log10(wavespec.validfreq),L4eventSpec.Pup.osci');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmin cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L4');

subplot(6,3,12);
imagesc(taxis,log10(wavespec.validfreq),L5aeventSpec.Pup.osci');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmin cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L5a');

subplot(6,3,15);
imagesc(taxis,log10(wavespec.validfreq),L56eventSpec.Pup.osci');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmin cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L5b');

subplot(6,3,18);
imagesc(taxis,log10(wavespec.validfreq),L6eventSpec.Pup.osci');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmin cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L6');

NiceSave('Laminar_PupSpec',figfolder,baseName)

%% Laminar comodulograms
i = 1;
L1comodspec = NaN(size(laminarspec,2),size(laminarspec,2),6);
L1comodfrac = NaN(size(laminarfrac,2),size(laminarfrac,2),6);
L1comodosci = NaN(size(laminarosci,2),size(laminarosci,2),6);
for ii = 1:size(laminarspec,3)
    L1comodspec(:,:,ii) = corr(laminarspec(:,:,i),...
        laminarspec(:,:,ii),'type','spearman',...
        'rows','complete');
    L1comodfrac(:,:,ii) = corr(laminarfrac(:,:,i),...
        laminarfrac(:,:,ii),'type','spearman',...
        'rows','complete');
    L1comodosci(:,:,ii) = corr(laminarosci(:,:,i),...
        laminarosci(:,:,ii),'type','spearman',...
        'rows','complete');
end

i = 2;
L23comodspec = NaN(size(laminarspec,2),size(laminarspec,2),6);
L23comodfrac = NaN(size(laminarfrac,2),size(laminarfrac,2),6);
L23comodosci = NaN(size(laminarosci,2),size(laminarosci,2),6);
for ii = 1:size(laminarspec,3)
    L23comodspec(:,:,ii) = corr(laminarspec(:,:,i),...
        laminarspec(:,:,ii),'type','spearman',...
        'rows','complete');
    L23comodfrac(:,:,ii) = corr(laminarfrac(:,:,i),...
        laminarfrac(:,:,ii),'type','spearman',...
        'rows','complete');
    L23comodosci(:,:,ii) = corr(laminarosci(:,:,i),...
        laminarosci(:,:,ii),'type','spearman',...
        'rows','complete');
end

i = 3;
L4comodspec = NaN(size(laminarspec,2),size(laminarspec,2),6);
L4comodfrac = NaN(size(laminarfrac,2),size(laminarfrac,2),6);
L4comodosci = NaN(size(laminarosci,2),size(laminarosci,2),6);
for ii = 1:size(laminarspec,3)
    L4comodspec(:,:,ii) = corr(laminarspec(:,:,i),...
        laminarspec(:,:,ii),'type','spearman',...
        'rows','complete');
    L4comodfrac(:,:,ii) = corr(laminarfrac(:,:,i),...
        laminarfrac(:,:,ii),'type','spearman',...
        'rows','complete');
    L4comodosci(:,:,ii) = corr(laminarosci(:,:,i),...
        laminarosci(:,:,ii),'type','spearman',...
        'rows','complete');
end

i = 4;
L5acomodspec = NaN(size(laminarspec,2),size(laminarspec,2),6);
L5acomodfrac = NaN(size(laminarfrac,2),size(laminarfrac,2),6);
L5acomodosci = NaN(size(laminarosci,2),size(laminarosci,2),6);
for ii = 1:size(laminarspec,3)
    L5acomodspec(:,:,ii) = corr(laminarspec(:,:,i),...
        laminarspec(:,:,ii),'type','spearman',...
        'rows','complete');
    L5acomodfrac(:,:,ii) = corr(laminarfrac(:,:,i),...
        laminarfrac(:,:,ii),'type','spearman',...
        'rows','complete');
    L5acomodosci(:,:,ii) = corr(laminarosci(:,:,i),...
        laminarosci(:,:,ii),'type','spearman',...
        'rows','complete');
end

i = 5;
L56comodspec = NaN(size(laminarspec,2),size(laminarspec,2),6);
L56comodfrac = NaN(size(laminarfrac,2),size(laminarfrac,2),6);
L56comodosci = NaN(size(laminarosci,2),size(laminarosci,2),6);
for ii = 1:size(laminarspec,3)
    L56comodspec(:,:,ii) = corr(laminarspec(:,:,i),...
        laminarspec(:,:,ii),'type','spearman',...
        'rows','complete');
    L56comodfrac(:,:,ii) = corr(laminarfrac(:,:,i),...
        laminarfrac(:,:,ii),'type','spearman',...
        'rows','complete');
    L56comodosci(:,:,ii) = corr(laminarosci(:,:,i),...
        laminarosci(:,:,ii),'type','spearman',...
        'rows','complete');
end

i = 6;
L6comodspec = NaN(size(laminarspec,2),size(laminarspec,2),6);
L6comodfrac = NaN(size(laminarfrac,2),size(laminarfrac,2),6);
L6comodosci = NaN(size(laminarosci,2),size(laminarosci,2),6);
for ii = 1:size(laminarspec,3)
    L6comodspec(:,:,ii) = corr(laminarspec(:,:,i),...
        laminarspec(:,:,ii),'type','spearman',...
        'rows','complete');
    L6comodfrac(:,:,ii) = corr(laminarfrac(:,:,i),...
        laminarfrac(:,:,ii),'type','spearman',...
        'rows','complete');
    L6comodosci(:,:,ii) = corr(laminarosci(:,:,i),...
        laminarosci(:,:,ii),'type','spearman',...
        'rows','complete');
end

% Saving to struct
LayerSpectral.L1comodspec = L1comodspec; 
LayerSpectral.L1comodfrac = L1comodfrac; 
LayerSpectral.L1comodosci = L1comodosci;
LayerSpectral.L23comodspec = L23comodspec; 
LayerSpectral.L23comodfrac = L23comodfrac; 
LayerSpectral.L23comodosci = L23comodosci;
LayerSpectral.L4comodspec = L4comodspec; 
LayerSpectral.L4comodfrac = L4comodfrac; 
LayerSpectral.L4comodosci = L4comodosci;
LayerSpectral.L5acomodspec = L5acomodspec; 
LayerSpectral.L5acomodfrac = L5acomodfrac; 
LayerSpectral.L5acomodosci = L5acomodosci;
LayerSpectral.L56comodspec = L56comodspec; 
LayerSpectral.L56comodfrac = L56comodfrac; 
LayerSpectral.L56comodosci = L56comodosci;
LayerSpectral.L6comodspec = L6comodspec; 
LayerSpectral.L6comodfrac = L6comodfrac; 
LayerSpectral.L6comodosci = L6comodosci;

% Finally saving all...
save(savefile,'-v7.3','LayerSpectral');

%% FIGURE 5: Laminar CoMOD Spec
cmin = min([min(min(min(L1comodspec))) min(min(min(L23comodspec)))...
    min(min(min(L4comodspec))) min(min(min(L5acomodspec)))...
    min(min(min(L56comodspec))) min(min(min(L6comodspec)))]);
cmax = max([max(max(max(L1comodspec))) max(max(max(L23comodspec)))...
    max(max(max(L4comodspec))) max(max(max(L5acomodspec)))...
    max(max(max(L56comodspec))) max(max(max(L6comodspec)))]);

figure;
for i = 1:size(L1comodspec,3)
    subplot(6,6,i*1);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L1comodspec(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
    xlim(log10([1 100])); ylim(log10([1 100]));
end

for i = 1:size(L23comodspec,3)
    subplot(6,6,i+6);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L23comodspec(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
    xlim(log10([1 100])); ylim(log10([1 100]));
end

for i = 1:size(L4comodspec,3)
    subplot(6,6,i+12);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L4comodspec(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
    xlim(log10([1 100])); ylim(log10([1 100]));
end

for i = 1:size(L5acomodspec,3)
    subplot(6,6,i+18);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L5acomodspec(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
    xlim(log10([1 100])); ylim(log10([1 100]));
end

for i = 1:size(L56comodspec,3)
    subplot(6,6,i+24);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L56comodspec(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
    xlim(log10([1 100])); ylim(log10([1 100]));
end

for i = 1:size(L6comodspec,3)
    subplot(6,6,i+30);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L6comodspec(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
    xlim(log10([1 100])); ylim(log10([1 100]));
end

NiceSave('LaminarCoMOD_Spec',figfolder,baseName)

%% FIGURE:
cmin = min([min(min(min(L1comodfrac))) min(min(min(L23comodfrac)))...
    min(min(min(L4comodfrac))) min(min(min(L5acomodfrac)))...
    min(min(min(L56comodfrac))) min(min(min(L6comodfrac)))]);
cmax = max([max(max(max(L1comodfrac))) max(max(max(L23comodfrac)))...
    max(max(max(L4comodfrac))) max(max(max(L5acomodfrac)))...
    max(max(max(L56comodfrac))) max(max(max(L6comodfrac)))]);

figure;
for i = 1:size(L1comodfrac,3)
    subplot(6,6,i*1);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L1comodfrac(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
    xlim(log10([1 100])); ylim(log10([1 100]));
end

for i = 1:size(L23comodfrac,3)
    subplot(6,6,i+6);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L23comodfrac(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
    xlim(log10([1 100])); ylim(log10([1 100]));
end

for i = 1:size(L4comodfrac,3)
    subplot(6,6,i+12);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L4comodfrac(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
    xlim(log10([1 100])); ylim(log10([1 100]));
end

for i = 1:size(L5acomodfrac,3)
    subplot(6,6,i+18);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L5acomodfrac(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
    xlim(log10([1 100])); ylim(log10([1 100]));
end

for i = 1:size(L56comodfrac,3)
    subplot(6,6,i+24);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L56comodfrac(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
    xlim(log10([1 100])); ylim(log10([1 100]));
end

for i = 1:size(L6comodfrac,3)
    subplot(6,6,i+30);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L6comodfrac(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
    xlim(log10([1 100])); ylim(log10([1 100]));
end

NiceSave('LaminarCoMOD_Frac',figfolder,baseName)

%% FIGURE:
cmin = min([min(min(min(L1comodosci))) min(min(min(L23comodosci)))...
    min(min(min(L4comodosci))) min(min(min(L5acomodosci)))...
    min(min(min(L56comodosci))) min(min(min(L6comodosci)))]);
cmax = max([max(max(max(L1comodosci))) max(max(max(L23comodosci)))...
    max(max(max(L4comodosci))) max(max(max(L5acomodosci)))...
    max(max(max(L56comodosci))) max(max(max(L6comodosci)))]);

figure;
for i = 1:size(L1comodosci,3)
    subplot(6,6,i*1);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L1comodosci(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
    xlim(log10([1 100])); ylim(log10([1 100]));
end

for i = 1:size(L23comodosci,3)
    subplot(6,6,i+6);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L23comodosci(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
    xlim(log10([1 100])); ylim(log10([1 100]));
end

for i = 1:size(L4comodosci,3)
    subplot(6,6,i+12);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L4comodosci(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
    xlim(log10([1 100])); ylim(log10([1 100]));
end

for i = 1:size(L5acomodosci,3)
    subplot(6,6,i+18);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L5acomodosci(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
    xlim(log10([1 100])); ylim(log10([1 100]));
end

for i = 1:size(L56comodosci,3)
    subplot(6,6,i+24);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L56comodosci(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
    xlim(log10([1 100])); ylim(log10([1 100]));
end

for i = 1:size(L6comodosci,3)
    subplot(6,6,i+30);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L6comodosci(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
    xlim(log10([1 100])); ylim(log10([1 100]));
end

NiceSave('LaminarCoMOD_Osci',figfolder,baseName)