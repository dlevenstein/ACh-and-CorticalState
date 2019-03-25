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
savefile = fullfile(basePath,[baseName,'.ColumnarSpectralAnalysis.mat']);
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
    if ~isempty(tidx)
        pup_on(i) = tidx;
    else
    end
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
ColumnSpectral.freqs = wavespec.freqs;
ColumnSpectral.validfreqs = wavespec.validfreq;
ColumnSpectral.columnspec_all = columnspec_all;
ColumnSpectral.columnspec_NWh = columnspec_NWh;
ColumnSpectral.columnspec_lNWh = columnspec_lNWh;
ColumnSpectral.columnspec_Wh = columnspec_Wh;
ColumnSpectral.columnspec_sWh = columnspec_sWh;
ColumnSpectral.columnspec_lWh = columnspec_lWh;
ColumnSpectral.columnspec_loWh = columnspec_loWh;
ColumnSpectral.columnspec_hiWh = columnspec_hiWh;
ColumnSpectral.columnspec_loP = columnspec_loP;
ColumnSpectral.columnspec_hiP = columnspec_hiP;

save(savefile,'-v7.3','ColumnSpectral');

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