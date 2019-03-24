basePath = pwd;
baseName = bz_BasenameFromBasepath(basePath);
savefolder = fullfile(basePath,'WaveSpec');
if (~exist(savefolder,'dir'))
    mkdir(savefolder)
end

%%
sessionInfo = bz_getSessionInfo(basePath, 'noPrompts', true);
datSampleRate = sessionInfo.rates.wideband;
datfilename = fullfile(basePath,[baseName,'.dat']);
channels = sessionInfo.channels;

%%
downfactor = 25;
datlfp.data = bz_LoadBinary(datfilename,...
              'frequency',datSampleRate,'nchannels',sessionInfo.nChannels,...
              'channels',channels+1,'downsample',downfactor);

datlfp.samplingRate = datSampleRate./downfactor;
datlfp.timestamps = [0:(length(datlfp.data)-1)]'/datlfp.samplingRate;  %To be overwritten later...
datlfp.channels = channels;

%%
load(fullfile(basePath,[baseName,'.MergePoints.events.mat']),'MergePoints');
sidx = find(startsWith(MergePoints.foldernames,"Spont"));
sponttimes = [MergePoints.timestamps(sidx(1),1) MergePoints.timestamps(sidx(end),2)];
%sponttimes = [MergePoints.timestamps(sidx(1),1) MergePoints.timestamps(sidx(1),2)/15];

spontidx = find(datlfp.timestamps < sponttimes(2));
datlfp.timestamps = datlfp.timestamps(spontidx);
datlfp.data = datlfp.data(spontidx,:);

%%
tempdatlfp = datlfp;
tempdatlfp.data = tempdatlfp.data(:,1:16);
tempdatlfp.channels = tempdatlfp.channels(1:16);
tempwavespec = bz_WaveSpec_GPU(tempdatlfp,'frange',[0.1 400],'nfreqs',100,'showprogress',true); %version 1: ncyc 5

%%

for i = 1:length(tempdatlfp.channels)
    i 
    % 
    lfpfilename = fullfile(savefolder,[baseName,'.',num2str(tempdatlfp.channels(i)),'.WaveSpec.lfp.mat']);
    wavespec.data = tempwavespec.data(:,:,i);
    wavespec.timestamps = tempwavespec.timestamps;
    wavespec.freqs = tempwavespec.freqs;
    wavespec.nfreqs = tempwavespec.nfreqs;
    wavespec.samplingRate = tempwavespec.samplingRate;
    wavespec.channels = tempdatlfp.channels(i);
    wavespec.filterparms = tempwavespec.filterparms;
       
    save(lfpfilename,'-v7.3','-nocompression','wavespec');
   
end

%%
tempdatlfp = datlfp;
tempdatlfp.data = tempdatlfp.data(:,17:32);
tempdatlfp.channels = tempdatlfp.channels(17:32);
tempwavespec = bz_WaveSpec_GPU(tempdatlfp,'frange',[0.1 400],'nfreqs',100,'showprogress',true); %version 1: ncyc 5

%%

for i = 1:length(tempdatlfp.channels)
    i 
    % 
    lfpfilename = fullfile(savefolder,[baseName,'.',num2str(tempdatlfp.channels(i)),'.WaveSpec.lfp.mat']);
    wavespec.data = tempwavespec.data(:,:,i);
    wavespec.timestamps = tempwavespec.timestamps;
    wavespec.freqs = tempwavespec.freqs;
    wavespec.nfreqs = tempwavespec.nfreqs;
    wavespec.samplingRate = tempwavespec.samplingRate;
    wavespec.channels = tempdatlfp.channels(i);
    wavespec.filterparms = tempwavespec.filterparms;
       
    save(lfpfilename,'-v7.3','-nocompression','wavespec');
   
end

%%
tempdatlfp = datlfp;
tempdatlfp.data = tempdatlfp.data(:,33:48);
tempdatlfp.channels = tempdatlfp.channels(33:48);
tempwavespec = bz_WaveSpec_GPU(tempdatlfp,'frange',[0.1 400],'nfreqs',100,'showprogress',true); %version 1: ncyc 5

%%

for i = 1:length(tempdatlfp.channels)
    i 
    % 
    lfpfilename = fullfile(savefolder,[baseName,'.',num2str(tempdatlfp.channels(i)),'.WaveSpec.lfp.mat']);
    wavespec.data = tempwavespec.data(:,:,i);
    wavespec.timestamps = tempwavespec.timestamps;
    wavespec.freqs = tempwavespec.freqs;
    wavespec.nfreqs = tempwavespec.nfreqs;
    wavespec.samplingRate = tempwavespec.samplingRate;
    wavespec.channels = tempdatlfp.channels(i);
    wavespec.filterparms = tempwavespec.filterparms;
       
    save(lfpfilename,'-v7.3','-nocompression','wavespec');
   
end

%%
tempdatlfp = datlfp;
tempdatlfp.data = tempdatlfp.data(:,49:64);
tempdatlfp.channels = tempdatlfp.channels(49:64);
tempwavespec = bz_WaveSpec_GPU(tempdatlfp,'frange',[0.1 400],'nfreqs',100,'showprogress',true); %version 1: ncyc 5

%%

for i = 1:length(tempdatlfp.channels)
    i 
    % 
    lfpfilename = fullfile(savefolder,[baseName,'.',num2str(tempdatlfp.channels(i)),'.WaveSpec.lfp.mat']);
    wavespec.data = tempwavespec.data(:,:,i);
    wavespec.timestamps = tempwavespec.timestamps;
    wavespec.freqs = tempwavespec.freqs;
    wavespec.nfreqs = tempwavespec.nfreqs;
    wavespec.samplingRate = tempwavespec.samplingRate;
    wavespec.channels = tempdatlfp.channels(i);
    wavespec.filterparms = tempwavespec.filterparms;
       
    save(lfpfilename,'-v7.3','-nocompression','wavespec');
   
end