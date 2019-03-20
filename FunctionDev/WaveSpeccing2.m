basePath = pwd;
baseName = bz_BasenameFromBasepath(basePath);

sessionInfo = bz_getSessionInfo(basePath, 'noPrompts', true);
datSampleRate = sessionInfo.rates.wideband;
datfilename = fullfile(basePath,[baseName,'.dat']);
channels = sessionInfo.channels;

%%
downfactor = 40;
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

spontidx = find(datlfp.timestamps < sponttimes(2));
datlfp.timestamps = datlfp.timestamps(spontidx);
datlfp.data = datlfp.data(spontidx,:);

%%
tempwavespec = bz_WaveSpec_GPU(datlfp,'frange',[0.1 250],'nfreqs',200,'showprogress',true); %version 1: ncyc 5

%%

for i = 1:length(lfp.channels)
    i 
    % 
    lfpfilename = fullfile(savefolder,[baseName,'.',num2str(lfp.channels(i)),'.WaveSpec2.lfp.mat']);
    wavespec.data = tempwavespec.data(:,:,i);
    wavespec.timestamps = tempwavespec.timestamps;
    wavespec.freqs = tempwavespec.freqs;
    wavespec.nfreqs = tempwavespec.nfreqs;
    wavespec.samplingRate = tempwavespec.samplingRate;
    wavespec.channels = lfp.channels(i);
    wavespec.filterparms = tempwavespec.filterparms;
    
    [wavespec.frac,wavespec.osci,wavespec.validfreq] = WaveIRASA(wavespec);
    
    save(lfpfilename,'-v7.3','wavespec');
    
end
