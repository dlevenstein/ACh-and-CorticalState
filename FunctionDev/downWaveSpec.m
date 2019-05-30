basePath = pwd;
baseName = bz_BasenameFromBasepath(basePath);
savefolder = fullfile(basePath,'WaveSpec_Downsampled');
if (~exist(savefolder,'dir'))
    mkdir(savefolder)
end

%%
downfactor = 8;
for x = 0:63
    x
    load(fullfile(fullfile(basePath,'WaveSpec'),[baseName,'.',num2str(x),'.WaveSpec.lfp.mat']));
    
    wavespec.data = downsample(wavespec.data,downfactor);
    wavespec.samplingRate = wavespec.samplingRate./downfactor;
    wavespec.timestamps = [0:(size(wavespec.data,1)-1)]'/wavespec.samplingRate;
    
    lfpfilename = fullfile(savefolder,[baseName,'.',num2str(x),'.WaveSpec.lfp.mat']);
    save(lfpfilename,'-v7.3','wavespec');
    
end