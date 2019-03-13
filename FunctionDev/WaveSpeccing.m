basePath = pwd;
baseName = bz_BasenameFromBasepath(basePath);
savefolder = fullfile(basePath,'WaveSpec');

lfp = bz_GetLFP('all','basepath',basePath,'noPrompts',true);
lfp = bz_DownsampleLFP(lfp,5);
% profile on
% tic
tempwavespec = bz_WaveSpec_GPU(lfp,'showprogress',true);
% toc
% profile off
% profile viewer

for i = 1:length(lfp.channels)
    
    lfpfilename = fullfile(savefolder,[baseName,'.',num2str(lfp.channels(i)),'.WaveSpec.lfp.mat']);
    wavespec.data = tempwavespec.data(:,:,i);
    wavespec.timestamps = tempwavespec.timestamps;
    wavespec.freqs = tempwavespec.freqs;
    wavespec.nfreqs = tempwavespec.nfreqs;
    wavespec.samplingRate = tempwavespec.samplingRate;
    wavespec.channels = lfp.channels(i);
    wavespec.filterparms = tempwavespec.filterparms;
    
    save(lfpfilename,'wavespec');
    
end