basePath = pwd;
baseName = bz_BasenameFromBasepath(basePath);
savefolder = fullfile(basePath,'WaveSpec');
if (~exist(savefolder,'dir'))
    mkdir(savefolder)
end
%%
lfp = bz_GetLFP('all','basepath',basePath,'noPrompts',true);
lfp = bz_DownsampleLFP(lfp,2.5); %version 1: downfactor 5

%%
load(fullfile(basePath,[baseName,'.MergePoints.events.mat']),'MergePoints');
sidx = find(startsWith(MergePoints.foldernames,"Spont"));
sponttimes = [MergePoints.timestamps(sidx(1),1) MergePoints.timestamps(sidx(end),2)];

spontidx = find(lfp.timestamps < sponttimes(2));
lfp.timestamps = lfp.timestamps(spontidx);
lfp.data = lfp.data(spontidx,:);

%%
% profile on
% tic
tempwavespec = bz_WaveSpec_GPU(lfp,'frange',[0.1 250],'nfreqs',100,'showprogress',true); %version 1: ncyc 5
% toc
% profile off
% profile viewer

%%
for i = 1:length(lfp.channels)
    lfpfilename = fullfile(savefolder,[baseName,'.',num2str(lfp.channels(i)),'.WaveSpec.lfp.mat']);
    wavespec.data = tempwavespec.data(:,:,i);
    wavespec.timestamps = tempwavespec.timestamps;
    wavespec.freqs = tempwavespec.freqs;
    wavespec.nfreqs = tempwavespec.nfreqs;
    wavespec.samplingRate = tempwavespec.samplingRate;
    wavespec.channels = lfp.channels(i);
    wavespec.filterparms = tempwavespec.filterparms;
    
    save(lfpfilename,'-v7.3','wavespec');
end