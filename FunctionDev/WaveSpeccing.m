basePath = pwd;
baseName = bz_BasenameFromBasepath(basePath);
savefolder = fullfile(basePath,'WaveSpec2');
if (~exist(savefolder,'dir'))
    mkdir(savefolder)
end
%%
lfp = bz_GetLFP('all','basepath',basePath,'noPrompts',true);
lfp = bz_DownsampleLFP(lfp,4); %version 1: downfactor 5

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
tempwavespec = bz_WaveSpec_GPU(lfp,'ncyc',15,'showprogress',true); %version 1: ncyc 5
% toc
% profile off
% profile viewer

%%
movingwin = round([0.5 0.1].*lfp.samplingRate);
nwin = floor((size(lfp.data,1) - movingwin(1))/movingwin(2));

for i = 1:length(lfp.channels)
    
    %
    sig = zeros(movingwin(1),nwin);
    timestamp = zeros(nwin,1);
    for ii = 1 : nwin
        idx = [ceil((ii-1)*movingwin(2))+1 : ceil((ii-1)*movingwin(2))+movingwin(1)];
        sig(:,ii) = double(lfp.data(idx,i));
        %figure out timestamp associated with window i
        timestamp(ii) = mean(lfp.timestamps(idx));
    end
    
    % 
    lfpfilename = fullfile(savefolder,[baseName,'.',num2str(lfp.channels(i)),'.WaveSpec2.lfp.mat']);
    wavespec.data = tempwavespec.data(:,:,i);
    wavespec.timestamps = tempwavespec.timestamps;
    wavespec.freqs = tempwavespec.freqs;
    wavespec.nfreqs = tempwavespec.nfreqs;
    wavespec.samplingRate = tempwavespec.samplingRate;
    wavespec.channels = lfp.channels(i);
    wavespec.filterparms = tempwavespec.filterparms;
    wavespec.Frac = amri_sig_fractal_gpu(sig,lfp.samplingRate,'detrend',1);
    wavespec.Frac.timestamps = timestamp;
    
    save(lfpfilename,'-v7.3','wavespec');
    
end