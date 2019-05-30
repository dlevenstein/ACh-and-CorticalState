function [ ] = AnalysisXXXXXXXX(basePath,figfolder)
% Date XX/XX/20XX
%
%Question: 
%
%Plots
%-
%-
%
%% Load Header
%Initiate Paths
reporoot = '/home/dlevenstein/ProjectRepos/ACh-and-CorticalState/';
%reporoot = '/Users/dlevenstein/Project Repos/ACh-and-CorticalState/';
basePath = '/mnt/proraidDL/Database/WMData/AChPupil/171209_WT_EM1M3/';
%basePath = '/mnt/proraidDL/Database/WMData/AChPupil/180706_WT_EM1M3/';
%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath,'noPrompts',true);


%% Load the LFP if needed

lfpchan = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID;
downsamplefactor = 1;
lfp = bz_GetLFP(lfpchan,...
    'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);
%Noralize the LFP
%lfp.data = NormToInt(single(lfp.data),'modZ', SleepState.ints.NREMstate,lfp.samplingRate);

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

% Filtered Pupil
lowfilter = [0.01 0.1];
%highfilter = [0.3 0.8];

pupil4filter = pupildilation;
pupilcycle = bz_Filter(pupil4filter,'passband',lowfilter,'filter' ,'fir1','order',3);
%highpupildata = bz_Filter(pupil4filter,'passband',highfilter,'filter' ,'fir1');
pupilcycle.pupthresh = -0.8;
pupilcycle.highpup = log10(pupilcycle.amp)>pupilcycle.pupthresh; 

% EMG
EMGwhisk = bz_LoadBehavior(basePath,'EMGwhisk');

%% Load PSS and aligning behavior
%load(fullfile(basePath,[baseName,'.PowerSpectrumSlope.mat']));
load(fullfile(basePath,[baseName,'.PowerSpectrumSlope.lfp.mat']));

depthinfo = rescaleCx(basePath);
L56idx = find(depthinfo.ndepth >= 0.6 & depthinfo.ndepth <= 0.9);
L56idx = depthinfo.channels(L56idx);

tempPSSEMGcorr = PSpecSlope.Shortwin.PSScorr.EMG(L56idx+1);
bestchan = find(tempPSSEMGcorr == max(tempPSSEMGcorr));

PSS.data = PSpecSlope.Shortwin.PSS(:,L56idx(bestchan)+1);
PSS.timestamps = PSpecSlope.Shortwin.timestamps;
PSS.samplingRate = 1/mean(diff(PSS.timestamps));
PSS.winsize = PSpecSlope.Shortwin.movingwin(1);
