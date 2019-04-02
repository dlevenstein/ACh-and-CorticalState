%%
baseGroup = ['/gpfs/data/rudylab/William/171006_WT_EM1M3';...
    '/gpfs/data/rudylab/William/171206_WT_EM1M3';...
    '/gpfs/data/rudylab/William/171209_WT_EM1M3';...
    '/gpfs/data/rudylab/William/180209_WT_EM1M3';...
    '/gpfs/data/rudylab/William/180211_WT_EM1M3';...
    '/gpfs/data/rudylab/William/180213_WT_EM1M3';...
    '/gpfs/data/rudylab/William/180214_WT_EM1M3';...
    '/gpfs/data/rudylab/William/180603_WT_EM1M3';...
    '/gpfs/data/rudylab/William/180605_WT_EM1M3';...
    '/gpfs/data/rudylab/William/180607_WT_EM1M3';...
    '/gpfs/data/rudylab/William/180703_WT_EM1M3';...
    '/gpfs/data/rudylab/William/180706_WT_EM1M3';...
    '/gpfs/data/rudylab/William/180708_WT_EM1M3';...
    '/gpfs/data/rudylab/William/180710_WT_EM1M3'];

baseGroup = ['/gpfs/data/rudylab/William/171211_KO_EM1M3';...
    '/gpfs/data/rudylab/William/180208_KO_EM1M3';...
    '/gpfs/data/rudylab/William/180210_KO_EM1M3';...
    '/gpfs/data/rudylab/William/180212_KO_EM1M3';...
    '/gpfs/data/rudylab/William/180525_KO_EM1M3';...
    '/gpfs/data/rudylab/William/180530_KO_EM1M3';...
    '/gpfs/data/rudylab/William/180601_KO_EM1M3';...
    '/gpfs/data/rudylab/William/180602_KO_EM1M3';...
    '/gpfs/data/rudylab/William/180608_KO_EM1M3';...
    '/gpfs/data/rudylab/William/180704_KO_EM1M3';...
    '/gpfs/data/rudylab/William/180705_KO_EM1M3';...
    '/gpfs/data/rudylab/William/180709_KO_EM1M3';...
    '/gpfs/data/rudylab/William/180711_KO_EM1M3'];

%%
basePath = pwd;
%baseName = 'WT_EM1M3';
baseName = 'KO_EM1M3';

savefile = fullfile(basePath,[baseName,'.BehaviorAnalysis2.mat']);
%savefile = fullfile(basePath,[baseName,'.ColLFPSpectralAnalysis.mat']);
%savefile = fullfile(basePath,[baseName,'.LamLFPSpectralAnalysis.mat']);
%savefile = fullfile(basePath,[baseName,'.PowerSpectrumSlope.lfp.mat']);
%savefile = fullfile(basePath,[baseName,'.PSSBehaviorAnalysis2.mat']);

%% Rescale and collect data
% For Behavior
groupBehavior = []; newedges = -20:0.1:20;
groupBehaviorInts.whints_pupil = [];
groupBehaviorInts.whints_pupildt = [];
groupBehaviorInts.whints_pupilphase = [];
groupBehaviorInts.whints_pupilamp = [];
groupBehaviorInts.whints_durs = [];
groupBehaviorInts.whints_maxamp = [];
groupBehaviorInts.whthresh = [];

% For Columnar Spectral
groupDepth = [];
groupColLFPSpectral = [];

% For Laminar Spectral
% groupLamLFPSpectral = [];

% For PSS
% groupDepth = [];
% groupPSS = [];

% For PSS Behavior
% groupPSSBehavior = [];
% puplockedPSS.data = [];
% puplockedPSS.timestamps = [];
% puplockedPSS.pupwhdiff = [];
% pupsorts.Pupwhdiff = [];
% pupsorts.dP = [];
% pupsorts.dur = [];
% pupsorts.peak = [];
% pupsorts.phase = [];
% pupsorts.pow = [];
% 
% timelockedPSS.data = [];
% timelockedPSS.timestamps = [];
% timelockedPSS.phases = [];
% timelockedPSS.highpupil = [];
% whisksorts.phase = [];
% whisksorts.dur = [];
% whisksorts.amp = [];
% whisksorts.maxpss = [];

for i = 1:size(baseGroup,1)
    
    basename = bz_BasenameFromBasepath(baseGroup(i,:));
    
    %For Behavior analysis
    load(fullfile(baseGroup(i,:),[basename,'.BehaviorAnalysis.mat']));
    
    groupBehaviorInts.whints_pupil = cat(1,groupBehaviorInts.whints_pupil,PupEMG.whints_pupil);
    groupBehaviorInts.whints_pupildt = cat(1,groupBehaviorInts.whints_pupildt,PupEMG.whints_pupildt);
    groupBehaviorInts.whints_pupilphase = cat(1,groupBehaviorInts.whints_pupilphase,PupEMG.whints_pupilphase);
    groupBehaviorInts.whints_pupilamp = cat(1,groupBehaviorInts.whints_pupilamp,PupEMG.whints_pupilamp);
    groupBehaviorInts.whints_durs = cat(1,groupBehaviorInts.whints_durs,PupEMG.whints_durs);
    groupBehaviorInts.whints_maxamp = cat(1,groupBehaviorInts.whints_maxamp,PupEMG.whints_maxamp);
    groupBehaviorInts.whthresh = cat(1,groupBehaviorInts.whthresh,PupEMG.EMGhist.threshold);
    
    PupEMG.pupACG.ACG = interp1(PupEMG.pupACG.tlag,PupEMG.pupACG.ACG,newedges,'nearest');
    PupEMG.pupACG.tlag = newedges;
    
    PupEMG.pupilEMGcorr.xcorr = interp1(PupEMG.pupilEMGcorr.corrlags,PupEMG.pupilEMGcorr.xcorr,newedges,'nearest');
    PupEMG.pupilEMGcorr.corrlags = newedges;
    
    PupEMG.pwCCG.pupil.WhOn = interp1(PupEMG.pwCCG.t_lag,PupEMG.pwCCG.pupil.WhOn,newedges,'nearest');
    PupEMG.pwCCG.pupilxy.WhOn = interp1(PupEMG.pwCCG.t_lag,PupEMG.pwCCG.pupilxy.WhOn,newedges,'nearest');
    PupEMG.pwCCG.EMG.WhOn = interp1(PupEMG.pwCCG.t_lag,PupEMG.pwCCG.EMG.WhOn,newedges,'nearest');
    PupEMG.pwCCG.t_lag = newedges;
    
    PupEMG.wpCCG.pupil.PupOn = interp1(PupEMG.wpCCG.t_lag,PupEMG.wpCCG.pupil.PupOn,newedges,'nearest');
    PupEMG.wpCCG.pupilxy.PupOn = interp1(PupEMG.wpCCG.t_lag,PupEMG.wpCCG.pupilxy.PupOn,newedges,'nearest');
    PupEMG.wpCCG.EMG.PupOn = interp1(PupEMG.wpCCG.t_lag,PupEMG.wpCCG.EMG.PupOn,newedges,'nearest');
    PupEMG.wpCCG.t_lag = newedges;
    
    PupEMG.PhaseAmpCoup = rmfield(PupEMG.PhaseAmpCoup,{'sig1amp','sig2amp'});
    PupEMG.EMGhist = rmfield(PupEMG.EMGhist,{'threshold'});
    
    PupEMG = rmfield(PupEMG,{'whints_pupil','whints_pupildt','whints_pupilphase','whints_pupilamp',...
        'whints_durs','whints_maxamp'});
    
    groupBehavior = bz_CollapseStruct([groupBehavior PupEMG],3,'justcat',true);
%     
%     % For Columnar LFP Spectral analysis
%         sessionInfo = bz_getSessionInfo(baseGroup(i,:),'noPrompts',true);
%         groupDepth = bz_CollapseStruct([groupDepth rescaleCx(baseGroup(i,:))],3,'justcat',true);
%         
%         load(fullfile(baseGroup(i,:),[basename,'.ColumnarSpectralAnalysis.mat']));
%         if isfield(ColumnSpectral,{'columnspec_NWh','columnspec_Wh','columnspec_loWh','columnspec_hiWh'})
%         ColumnSpectral = rmfield(ColumnSpectral,{'columnspec_NWh','columnspec_Wh','columnspec_loWh','columnspec_hiWh'});
%         else
%         end
%         groupColLFPSpectral = bz_CollapseStruct([groupColLFPSpectral ColumnSpectral],3,'justcat',true);
%     
    % For Laminar LFP Spectral analysis
    %     load(fullfile(baseGroup(i,:),[basename,'.LaminarSpectralAnalysis.lfp.mat']));
    %         load(fullfile(baseGroup(i,:),[basename,'.LaminarSpectralAnalysis.mat']));
    %         groupLamLFPSpectral = bz_CollapseStruct([groupLamLFPSpectral LayerSpectral],3,'justcat',true);
    
    % For PSS analysis
    %     sessionInfo = bz_getSessionInfo(baseGroup(i,:),'noPrompts',true);
    %     groupDepth = bz_CollapseStruct([groupDepth rescaleCx(baseGroup(i,:))],3,'justcat',true);
    %     load(fullfile(baseGroup(i,:),[basename,'.PowerSpectrumSlope.lfp.mat']));
    %     PSpecSlope.Shortwin = rmfield(PSpecSlope.Shortwin,{'timestamps','PSS','Osci','Deltaosci','Thetaosci','Gammaosci'});
    %     PSpecSlope.Longwin = rmfield(PSpecSlope.Longwin,{'timestamps','PSS','Osci','Deltaosci','Thetaosci','Gammaosci'});
    %     groupPSS = bz_CollapseStruct([groupPSS PSpecSlope],3,'justcat',true);
    
    % For PSS-Behavior analysis
    %     load(fullfile(baseGroup(i,:),[basename,'.PSSBehaviorAnalysis.mat']));
%             load(fullfile(baseGroup(i,:),[basename,'.PSSBehaviorAnalysis2.mat']));
%     
%             puplockedPSS.data = cat(2,puplockedPSS.data,PSSBehavior.puplockedPSS.data);
%             puplockedPSS.timestamps = cat(2,puplockedPSS.timestamps,PSSBehavior.puplockedPSS.timestamps);
%             puplockedPSS.pupwhdiff = cat(2,puplockedPSS.pupwhdiff,PSSBehavior.puplockedPSS.pupwhdiff);
%             pupsorts.Pupwhdiff = cat(1,pupsorts.Pupwhdiff,PSSBehavior.pupidx_pupwhdiff);
%             pupsorts.dP = cat(1,pupsorts.dP,PSSBehavior.pupidx_dPpeak);
%             pupsorts.dur = cat(1,pupsorts.dur,PSSBehavior.pupidx_dur);
%             pupsorts.peak = cat(2,pupsorts.peak,PSSBehavior.pupidx_amp);
%             pupsorts.phase = cat(1,pupsorts.phase,PSSBehavior.pupidx_phase);
%             pupsorts.pow = cat(1,pupsorts.pow,PSSBehavior.pupidx_pow);
%     
%             timelockedPSS.data = cat(2,timelockedPSS.data,PSSBehavior.timelockedPSS.data);
%             timelockedPSS.timestamps = cat(2,timelockedPSS.timestamps,PSSBehavior.timelockedPSS.timestamps);
%             timelockedPSS.phases = cat(2,timelockedPSS.phases,PSSBehavior.timelockedPSS.phases);
%             timelockedPSS.highpupil = cat(2,timelockedPSS.highpupil,PSSBehavior.timelockedPSS.highpupil);
%             whisksorts.phase = cat(1,whisksorts.phase,PSSBehavior.whiskidx_phase);
%             whisksorts.dur = cat(1,whisksorts.dur,PSSBehavior.whiskidx_dur);
%             whisksorts.amp = cat(1,whisksorts.amp,PSSBehavior.whiskidx_amp);
%             whisksorts.maxpss = cat(1,whisksorts.maxpss,PSSBehavior.whiskidx_maxpss);
%     
%             PSSBehavior = rmfield(PSSBehavior,{'timelockedPSS','whiskidx_phase','whiskidx_dur','whiskidx_amp',...
%                 'whiskidx_maxpss','puplockedPSS','pupidx_pupwhdiff','pupidx_dPpeak','pupidx_dur',...
%                 'pupidx_amp','pupidx_phase','pupidx_pow'});
%             groupPSSBehavior = bz_CollapseStruct([groupPSSBehavior PSSBehavior],3,'justcat',true);
%     
    % For Layer ID analysis
    % For UP/DOWN analysis
    % For Unit analysis
end

% Saving group structs
save(savefile,'groupBehavior','groupBehaviorInts');
%save(savefile,'groupPSS','groupDepth');
%save(savefile,'groupPSSBehavior','puplockedPSS','pupsorts','timelockedPSS','whisksorts');
%save(savefile,'groupColLFPSpectral','groupDepth');
%save(savefile,'groupLamLFPSpectral');
