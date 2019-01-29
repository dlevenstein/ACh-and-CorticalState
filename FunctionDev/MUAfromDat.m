function [ output ] = MUAfromDat( basePath,varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%INPUTS
%   basePath
%
%   (options)
%       'channels'      (default: all)   (0-indexed, like neuroscope)
%       'saveMat'       (default: true) save a .mat file
%                                   basePath/baseName.MUAGamma.lfp.mat
%       'filterbounds'  (default: [500 5000])
%       'movingNorm'    (default: false)
%       'MUAsmoothwin'  (default: 0.005)
%       'usepeaks'      (default: false)
%       'SHOWFIG'       (default: false)
%       'compareEMG'    (default: false)
%
%OUTPUTS

%%
p = inputParser;
addParameter(p,'channels','all');
addParameter(p,'saveMat',true);  %DOESN"T YET WORK
addParameter(p,'filterbounds',[500 5000]);
addParameter(p,'movingNorm',false);
addParameter(p,'MUAsmoothwin',0.005);
addParameter(p,'usepeaks',false);
addParameter(p,'SHOWFIG',false);
addParameter(p,'compareEMG',false);

parse(p,varargin{:})

channels = p.Results.channels;
saveMat = p.Results.saveMat;
MUAfilter = p.Results.filterbounds;
movingNorm = p.Results.movingNorm;
MUAsmoothwin = p.Results.MUAsmoothwin;
usepeaks = p.Results.usepeaks;
SHOWFIG = p.Results.SHOWFIG;
compareEMG = p.Results.compareEMG;

%% DEV
%basePath = '/mnt/proraidDL/Database/WMProbeData/180211_WT_M1M3_LFP_Layers_Pupil_EMG_Pole/180211_WT_M1M3_LFP_Layers_Pupil_EMG_180211_130605';
%channels = 10;
%note: after 180211 good grounding

%%
sessionInfo = bz_getSessionInfo(basePath, 'noPrompts', true);
datSampleRate = sessionInfo.rates.wideband;
lfpSampleRate = sessionInfo.lfpSampleRate;
downfactor = datSampleRate./lfpSampleRate;

%%
baseName = bz_BasenameFromBasepath(basePath);

%%
datfilename = fullfile(basePath,[baseName,'.dat']);
savefile = fullfile(basePath,[baseName,'.mua.mat']);

%%
if isequal(channels,'all')
    channels = sessionInfo.channels;
end

%%
lfp = bz_GetLFP(1,'basepath',basePath,'noPrompts',true);

%%
downfactor = 1;
MUA = NaN(size(lfp.data,1),length(channels));
for i = 1:length(channels)
    
    datlfp.data = bz_LoadBinary(datfilename,...
        'frequency',datSampleRate,'nchannels',sessionInfo.nChannels,...
        'channels',channels(i)+1,'downsample',downfactor);
    
    datlfp.samplingRate = datSampleRate./downfactor;
    datlfp.timestamps = [0:(length(datlfp.data)-1)]'/datlfp.samplingRate;  %To be overwritten later...
%     datlfp.channels = channels;
    
    %%
    MUALFP = bz_Filter(datlfp,'passband',MUAfilter,'filter','fir1','order',3);
%     MUALFP.channels = channels;
    
    if movingNorm
        %Normalize and re-hilbert
        MUALFP.data = NormToInt(MUALFP.data,'modZ',[0 Inf],MUALFP.samplingRate,'moving',0.2);
        MUALFP.hilb = hilbert(MUALFP.data);
        MUALFP.amp = abs(MUALFP.hilb);
        
        % MUALFP.data_mov1000 = NormToInt(MUALFP.data,'modZ',[0 Inf],MUALFP.samplingRate,'moving',1);
        % MUALFP.hilb = hilbert(MUALFP.data_mov1000);
        % MUALFP.amp_mov1000 = abs(MUALFP.hilb);
    else
        
        %MUALFP.data = NormToInt(MUALFP.data,'modZ',[0 Inf],MUALFP.samplingRate);
        %MUALFP.hilb = hilbert(MUALFP.data);
        %MUALFP.amp = abs(MUALFP.hilb);
    end
    
    %% Smooth the MUA
    MUALFP.smoothamp = smooth(MUALFP.amp,round(MUAsmoothwin.*MUALFP.samplingRate),'moving' );
    % MUALFP.smoothamp_mov200 = smooth(MUALFP.amp_mov200,round(MUAsmoothwin.*MUALFP.samplingRate),'moving' );
    % MUALFP.smoothamp_mov1000 = smooth(MUALFP.amp_mov1000,round(MUAsmoothwin.*MUALFP.samplingRate),'moving' );
    
    %% Threshold crossings
    % figure
    % plot(MUALFP.smoothamp,MUALFP.data,'.')
    
    smMUA = MUALFP.smoothamp;
    
    %% Interpolate to lfp timestamps
    
    MUA(:,i) = interp1(datlfp.timestamps,smMUA,lfp.timestamps);
      
end

%% Saving onto

MUA_lfp.data = MUA;
MUA_lfp.timestamps = lfp.timestamps;
MUA_lfp.channels = channels;

save(savefile,'MUA_lfp');

