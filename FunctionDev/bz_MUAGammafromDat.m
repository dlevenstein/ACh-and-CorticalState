function [ MUAGamma ] = bz_MUAGammafromDat( basePath,varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%INPUTS
%   basePath
%   
%   (options)
%       'channels'  (default: all)   (0-indexed, like neuroscope)
%       'saveMat'   (default: true) save a .mat file
%                       basePath/baseName.MUAGamma.lfp.mat
%
%OUTPUTS

%%
p = inputParser;
addParameter(p,'channels','all');
addParameter(p,'saveMat',true);
parse(p,varargin{:})

channels = p.Results.channels;
saveMat = p.Results.saveMat;

%%
sessionInfo = bz_getSessionInfo(basePath, 'noPrompts', true);
datSampleRate = sessionInfo.rates.wideband;
lfpSampleRate = sessionInfo.lfpSampleRate;
downfactor = datSampleRate./lfpSampleRate;

%%
baseName = bz_BasenameFromBasepath(basePath);
%%
datfilename = fullfile(basePath,[baseName,'.dat']);
%%

if isequal(channels,'all')
    channels = 0:sessionInfo.nChannels-1;
end

%%
% Load data and put into struct
% we assume 0-indexing like neuroscope, but bz_LoadBinary uses 1-indexing to
% load....
datlfp.data = bz_LoadBinary(datfilename,...
              'frequency',datSampleRate,'nchannels',sessionInfo.nChannels,...
              'channels',channels+1,'downsample',downfactor);
datlfp.timestamps = [0:(length(datlfp.data)-1)]'/lfpSampleRate;
datlfp.channels = channels;
datlfp.samplingRate = lfpSampleRate;

%%
oldlfp = bz_GetLFP(channels,'basepath',basePath);

%%
gammafilter = [100 400];
gammaLFP = bz_Filter(datlfp,'passband',gammafilter,'filter','fir1','order',4);

%%
gammasmoothwin = 0.08; %s
gammaLFP.smoothamp = smooth(gammaLFP.amp,round(gammasmoothwin.*gammaLFP.samplingRate),'moving' );
%%
showwin = bz_RandomWindowInIntervals(trygammaLFP.timestamps([1 end]),1);
figure
subplot(4,1,1)
plot(datlfp.timestamps,datlfp.data(:,1),'k')
hold on
plot(oldlfp.timestamps,oldlfp.data(:,1),'r')
xlim(showwin)
subplot(4,1,2)
plot(gammaLFP.timestamps,gammaLFP.amp(:,1),'k')
hold on
plot(gammaLFP_fromlfp.timestamps,gammaLFP_fromlfp.amp(:,1),'r')
plot(gammaLFP_cheb.timestamps,gammaLFP_cheb.amp(:,1),'k--')
plot(gammaLFP_fromlfpcheb.timestamps,gammaLFP_fromlfpcheb.amp(:,1),'r--')
plot(gammaLFP.timestamps,gammaLFP.smoothamp(:,1),'k','linewidth',2)
xlim(showwin)
end

