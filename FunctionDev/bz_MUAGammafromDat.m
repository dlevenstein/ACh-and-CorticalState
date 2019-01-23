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

%% DEV
basePath = '/mnt/proraidDL/Database/WMProbeData/180211_WT_M1M3_LFP_Layers_Pupil_EMG_Pole/180211_WT_M1M3_LFP_Layers_Pupil_EMG_180211_130605';
channels = 10;
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
%%

if isequal(channels,'all')
    channels = sessionInfo.channels;
end

%%
% Load data and put into struct
% we assume 0-indexing like neuroscope, but bz_LoadBinary uses 1-indexing to
% load....
downfactor = 1;
datlfp.data = bz_LoadBinary(datfilename,...
              'frequency',datSampleRate,'nchannels',sessionInfo.nChannels,...
              'channels',channels+1,'downsample',downfactor);

datlfp.samplingRate = datSampleRate./downfactor;
datlfp.timestamps = [0:(length(datlfp.data)-1)]'/datlfp.samplingRate;  %To be overwritten later...
datlfp.channels = channels;

%%
%oldlfp = bz_GetLFP(channels,'basepath',basePath);

%%
gammafilter = [100 400];
gammaLFP = bz_Filter(datlfp,'passband',gammafilter,'filter','fir1','order',4);

%%
MUAfilter = [500 5000];
MUALFP = bz_Filter(datlfp,'passband',MUAfilter,'filter','fir1','order',3);
MUALFP.channels = channels;

%Normalize and re-hilbert
MUALFP.data_mov200 = NormToInt(MUALFP.data,'modZ',[0 Inf],MUALFP.samplingRate,'moving',0.2);
MUALFP.hilb = hilbert(MUALFP.data_mov200);
MUALFP.amp_mov200 = abs(MUALFP.hilb);

MUALFP.data_mov1000 = NormToInt(MUALFP.data,'modZ',[0 Inf],MUALFP.samplingRate,'moving',1);
MUALFP.hilb = hilbert(MUALFP.data_mov1000);
MUALFP.amp_mov1000 = abs(MUALFP.hilb);

MUALFP.data = NormToInt(MUALFP.data,'modZ',[0 Inf],MUALFP.samplingRate,'moving',0);
MUALFP.hilb = hilbert(MUALFP.data);
MUALFP.amp = abs(MUALFP.hilb);
%% Smooth the MUA
MUAsmoothwin = 0.005; %s
MUALFP.smoothamp = smooth(MUALFP.amp,round(MUAsmoothwin.*MUALFP.samplingRate),'moving' );
MUALFP.smoothamp_mov200 = smooth(MUALFP.amp_mov200,round(MUAsmoothwin.*MUALFP.samplingRate),'moving' );
MUALFP.smoothamp_mov1000 = smooth(MUALFP.amp_mov1000,round(MUAsmoothwin.*MUALFP.samplingRate),'moving' );

%% Threshold crossings
% figure
% plot(MUALFP.smoothamp,MUALFP.data,'.')


%%
thresh = 3;
[MUALFP.peakmags,MUALFP.peaktimes] = findpeaks(-MUALFP.data,'MinPeakHeight',thresh);
MUALFP.peaktimes = MUALFP.timestamps(MUALFP.peaktimes);

[MUALFP.peakmags_mov200,MUALFP.peaktimes_mov200] = findpeaks(-MUALFP.data_mov200,'MinPeakHeight',thresh);
MUALFP.peaktimes_mov200 = MUALFP.timestamps(MUALFP.peaktimes_mov200);
%%
figure
hist(MUALFP.peakmags)

%%
MUAspks = bz_SpktToSpkmat({MUALFP.peaktimes},'binsize',0.02);
MUAspks_mov200 = bz_SpktToSpkmat({MUALFP.peaktimes_mov200},'binsize',0.02);

%%
%%
% figure
% subplot(2,2,1)
% hist(MUALFP.data,25)
% subplot(2,2,2)
% hist(MUALFP.smoothamp,25)
%%
showwin = bz_RandomWindowInIntervals(MUALFP.timestamps([1 end]),0.5);

timewin = MUALFP.timestamps >showwin(1) & MUALFP.timestamps <showwin(2);

figure
subplot(3,1,2)
plot(MUALFP.timestamps(timewin),MUALFP.data(timewin),'k')
hold on
plot(MUALFP.timestamps(timewin),MUALFP.amp(timewin),'b')
plot(MUALFP.timestamps(timewin),MUALFP.smoothamp(timewin),'r')
plot(MUALFP.peaktimes,-MUALFP.peakmags,'.')
plot(MUAspks.timestamps,MUAspks.data,'g')
legend(['MUA (',num2str(MUAfilter(1)),'-',num2str(MUAfilter(2)),'Hz)'],...
    'Amplitude',['Smooth: ',num2str(MUAsmoothwin.*1000),'ms'])
xlim(showwin)


subplot(3,1,3)
plot(MUALFP.timestamps(timewin),MUALFP.data_mov200(timewin),'k')
hold on
plot(MUALFP.timestamps(timewin),MUALFP.amp_mov200(timewin),'b')
plot(MUALFP.timestamps(timewin),MUALFP.smoothamp_mov200(timewin),'r')
plot(MUALFP.peaktimes_mov200,-MUALFP.peakmags_mov200,'.')
plot(MUAspks_mov200.timestamps,MUAspks_mov200.data,'g')
legend(['MUA (',num2str(MUAfilter(1)),'-',num2str(MUAfilter(2)),'Hz)'],...
    'Amplitude',['Smooth: ',num2str(MUAsmoothwin.*1000),'ms'])
xlim(showwin)


subplot(3,1,1)
plot(datlfp.timestamps(timewin),datlfp.data(timewin),'k')
hold on
bz_ScaleBar('s')
xlim(showwin)

%% Get the EMG
[EMGFromLFP] = bz_EMGFromLFP(basePath,'samplingFrequency',10,'fromDat',true);

%%
EMGFromLFP.MUA = interp1(MUALFP.timestamps,MUALFP.smoothamp,EMGFromLFP.timestamps);
smoothcorr = corr(EMGFromLFP.MUA,EMGFromLFP.data)
EMGFromLFP.MUA_smooth200 = interp1(MUALFP.timestamps,MUALFP.smoothamp_mov200,EMGFromLFP.timestamps);
smoothcorr_smooth200 = corr(EMGFromLFP.MUA_smooth200,EMGFromLFP.data)
EMGFromLFP.MUA_smooth1000 = interp1(MUALFP.timestamps,MUALFP.smoothamp_mov1000,EMGFromLFP.timestamps);
smoothcorr_smooth1000 = corr(EMGFromLFP.MUA_smooth1000,EMGFromLFP.data)
EMGFromLFP.MUAspks = interp1(MUAspks.timestamps,MUAspks.data,EMGFromLFP.timestamps);
bincorr = corr(EMGFromLFP.MUAspks,EMGFromLFP.data)
%%
showwin = bz_RandomWindowInIntervals(MUALFP.timestamps([1 end]),10);
timewin = MUALFP.timestamps >showwin(1) & MUALFP.timestamps <showwin(2);

figure
subplot(4,1,1)
plot(EMGFromLFP.timestamps,EMGFromLFP.data)
xlim(showwin)
ylabel('EMGfromLFP')

subplot(4,1,2)
plot(MUALFP.timestamps(timewin),MUALFP.smoothamp(timewin),'r')
xlim(showwin)
ylabel('MUA')

subplot(4,1,3)
plot(MUALFP.timestamps(timewin),MUALFP.smoothamp_mov200(timewin),'r')
xlim(showwin)
ylabel('MUA')

subplot(4,1,4)
plot(MUALFP.timestamps(timewin),MUALFP.smoothamp_mov1000(timewin),'r')
xlim(showwin)
ylabel('MUA')

% subplot(3,1,3)
% plot(MUAspks.timestamps,MUAspks.data,'r')
% xlim(showwin)
% ylabel('MUA Binned')
%%
%TO DO: vary lower and upper - quantify correlation with EMGfromLFP
%(calculate EMGfromLFP)

%%
%Threshold crossings?  Time window?
%
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

