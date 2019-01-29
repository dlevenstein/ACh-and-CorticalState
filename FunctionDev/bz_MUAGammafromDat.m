function [ lof, lospec, hif, hispec, MUA ] = bz_MUAGammafromDat( basePath,varargin )
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
%       'filterbounds'  (default: [500 5000])
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
              'channels',channels,'downsample',downfactor);

datlfp.samplingRate = datSampleRate./downfactor;
datlfp.timestamps = [0:(length(datlfp.data)-1)]'/datlfp.samplingRate;  %To be overwritten later...
datlfp.channels = channels;

%%
%oldlfp = bz_GetLFP(channels,'basepath',basePath);

%% 
lfp = bz_GetLFP(1,'basepath',basePath,'noPrompts',true);

%%
freqlist = logspace(log10(0.5),log10(100),100);
window = 1;
noverlap = 0.8;
window = window*datlfp.samplingRate;
noverlap = noverlap*datlfp.samplingRate;

[spec,lof,t_FFT] = spectrogram(single(datlfp.data),window,noverlap,freqlist,datlfp.samplingRate);
lospec = log10(abs(spec));
lospec = mean(lospec,2);

freqlist = logspace(log10(100),log10(10000),1000);
window = 10;
noverlap = 5;
window = window*datlfp.samplingRate;
noverlap = noverlap*datlfp.samplingRate;

[spec,hif,t_FFT] = spectrogram(single(datlfp.data),window,noverlap,freqlist,datlfp.samplingRate);
hispec = log10(abs(spec));
hispec = mean(hispec,2);

%%
% viewwin = bz_RandomWindowInIntervals(lfp.timestamps([1 end]),500);
% figure
% subplot(4,1,1)
% imagesc(T,log10(F),amp)
% axis xy
% LogScale('y',10)
% xlim(viewwin)
% subplot(4,1,2)
% plot(lfp.timestamps,lfp.data,'k')
% xlim(viewwin)
% 
% subplot(2,2,3)
% plot(log10(freqs),PSD,'k')
% LogScale('x',10)

%%
%pxx = pspectrum(single(datlfp.data),datlfp.samplingRate,'FrequencyLimits',[0.5 100],'FrequencyResolution',1);
%[pxx,f] = pspectrum(single(datlfp.data),datlfp.samplingRate,'FrequencyLimits',[100 10000],'FrequencyResolution',25);
%[pxx,f] = pspectrum(single(datlfp.data),datlfp.samplingRate,'FrequencyLimits',[0.5 10000],'FrequencyResolution',25);
%pspec = pow2db(pxx);

%%
% gammafilter = [100 400];
% gammaLFP = bz_Filter(datlfp,'passband',gammafilter,'filter','fir1','order',4);

%%
MUALFP = bz_Filter(datlfp,'passband',MUAfilter,'filter','fir1','order',3);
MUALFP.channels = channels;

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

%% Interpolate to lfp timestamps

MUA.data = interp1(datlfp.timestamps,MUALFP.smoothamp,lfp.timestamps);
MUA.timestamps = lfp.timestamps;

%%
if usepeaks
    thresh = 3;
    [MUALFP.peakmags,MUALFP.peaktimes] = findpeaks(-MUALFP.data,'MinPeakHeight',thresh);
    MUALFP.peaktimes = MUALFP.timestamps(MUALFP.peaktimes);

    MUAspks = bz_SpktToSpkmat({MUALFP.peaktimes},'binsize',0.02);

end

%%
% figure
% hist(MUALFP.peakmags)

%%
if SHOWFIG
    showwin = bz_RandomWindowInIntervals(MUALFP.timestamps([1 end]),0.5);

    timewin = MUALFP.timestamps >showwin(1) & MUALFP.timestamps <showwin(2);

    figure
    subplot(3,1,2)
    plot(MUALFP.timestamps(timewin),MUALFP.data(timewin),'k')
    hold on
    plot(MUALFP.timestamps(timewin),MUALFP.amp(timewin),'b')
    plot(MUALFP.timestamps(timewin),MUALFP.smoothamp(timewin),'r')
   % plot(MUALFP.peaktimes,-MUALFP.peakmags,'.')
    %plot(MUAspks.timestamps,MUAspks.data,'g')
    legend(['MUA (',num2str(MUAfilter(1)),'-',num2str(MUAfilter(2)),'Hz)'],...
        'Amplitude',['Smooth: ',num2str(MUAsmoothwin.*1000),'ms'])
    xlim(showwin)


    subplot(3,1,1)
    plot(datlfp.timestamps(timewin),datlfp.data(timewin),'k')
    hold on
    bz_ScaleBar('s')
    xlim(showwin)
end

%% Get the EMG
if compareEMG
    [EMGFromLFP] = bz_EMGFromLFP(basePath,'samplingFrequency',10,'fromDat',true);

    %%
    EMGFromLFP.MUA = interp1(MUALFP.timestamps,MUALFP.smoothamp,EMGFromLFP.timestamps);
    smoothcorr = corr(EMGFromLFP.MUA,EMGFromLFP.data)

    if SHOWFIG
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

    end 
end

%%
%TO DO: vary lower and upper - quantify correlation with EMGfromLFP
%(calculate EMGfromLFP)

%%
%Threshold crossings?  Time window?
%
% %%
% gammasmoothwin = 0.08; %s
% gammaLFP.smoothamp = smooth(gammaLFP.amp,round(gammasmoothwin.*gammaLFP.samplingRate),'moving' );
% %%
% showwin = bz_RandomWindowInIntervals(trygammaLFP.timestamps([1 end]),1);
% figure
% subplot(4,1,1)
% plot(datlfp.timestamps,datlfp.data(:,1),'k')
% hold on
% plot(oldlfp.timestamps,oldlfp.data(:,1),'r')
% xlim(showwin)
% subplot(4,1,2)
% plot(gammaLFP.timestamps,gammaLFP.amp(:,1),'k')
% hold on
% plot(gammaLFP_fromlfp.timestamps,gammaLFP_fromlfp.amp(:,1),'r')
% plot(gammaLFP_cheb.timestamps,gammaLFP_cheb.amp(:,1),'k--')
% plot(gammaLFP_fromlfpcheb.timestamps,gammaLFP_fromlfpcheb.amp(:,1),'r--')
% plot(gammaLFP.timestamps,gammaLFP.smoothamp(:,1),'k','linewidth',2)
% xlim(showwin)

end

