function [ output_args ] = LFPSlopeAndPupilAnalysis( basePath,figfolder )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
basePath = '/mnt/proraidDL/Database/WMProbeData/180213_WT_M1M3_LFP_Layers_Pupil_EMG_Pole/180213_WT_M1M3_LFP_Layers_Pupil_EMG_180213_113045';
%basePath = pwd;
%figfolder = '/mnt/data1/Dropbox/research/Current Projects/S1State/AnalysisScripts/figures/UPDOWNandPupilAnalysis';
figfolder = '/home/dlevenstein/ProjectRepos/ACh-and-CorticalState/AnalysisScripts/AnalysisFigs/LFPSlopeAndPupilAnalysis';

%%
sessionInfo = bz_getSessionInfo(basePath,'noPrompts',true);

%%
pupildilation = bz_LoadBehavior(basePath,'pupildiameter');

%pupildilation.dpdt = 
smoothwin =2;%s
pupildilation.dpdt = diff(smooth(pupildilation.data,smoothwin.*pupildilation.samplingRate,'moving')).*pupildilation.samplingRate;
pupildilation.dpdt = smooth(pupildilation.dpdt,smoothwin.*pupildilation.samplingRate,'moving');

nantimes = isnan(pupildilation.data);
pupildilation.interpdata = interp1(pupildilation.timestamps(~nantimes),...
    pupildilation.data(~nantimes),pupildilation.timestamps);
%%
slopepupilcorr.p = zeros(size(sessionInfo.SpkGrps.Channels));
slopepupilcorr.dpdt = zeros(size(sessionInfo.SpkGrps.Channels));
for cc = 1:length(sessionInfo.SpkGrps.Channels)
    cc
    channum = sessionInfo.SpkGrps.Channels(cc);
    %channum = 31;
    %%
    lfp = bz_GetLFP(channum,'basepath',basePath,'noPrompts',true);

    %%
    dt = 0.5;
    winsize = 2;
    [specslope,spec] = bz_PowerSpectrumSlope(lfp,winsize,dt,'showfig',false);

    %%
    specslope.pupilsize = interp1(pupildilation.timestamps,pupildilation.data,...
        specslope.timestamps,'nearest');
    specslope.dpdt = interp1(pupildilation.timestamps(1:end-1),pupildilation.dpdt,...
        specslope.timestamps,'nearest');

    %%
    slopepupilcorr.p(cc) = corr(log10(specslope.pupilsize),specslope.data,'type','spearman');
    slopepupilcorr.dpdt(cc) = corr(specslope.dpdt,specslope.data,'type','spearman');
end
%%
figure
subplot(2,2,1)
    plot(log10(specslope.pupilsize),specslope.data,'.')
    xlabel('Pupil Area (med^-^1)');ylabel('Spectrum Slope')
subplot(2,2,2)
    plot((specslope.dpdt),specslope.data,'.')
    xlabel('dpdt (med^-^1s^-^1)');ylabel('Spectrum Slope')
end

