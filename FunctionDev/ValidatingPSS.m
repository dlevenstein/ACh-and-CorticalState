basePath = pwd;
baseName = bz_BasenameFromBasepath(basePath);
sessionInfo = bz_getSessionInfo(basePath,'noPrompts',true);

figfolder = fullfile(basePath,'DetectionFigures');

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
pupildilation.timestamps = pupildilation.timestamps(~nantimes);

% EMG
EMGwhisk = bz_LoadBehavior(basePath,'EMGwhisk');

%Specifying SPONT whisking
load(fullfile(basePath,[baseName,'.MergePoints.events.mat']),'MergePoints');
sidx = find(startsWith(MergePoints.foldernames,"Spont"));
sponttimes = [MergePoints.timestamps(sidx(1),1) MergePoints.timestamps(sidx(end),2)];

spontidx = find(EMGwhisk.ints.Wh(:,2) < sponttimes(2));
EMGwhisk.ints.Wh = EMGwhisk.ints.Wh(spontidx,:);

spontidx = find(EMGwhisk.ints.NWh(:,2) < sponttimes(2));
EMGwhisk.ints.NWh = EMGwhisk.ints.NWh(spontidx,:);

spontidx = find(EMGwhisk.timestamps < sponttimes(2));
EMGwhisk.timestamps = EMGwhisk.timestamps(spontidx);
EMGwhisk.EMGenvelope = EMGwhisk.EMGenvelope(spontidx);
EMGwhisk.EMG = EMGwhisk.EMG(spontidx);
EMGwhisk.EMGsm = EMGwhisk.EMGsm(spontidx);

%% Optimizing lower/upper frequency bounds as a function of window/dt
% Loading LFP for SlowWave detection channel
%channel = 2; specify slow detection channel
lfp = bz_GetLFP(channel,'basepath',basePath,'noPrompts',true);


% Separate fractal and oscillatory components using sliding window
%srate = 1250; % sampling frequency SPECIFY...
movingwin = [0.5 0.125; 2 0.5; 5 1.25].*srate; % [window size, sliding step]
lowerbound = [1:30];
upperbound = [40:120];

Frac = [];
Osci = [];
Frac_times = [];
FracEMGcorr = zeros(length(lowerbound),length(upperbound));
FracPupilcorr = zeros(length(lowerbound),length(upperbound));
fracttot = [];
for x = 1:size(movingwin,1)
    nwin = floor((length(lfp.data) - win)/step);
    sig = zeros(win,nwin);
    for i = 1 : nwin
        sig(:,i) = lfp.data(ceil((i-1)*step)+1 : ceil((i-1)*step)+win);
    end
    
    Frac = amri_sig_fractal(sig,srate,'detrend',1);
    
    Frac.time = (1:step/srate:step*(nwin+(movingwin(1)/2))/srate)';
    
    % Interpolating...
PSS.EMG = interp1(EMGwhisk.timestamps,EMGwhisk.EMGenvelope,PSS.timestamps);
PSS.pupilsize = interp1(pupildilation.timestamps,pupildilation.data,...
    PSS.timestamps,'nearest');
PSS.dpdt = interp1(pupildilation.timestamps(1:end-1),pupildilation.dpdt,...
    PSS.timestamps,'nearest');
PSS.pupilphase = interp1(lowpupildata.timestamps,lowpupildata.phase,...
    PSS.timestamps,'nearest');

for f = 1:length(lowerbound)
    for ff = 1:length(upperbound)
        Frange = [lowerbound(f), upperbound(ff)]; % define frequency range for power-law fitting
        Frac_short = amri_sig_plawfit(Frac_short,Frange);
        
        FracEMGcorr(f,ff) = corr(Frac_short.Beta.*-1,PSS.EMG(100:150));
        FracPupilcorr(f,ff) = corr(Frac_short.Beta.*-1,PSS.pupilsize(100:150));
        FracdpPdtcorr(f,ff) = corr(Frac_short.Beta.*-1,PSS.dpdt(100:150));
        FracPphasecorr(f,ff) = corr(Frac_short.Beta.*-1,PSS.pupilphase(100:150));
    end
end

end

%% FIGURE
rwbcolormap = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);

figure;

subplot(2,2,1);
imagesc(lowerbound,upperbound,FracEMGcorr')
colormap(gca,rwbcolormap)
axis xy
axis tight
ColorbarWithAxis([min(min(FracEMGcorr)) max(max(FracEMGcorr))],['Spearman corr'])
caxis([min(min(FracEMGcorr)) max(max(FracEMGcorr))])
xlabel('lower f bound (Hz)');ylabel('upper f bound (Hz)');
title('Frac-EMG corr');

subplot(2,2,2);
imagesc(lowerbound,upperbound,FracPupilcorr')
colormap(gca,rwbcolormap)
axis xy
axis tight
ColorbarWithAxis([min(min(FracPupilcorr)) max(max(FracPupilcorr))],['Spearman corr'])
caxis([min(min(FracPupilcorr)) max(max(FracPupilcorr))])
xlabel('lower f bound (Hz)');ylabel('upper f bound (Hz)');
title('Frac-Pupil diameter corr');

subplot(2,2,3);
imagesc(lowerbound,upperbound,FracdpPdtcorr')
colormap(gca,rwbcolormap)
axis xy
axis tight
ColorbarWithAxis([min(min(FracdpPdtcorr)) max(max(FracdpPdtcorr))],['Spearman corr'])
caxis([min(min(FracdpPdtcorr)) max(max(FracdpPdtcorr))])
xlabel('lower f bound (Hz)');ylabel('upper f bound (Hz)');
title('Frac-dpPdt corr');

subplot(2,2,4);
imagesc(lowerbound,upperbound,FracPphasecorr')
colormap(gca,rwbcolormap)
axis xy
axis tight
ColorbarWithAxis([min(min(FracPphasecorr)) max(max(FracPphasecorr))],['Spearman corr'])
caxis([min(min(FracPphasecorr)) max(max(FracPphasecorr))])
xlabel('lower f bound (Hz)');ylabel('upper f bound (Hz)');
title('Frac-Pupil phase corr');

NiceSave('Frac_EMG_Pupil_lower_upperbound_Optimization_win05_dt0125',figfolder,baseName);

%% Correlating PSS/Oscillatory component to EMG/Pupil diameter by depth

%% FIGURE

%% Xcorr PSS/Oscillatory component by depth  

%% FIGURE


