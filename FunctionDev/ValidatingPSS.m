basePath = pwd;
baseName = bz_BasenameFromBasepath(basePath);
sessionInfo = bz_getSessionInfo(basePath,'noPrompts',true);

savefile = fullfile(basePath,[baseName,'.PSS.lfp.mat']);
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
load(fullfile(basePath,[baseName,'.SlowWaves.events.mat']),'SlowWaves');
channel = SlowWaves.detectorinfo.detectionchannel; 
clear SlowWaves;

lfp = bz_GetLFP(channel,'basepath',basePath,'noPrompts',true);

% Separate fractal and oscillatory components using sliding window
srate = lfp.samplingRate; % sampling frequency SPECIFY...
movingwin = [0.5 0.125; 2 0.5; 5 1.25].*srate; % [window size, sliding step]
lowerbound = [1:30];
upperbound = [40:120];

CorrFracEMG = [];
CorrFracPupil = [];
for x = 1:size(movingwin,1)
    nwin = floor((length(lfp.data) - movingwin(x,1))/movingwin(x,2));
    sig = zeros(movingwin(x,1),nwin);
    
    for i = 1 : nwin
        sig(:,i) = lfp.data(ceil((i-1)*movingwin(x,2))+1 : ceil((i-1)*movingwin(x,2))+movingwin(x,1));
    end
    
    Frac = amri_sig_fractal_gpu(sig,srate,'detrend',1);
    st = movingwin(x,1)/(srate*2);
    Frac.timestamps = st + ((0:nwin-1) * movingwin(x,2)/srate);
    
    % Interpolating...
    PSS.EMG = interp1(EMGwhisk.timestamps,EMGwhisk.EMGenvelope,Frac.timestamps);
    PSS.pupilsize = interp1(pupildilation.timestamps,pupildilation.data,...
        Frac.timestamps,'nearest');
    
    % Obtaining corr
    FracEMGrho = zeros(length(lowerbound),length(upperbound));
    FracPupilrho = zeros(length(lowerbound),length(upperbound));
    
    for f = 1:length(lowerbound)
        for ff = 1:length(upperbound)
            Frange = [lowerbound(f), upperbound(ff)]; % define frequency range for power-law fitting
            Frac = amri_sig_plawfit(Frac,Frange);
            
            FracEMGrho(f,ff) = corr(log10(PSS.EMG),Frac.Beta.*-1,'Type','Spearman');
            FracPupilrho(f,ff) = corr(log10(PSS.pupilsize),Frac.Beta.*-1,'Type','Spearman');
        end
    end
    
    CorrFracEMG = cat(3,CorrFracEMG,FracEMGrho);
    CorrFracPupil = cat(3,CorrFracPupil,FracPupilrho);
end

%% FIGURE
rwbcolormap = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);

figure;
for x = 1:size(movingwin,1)
    subplot(x,2,(x*2)-1);
    imagesc(lowerbound,upperbound,CorrFracEMG(:,:,x)')
    colormap(gca,rwbcolormap)
    axis xy
    axis tight
    ColorbarWithAxis([min(min(CorrFracEMG(:,:,x))) max(max(CorrFracEMG(:,:,x)))],['Spearman corr'])
    caxis([min(min(CorrFracEMG(:,:,x))) max(max(CorrFracEMG(:,:,x)))])
    xlabel('lower f bound (Hz)');ylabel('upper f bound (Hz)');
    title(['Frac-EMG correlation win =',num2str(movingwin(x,1)),'dt=',num2str(movingwin(x,2))]);
    
    subplot(x,2,x*2);
    imagesc(lowerbound,upperbound,CorrFracPupil(:,:,x)')
    colormap(gca,rwbcolormap)
    axis xy
    axis tight
    ColorbarWithAxis([min(min(CorrFracPupil(:,:,x))) max(max(CorrFracPupil(:,:,x)))],['Spearman corr'])
    caxis([min(min(CorrFracPupil(:,:,x))) max(max(CorrFracPupil(:,:,x)))])
    xlabel('lower f bound (Hz)');ylabel('upper f bound (Hz)');
    title('Frac-Pupil diameter correlation');
end

NiceSave('Frac_EMG_Pupil_lower_upperbound_Optimization',figfolder,baseName);

%% Correlating PSS/Oscillatory component to EMG/Pupil diameter by depth
% Take only good channels

%
PSS = [];
Osci = [];

PSScorr.EMG = zeros(size(sessionInfo.AnatGrps.Channels));
PSScorr.Pup = zeros(size(sessionInfo.AnatGrps.Channels));

Oscicorr.EMG = zeros(size(sessionInfo.AnatGrps.Channels));
PSScorr.Pup = zeros(size(sessionInfo.AnatGrps.Channels));

nbins = 50;
PSSstatsdepth.bins = linspace(-2,0,nbins);
PSSstatsdepth.dist = zeros(length(sessionInfo.AnatGrps.Channels),nbins);

for cc = 1:length(sessionInfo.AnatGrps.Channels)
    cc
    channum = sessionInfo.AnatGrps.Channels(cc);
    %channum = 31;
    WhiskPSScorr.channum(cc) = channum;
    WhiskPSScorr.chanpos(cc) = cc;
    
    PSSstatsdepth.channum(cc) = channum;
    PSSstatsdepth.chanpos(cc) = cc;
    %%
    lfp = bz_GetLFP(channum,'basepath',basePath,'noPrompts',true);
    
    %%
    dt = 0.2;
    winsize = 1;
    [PSS] = bz_PowerSpectrumSlope(lfp,winsize,dt,'showfig',false);
    
    PSSstatsdepth.dist(cc,:) = hist(PSS.data,PSSstatsdepth.bins);
    PSSstatsdepth.dist(cc,:)./sum(PSSstatsdepth.dist(cc,:));
    
    %%
    PSS.EMG = interp1(EMGwhisk.timestamps,EMGwhisk.EMGenvelope,PSS.timestamps);
    PSS.pupilsize = interp1(pupildilation.timestamps,pupildilation.data,...
        PSS.timestamps,'nearest');
    PSS.dpdt = interp1(pupildilation.timestamps(1:end-1),pupildilation.dpdt,...
        PSS.timestamps,'nearest');
    PSS.pupilphase = interp1(lowpupildata.timestamps,lowpupildata.phase,...
        PSS.timestamps,'nearest');
    
    %%
    [WhiskPSScorr.EMG(cc),WhiskPSScorr.EMG_p(cc)] =...
        corr(log10(PSS.EMG),PSS.data,...
        'type','spearman','rows','complete');
    [WhiskPSScorr.pup(cc),WhiskPSScorr.pup_p(cc)] =...
        corr(log10(PSS.pupilsize),PSS.data,...
        'type','spearman','rows','complete');
    [WhiskPSScorr.dpdt(cc),WhiskPSScorr.dpdt_p(cc)] =...
        corr(PSS.dpdt,PSS.data,...
        'type','spearman','rows','complete');
    WhiskPSScorr.phasecoupling(cc) = ...
        abs(nanmean((PSS.data./nanmean(PSS.data)).*exp(1i.*PSS.pupilphase)));
    
    clear lfp
    
end

%
% Saving structs

%% FIGURE
figure;
subplot(1,3,1);
imagesc(PSSstatsdepth.bins,PSSstatsdepth.chanpos,PSSstatsdepth.dist)
ColorbarWithAxis([min(min(CorrFracEMG(:,:,x))) max(max(CorrFracEMG(:,:,x)))],['Spearman corr'])
xlabel('PSS')
ylabel('Channel by Depth')
colormap(gca,'jet')
axis tight

subplot(1,3,2);
plot(WhiskPSScorr.EMG,-WhiskPSScorr.chanpos,'b','linewidth',2)
hold on
plot(WhiskPSScorr.pup,-WhiskPSScorr.chanpos,'k','linewidth',2)
legend('EMG','Pupil Area','dpdt','location','eastoutside')
xlabel('PSS Correlation');
axis tight
%xlim([0 0.6])

subplot(1,3,3);
imagesc(lowerbound,upperbound,CorrFracPupil(:,:,x)')
colormap(gca,rwbcolormap)
axis tight
ColorbarWithAxis([min(min(CorrFracEMG(:,:,x))) max(max(CorrFracEMG(:,:,x)))],['Spearman corr'])
caxis([min(min(CorrFracEMG(:,:,x))) max(max(CorrFracEMG(:,:,x)))])
xlabel('lower f bound (Hz)');ylabel('upper f bound (Hz)');
title(['Frac-EMG correlation win =',num2str(movingwin(x,1)),'dt=',num2str(movingwin(x,2))]);

NiceSave('PSSCorrbyDepth',figfolder,baseName)

%% Xcorr PSS/Oscillatory component by depth
PSSxcorr = zeros(size(PSS,1),size(PSS,1));
for x = 1:size(PSS,1)
    for y = 1:size(PSS,1)
        PSSxcorr(x,y) = corr(PSS(x,:),PSS(y,:),'Type','Spearman');
    end
end

Oscixcorr = zeros(size(PSS,1),size(PSS,1));
for x = 1:size(Osci,1)
    for y = 1:size(Osci,1)
        Oscixcorr(x,y) = corr(Osci(x,:),Osci(y,:),'Type','Spearman');
    end
end

%% FIGURE
figure;
subplot(1,2,1);
imagesc(lowerbound,upperbound,CorrFracEMG(:,:,x)')
colormap(gca,rwbcolormap)
axis xy
axis tight
ColorbarWithAxis([min(min(CorrFracEMG(:,:,x))) max(max(CorrFracEMG(:,:,x)))],['Spearman corr'])
caxis([min(min(CorrFracEMG(:,:,x))) max(max(CorrFracEMG(:,:,x)))])
xlabel('lower f bound (Hz)');ylabel('upper f bound (Hz)');
title(['Frac-EMG correlation win =',num2str(movingwin(x,1)),'dt=',num2str(movingwin(x,2))]);

subplot(1,2,2);
imagesc(lowerbound,upperbound,CorrFracPupil(:,:,x)')
colormap(gca,rwbcolormap)
axis xy
axis tight
ColorbarWithAxis([min(min(CorrFracPupil(:,:,x))) max(max(CorrFracPupil(:,:,x)))],['Spearman corr'])
caxis([min(min(CorrFracPupil(:,:,x))) max(max(CorrFracPupil(:,:,x)))])
xlabel('lower f bound (Hz)');ylabel('upper f bound (Hz)');
title('Frac-Pupil diameter correlation');

NiceSave('Frac_Osci_Xcorr',figfolder,baseName);
