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
load(fullfile(basePath,[baseName,'.SlowWaves.events.mat']),'SlowWaves');
channel = SlowWaves.detectorinfo.detectionchannel; 
clear SlowWaves;

lfp = bz_GetLFP(channel,'basepath',basePath,'noPrompts',true);

% Separate fractal and oscillatory components using sliding window
srate = lfp.samplingRate; % sampling frequency SPECIFY...
movingwin = [0.5 0.125; 2 0.5; 5 1.25].*srate; % [window size, sliding step]
lowerbound = logspace(log10(1),log10(30),20);
upperbound = logspace(log10(40),log10(120),20);

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
           
            FracEMGrho(f,ff) = corr(log10(PSS.EMG)',...
                Frac.Beta.*-1,'type','spearman','rows','complete');
            
            FracPupilrho(f,ff) = corr(log10(PSS.pupilsize)',...
                Frac.Beta.*-1,'type','spearman','rows','complete');
        end
    end
    
    CorrFracEMG = cat(3,CorrFracEMG,FracEMGrho);
    CorrFracPupil = cat(3,CorrFracPupil,FracPupilrho);
end

% Saving data...
savefile = fullfile(basePath,[baseName,'.FracCorr.Optimization.PSS.lfp.mat']);
save(savefile,'CorrFracEMG','CorrFracPupil');

%% FIGURE
rwbcolormap = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);

figure;
for x = 1:size(movingwin,1)
    subplot(x,2,(x*2)-1);
    imagesc(lowerbound,upperbound,CorrFracEMG(:,:,x)')
    colormap(gca,rwbcolormap)
    LogScale('x',10);
    axis tight
    ColorbarWithAxis([min(min(CorrFracEMG(:,:,x))) max(max(CorrFracEMG(:,:,x)))],['Spearman corr'])
    caxis([min(min(CorrFracEMG(:,:,x))) max(max(CorrFracEMG(:,:,x)))])
    xlabel('lower f bound (Hz)');ylabel('upper f bound (Hz)');
    title(['Frac-EMG correlation win =',num2str(movingwin(x,1)),'dt=',num2str(movingwin(x,2))]);
    
    subplot(x,2,x*2);
    imagesc(lowerbound,upperbound,CorrFracPupil(:,:,x)')
    colormap(gca,rwbcolormap)
    LogScale('x',10);
    axis tight
    ColorbarWithAxis([min(min(CorrFracPupil(:,:,x))) max(max(CorrFracPupil(:,:,x)))],['Spearman corr'])
    caxis([min(min(CorrFracPupil(:,:,x))) max(max(CorrFracPupil(:,:,x)))])
    xlabel('lower f bound (Hz)'); ylabel('upper f bound (Hz)');
    title('Frac-Pupil diameter correlation');
end

NiceSave('Frac_EMG_Pupil_lower_upperbound_Optimization',figfolder,baseName);

%% Correlating PSS/Oscillatory component to EMG/Pupil diameter by depth
badchannels = sessionInfo.badchannels;
usechannels = sessionInfo.AnatGrps.Channels;
usechannels(ismember(usechannels,badchannels))=[];

% Assuming that LFP still remains loaded from prior analysis...
movingwin = [2 0.5].*srate;
nwin = floor((length(lfp.data) - movingwin(1))/movingwin(2));
st = movingwin(1)/(srate*2);

% Selected from prior analysis
Frange = [10, 100]; % define frequency range for power-law fitting
            
% Now calculating corrs...
PSS = []; Osci = [];

PSScorr.EMG = zeros(length(usechannels));
PSScorr.EMG_p = zeros(length(usechannels));
PSScorr.Pup = zeros(length(usechannels));
PSScorr.Pup_p = zeros(length(usechannels));

Oscicorr.EMG = zeros(size(Frac.osci,1),length(usechannels));
Oscicorr.EMG_p = zeros(size(Frac.osci,1),length(usechannels));
Oscicorr.Pup = zeros(size(Frac.osci,1),length(usechannels));
Oscicorr.Pup_p = zeros(size(Frac.osci,1),length(usechannels));

nbins = 50;
PSSstatsdepth.bins = linspace(-2,0,nbins);
PSSstatsdepth.dist = zeros(length(usechannels),nbins);

for cc = 1:length(usechannels)
    cc
    channum = usechannels(cc);
    lfp = bz_GetLFP(channum,'basepath',basePath,'noPrompts',true);
    
    %% Deconstruction
    sig = zeros(movingwin(1),nwin);
    for i = 1 : nwin
        sig(:,i) = lfp.data(ceil((i-1)*movingwin(2))+1 : ceil((i-1)*movingwin(2))+movingwin(1));
    end
    clear lfp
    
    Frac = amri_sig_fractal_gpu(sig,srate,'detrend',1);
    Frac.timestamps = st + ((0:nwin-1) * movingwin(2)/srate);
    Frac = amri_sig_plawfit(Frac,Frange);
    
    PSS = cat(2,PSS,Frac.Beta.*-1);
    Osci = cat(2,Osci,mean(Frac.osci,2));
    
    %% Histos
    PSSstatsdepth.dist(cc,:) = hist(Frac.Beta.*-1,PSSstatsdepth.bins);
    PSSstatsdepth.dist(cc,:)./sum(PSSstatsdepth.dist(cc,:));
    
    %% Interpolating...
    temp_EMG = interp1(EMGwhisk.timestamps,EMGwhisk.EMGenvelope,Frac.timestamps);
    temp_Pup = interp1(pupildilation.timestamps,pupildilation.data,...
        Frac.timestamps,'nearest');
    
    %% Corrsss
    [PSScorr.EMG(cc),PSScorr.EMG_p(cc)] =...
        corr(log10(temp_EMG)',Frac.Beta.*-1,...
        'type','spearman','rows','complete');
    [PSScorr.Pup(cc),PSScorr.Pup_p(cc)] =...
        corr(log10(temp_Pup)',Frac.Beta.*-1,...
        'type','spearman','rows','complete');
    
    for x = 1:size(Frac.osci,1)
        [Oscicorr.EMG(x,cc),Oscicorr.EMG_p(x,cc)] =...
            corr(log10(temp_EMG)',Frac.osci(x,:),...
            'type','spearman','rows','complete');
        [Oscicorr.Pup(x,cc),Oscicorr.Pup_p(x,cc)] =...
            corr(log10(temp_Pup)',Frac.osci(x,:),...
            'type','spearman','rows','complete');
    end
end

% Saving data...
savefile = fullfile(basePath,[baseName,'.PSS.lfp.mat']);
save(savefile,'PSS','Osci','PSScorr','Oscicorr','usechannels','PSSstatsdepth','-v7.3');

%% FIGURE
figure;
subplot(1,4,1);
imagesc(PSSstatsdepth.bins,1:length(usechannels),PSSstatsdepth.dist)
ColorbarWithAxis([min(min(PSSstatsdepth.dist)) max(max(PSSstatsdepth.dist))],['counts'])
caxis([min(min(PSSstatsdepth.dist)) max(max(PSSstatsdepth.dist))])
xlabel('PSS distributions (au)')
ylabel('channel no. (depth-aligned)')
colormap(gca,'jet')
axis tight

subplot(1,4,2);
plot(PSScorr.EMG,1:length(usechannels),'k','linewidth',2)
hold on;
plot(PSScorr.pup,1:length(usechannels),'r','linewidth',2)
legend('PSS-EMG corr','PSS-pupil diameter corr','location','eastoutside')
xlabel('PSS-Behavior correlation');
axis tight

subplot(1,4,3);
imagesc(1:length(usechannels),log10(Frac.freq),Oscicorr.EMG')
LogScale('x',10);
colormap(gca,rwbcolormap)
ColorbarWithAxis([min(min(Oscicorr.EMG)) max(max(Oscicorr.EMG))],['Spearman corr'])
caxis([min(min(Oscicorr.EMG)) max(max(Oscicorr.EMG))])
xlabel('f (Hz)'); ylabel('channel no. (depth-aligned)');
title('Oscillatory LFP-EMG correlation');

subplot(1,4,4);
imagesc(1:length(usechannels),log10(Frac.freq),Oscicorr.Pup')
LogScale('x',10);
colormap(gca,rwbcolormap)
ColorbarWithAxis([min(min(Oscicorr.Pup)) max(max(Oscicorr.Pup))],['Spearman corr'])
caxis([min(min(Oscicorr.Pup)) max(max(Oscicorr.Pup))])
xlabel('f (Hz)'); ylabel('channel no. (depth-aligned)');
title('Oscillatory LFP-Pupil diameter correlation');

NiceSave('PSS_Behavior_CorrbyDepth',figfolder,baseName)

%% Xcorr PSS/Oscillatory component by depth
PSSxcorr = zeros(size(PSS,2),size(PSS,2));
PSSxcorr_p = zeros(size(PSS,2),size(PSS,2));
for x = 1:size(PSS,2)
    for y = 1:size(PSS,2)
        [PSSxcorr(x,y),PSSxcorr_p(x,y)] = corr(PSS(:,x),PSS(:,y),...
            'type','spearman','rows','complete');
    end
end

Oscixcorr = zeros(size(Osci,2),size(Osci,2));
Oscixcorr_p = zeros(size(Osci,2),size(Osci,2));
for x = 1:size(Osci,2)
    for y = 1:size(Osci,2)
        [Oscixcorr(x,y),Oscixcorr_p(x,y)] = corr(Osci(:,x),Osci(:,y),...
            'type','spearman','rows','complete');
    end
end

% Saving data...
savefile = fullfile(basePath,[baseName,'.PSSOsciXCorr.PSS.lfp.mat']);
save(savefile,'PSSxcorr','PSSxcorr_p','Oscixcorr','Oscixcorr_p');

%% FIGURE
figure;

subplot(1,2,1);
imagesc(1:length(usechannels),1:length(usechannels),PSSxcorr);
colormap(gca,rwbcolormap)
axis xy
axis tight
ColorbarWithAxis([min(min(PSSxcorr)) max(max(PSSxcorr))],['Spearman corr'])
caxis([min(min(PSSxcorr)) max(max(PSSxcorr))])
xlabel('channel no. (depth-aligned)');ylabel('channel no. (depth-aligned)');
title('PSS-PSS xcorr by depth');

subplot(1,2,1);
imagesc(1:length(usechannels),1:length(usechannels),Oscixcorr);
colormap(gca,rwbcolormap)
axis xy
axis tight
ColorbarWithAxis([min(min(Oscixcorr)) max(max(Oscixcorr))],['Spearman corr'])
caxis([min(min(Oscixcorr)) max(max(Oscixcorr))])
xlabel('channel no. (depth-aligned)');ylabel('channel no. (depth-aligned)');
title('Oscillatory-Oscillatory LFP xcorr by depth');

NiceSave('Frac_Osci_Xcorrbydepth',figfolder,baseName);
