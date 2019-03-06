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

%% Optimizing lower frequency bounds and window/dt
% Loading LFP for SlowWave detection channel
load(fullfile(basePath,[baseName,'.SlowWaves.events.mat']),'SlowWaves');
%SlowWaves = bz_LoadEvents(basePath,'SlowWaves');
channel = SlowWaves.detectorinfo.detectionchannel;
clear SlowWaves;

lfp = bz_GetLFP(channel,'basepath',basePath,'noPrompts',true);

spontidx = find(lfp.timestamps < sponttimes(2));
lfp.data = lfp.data(spontidx);
lfp.timestamps = lfp.timestamps(spontidx);

downsamplefactor = 3;
%lfp = bz_DownsampleLFP(lfp,downsamplefactor);
lfp.samplingRate = lfp.samplingRate./downsamplefactor;
lfp.data = downsample(double(lfp.data),downsamplefactor);
lfp.timestamps = downsample(lfp.timestamps,downsamplefactor);

% Separate fractal and oscillatory components using sliding window
srate = lfp.samplingRate; % sampling frequency SPECIFY...
movingwin = round(([logspace(log10(0.25),log10(100),10); logspace(log10(0.25),log10(100),10)./4].*srate)'); % [window size, sliding step]
lowerbound = logspace(log10(0.5),log10(80),10);

FracEMGrho = zeros(length(lowerbound),size(movingwin,1));
FracEMGp = zeros(length(lowerbound),size(movingwin,1));
FracPupilrho = zeros(length(lowerbound),size(movingwin,1));
FracPupilp = zeros(length(lowerbound),size(movingwin,1));

Fracmeanrsq = zeros(length(lowerbound),size(movingwin,1));
Fracrsqcorr = zeros(length(lowerbound),size(movingwin,1));
Fracrsqp = zeros(length(lowerbound),size(movingwin,1));

for x = 1:size(movingwin,1)
    
    nwin = floor((length(lfp.data) - movingwin(x,1))/movingwin(x,2));
    sig = zeros(movingwin(x,1),nwin);
    timestamp = zeros(nwin,1);
    for i = 1 : nwin
        idx = [ceil((i-1)*movingwin(x,2))+1 : ceil((i-1)*movingwin(x,2))+movingwin(x,1)];
        sig(:,i) = lfp.data(idx);
        %figure out timestamp associated with window i 
        timestamp(i) = mean(lfp.timestamps(idx));
    end
    
    Frac = amri_sig_fractal_gpu(sig,srate,'detrend',1);
    %st = movingwin(x,1)/(srate*2);
    %Frac.timestamps = st + ((0:nwin-1) * movingwin(x,2)/srate);
    Frac.timestamps = timestamp;
    
    % Interpolating...
    PSS.EMG = interp1(EMGwhisk.timestamps,EMGwhisk.EMGenvelope,...
        Frac.timestamps,'nearest');
    PSS.pupilsize = interp1(pupildilation.timestamps,pupildilation.data,...
        Frac.timestamps,'nearest');
    
    for f = 1:length(lowerbound)
        Frange = [lowerbound(f), 100]; % define frequency range for power-law fitting
        Frac = amri_sig_plawfit(Frac,Frange);
        
        %
        [FracEMGrho(f,x),FracEMGp(f,x)] = corr(log10(PSS.EMG)',...
            Frac.Beta.*-1,'type','spearman','rows','complete');
        
        [FracPupilrho(f,x),FracPupilp(f,x)] = corr(log10(PSS.pupilsize)',...
            Frac.Beta.*-1,'type','spearman','rows','complete');
        
        %
        Fracmeanrsq(f,x) = nanmean(Frac.rsq);
        
        [Fracrsqcorr(f,x),Fracrsqp(f,x)] = corr(Frac.rsq,...
            Frac.Beta.*-1,'type','spearman','rows','complete');
    end
end

%% FIGURE
figure;

subplot(2,2,1); hold on;
winds = movingwin(:,1)./srate;
imagesc(log10(lowerbound),log10(winds),FracEMGrho')
[row,col] = find(FracEMGp > 0.05);
plot(log10(lowerbound(row)),log10(winds(col)),'.w')
colormap(gca,'jet')
LogScale('x',10); LogScale('y',10);
xticks(log10([1 2.5 5 10 20 40 80]));
xticklabels({'1','2.5','5','10','20','40','80'});
yticks(log10([0.25 0.5 1 2.5 5 15 30 60 90]));
yticklabels({'0.25','0.5','1','2.5','5','15','30','60','90'});
axis square
axis tight
ColorbarWithAxis([min(min(FracEMGrho)) max(max(FracEMGrho))],['Spearman corr'])
caxis([min(min(FracEMGrho)) max(max(FracEMGrho))])
xlabel('lower f bound (Hz)');ylabel('interval window (s)');
title(['Frac-EMG correlation']);

subplot(2,2,2); hold on;
imagesc(log10(lowerbound),log10(winds),FracPupilrho')
[row,col] = find(FracPupilp > 0.05);
plot(log10(lowerbound(row)),log10(winds(col)),'.w')
colormap(gca,'jet')
LogScale('x',10); LogScale('y',10);
xticks(log10([1 2.5 5 10 20 40 80]));
xticklabels({'1','2.5','5','10','20','40','80'});
yticks(log10([0.25 0.5 1 2.5 5 15 30 60 90]));
yticklabels({'0.25','0.5','1','2.5','5','15','30','60','90'});
axis square
axis tight
ColorbarWithAxis([min(min(FracPupilrho)) max(max(FracPupilrho))],['Spearman corr'])
caxis([min(min(FracPupilrho)) max(max(FracPupilrho))])
xlabel('lower f bound (Hz)'); ylabel('interval window (s)');
title('Frac-Pupil diameter correlation');

subplot(2,2,3); hold on;
imagesc(log10(lowerbound),log10(winds),Fracmeanrsq')
colormap(gca,'jet')
LogScale('x',10); LogScale('y',10);
xticks(log10([1 2.5 5 10 20 40 80]));
xticklabels({'1','2.5','5','10','20','40','80'});
yticks(log10([0.25 0.5 1 2.5 5 15 30 60 90]));
yticklabels({'0.25','0.5','1','2.5','5','15','30','60','90'});
axis square
axis tight
ColorbarWithAxis([min(min(Fracmeanrsq)) max(max(Fracmeanrsq))],['au'])
caxis([min(min(Fracmeanrsq)) max(max(Fracmeanrsq))])
xlabel('lower f bound (Hz)');ylabel('interval window (s)');
title(['Mean RSQ slope fit']);

subplot(2,2,4); hold on;
imagesc(log10(lowerbound),log10(winds),Fracrsqcorr')
[row,col] = find(Fracrsqp > 0.05);
plot(log10(lowerbound(row)),log10(winds(col)),'.w')
colormap(gca,'jet')
LogScale('x',10); LogScale('y',10);
xticks(log10([1 2.5 5 10 20 40 80]));
xticklabels({'1','2.5','5','10','20','40','80'});
yticks(log10([0.25 0.5 1 2.5 5 15 30 60 90]));
yticklabels({'0.25','0.5','1','2.5','5','15','30','60','90'});
axis square
axis tight
ColorbarWithAxis([-1 1],['Spearman corr'])
caxis([-1 1])
xlabel('lower f bound (Hz)'); ylabel('interval window (s)');
title('Frac-RSQ correlation');

NiceSave('Frac_EMG_Pupil_RSQ_Corr_p_win_lof',figfolder,baseName);

%% Example slow/fast PSS, RSQ, EMG and Pupil...
% Slow PSS
movingwin = round([15 3.75].*srate);
Frange = [2.5,100];

nwin = floor((length(lfp.data) - movingwin(1))/movingwin(2));
sig = zeros(movingwin(1),nwin);
timestamp = zeros(nwin,1);
for i = 1 : nwin
    idx = [ceil((i-1)*movingwin(2))+1 : ceil((i-1)*movingwin(2))+movingwin(1)];
    sig(:,i) = lfp.data(idx);
    %figure out timestamp associated with window i
    timestamp(i) = mean(lfp.timestamps(idx));
end
Frac_slow = amri_sig_fractal_gpu(sig,srate,'detrend',1);
Frac_slow.timestamps = timestamp;
Frac_slow = amri_sig_plawfit(Frac_slow,Frange);

% Fast PSS
movingwin = round([1 0.25].*srate);
Frange = [5,100];

nwin = floor((length(lfp.data) - movingwin(1))/movingwin(2));
sig = zeros(movingwin(1),nwin);
timestamp = zeros(nwin,1);
for i = 1 : nwin
    idx = [ceil((i-1)*movingwin(2))+1 : ceil((i-1)*movingwin(2))+movingwin(1)];
    sig(:,i) = lfp.data(idx);
    %figure out timestamp associated with window i
    timestamp(i) = mean(lfp.timestamps(idx));
end
Frac_fast = amri_sig_fractal_gpu(sig,srate,'detrend',1);
Frac_fast.timestamps = timestamp;
Frac_fast = amri_sig_plawfit(Frac_fast,Frange);

% Interpolating...
PSS.EMG = interp1(EMGwhisk.timestamps,EMGwhisk.EMGenvelope,...
    Frac_fast.timestamps,'nearest');
PSS.pupilsize = interp1(pupildilation.timestamps,pupildilation.data,...
    Frac_fast.timestamps,'nearest');

%% FIGURE
figure;

h(1) = subplot(9,1,1); hold on;
plot(Frac_fast.timestamps,Frac_fast.rsq,'k');
plot(Frac_slow.timestamps,Frac_slow.rsq,'r');
ylim([0 1]);
legend({'RSQ fast','RSQ slow'},'location','southeast')

h(2) = subplot(9,1,2:3); hold on;
plot(Frac_fast.timestamps,Frac_fast.Beta.*-1,'k');
plot(Frac_slow.timestamps,Frac_slow.Beta.*-1,'r');
legend({'Frac fast','Frac slow'},'location','southeast')

h(3) = subplot(9,1,4:5); hold on;
plot(Frac_fast.timestamps,PSS.EMG,'k');
legend({'EMG'},'location','southeast')

h(4) = subplot(9,1,6:7); hold on;
plot(Frac_fast.timestamps,PSS.pupilsize,'k');
legend({'Pupil diameter'},'location','southeast')

h(5) = subplot(9,1,8:9); hold on;
imagesc(Frac_fast.timestamps,Frac_fast.freq,Frac_fast.osci);
caxis([min(min(Frac_fast.osci)) max(max(Frac_fast.osci))])
ylabel('f (hz)'); ylim([1 120]);
colormap(gca,'jet');
legend({'Oscillatory PSpec'},'location','southeast')

linkaxes(h,'x');
xlim([250 650]);
xlabel('time (s)');

NiceSave('Example_PSS_rsq_osci_behavior',figfolder,baseName);

%% Correlating PSS/Oscillatory component to EMG/Pupil diameter by depth
badchannels = sessionInfo.badchannels;
usechannels = sessionInfo.AnatGrps.Channels;
usechannels(ismember(usechannels,badchannels))=[];

% Assuming that LFP still remains loaded from prior analysis...
movingwin = round([1 0.25].*srate);

nwin = floor((length(lfp.data) - movingwin(1))/movingwin(2));
sig = zeros(movingwin(1),nwin);
timestamp = zeros(nwin,1);
for i = 1 : nwin
    idx = [ceil((i-1)*movingwin(2))+1 : ceil((i-1)*movingwin(2))+movingwin(1)];
    sig(:,i) = lfp.data(idx);
    %figure out timestamp associated with window i
    timestamp(i) = mean(lfp.timestamps(idx));
end
clear lfp

Frac = amri_sig_fractal_gpu(sig,srate,'detrend',1);
Frac.timestamps = timestamp;
Frac = amri_sig_plawfit(Frac,Frange);

% Selected from prior analysis
Frange = [5, 100]; % define frequency range for power-law fitting

% Now calculating corrs...
PSS = []; Osci = [];

PSScorr.EMG = zeros(length(usechannels),1);
PSScorr.EMG_p = zeros(length(usechannels),1);
PSScorr.Pup = zeros(length(usechannels),1);
PSScorr.Pup_p = zeros(length(usechannels),1);

dcorr.EMG = zeros(length(usechannels),1);
dcorr.EMG_p = zeros(length(usechannels),1);
dcorr.Pup = zeros(length(usechannels),1);
dcorr.Pup_p = zeros(length(usechannels),1);

tcorr.EMG = zeros(length(usechannels),1);
tcorr.EMG_p = zeros(length(usechannels),1);
tcorr.Pup = zeros(length(usechannels),1);
tcorr.Pup_p = zeros(length(usechannels),1);

gcorr.EMG = zeros(length(usechannels),1);
gcorr.EMG_p = zeros(length(usechannels),1);
gcorr.Pup = zeros(length(usechannels),1);
gcorr.Pup_p = zeros(length(usechannels),1);

Oscicorr.EMG = zeros(size(Frac.osci,1),length(usechannels));
Oscicorr.EMG_p = zeros(size(Frac.osci,1),length(usechannels));
Oscicorr.Pup = zeros(size(Frac.osci,1),length(usechannels));
Oscicorr.Pup_p = zeros(size(Frac.osci,1),length(usechannels));

nbins = 100;
PSSstatsdepth.bins = linspace(-5,0,nbins);
PSSstatsdepth.dist = zeros(length(usechannels),nbins);

for cc = 1:length(usechannels)
    cc
    channum = usechannels(cc);
    lfp = bz_GetLFP(channum,'basepath',basePath,'noPrompts',true);
    
    lfp.data = lfp.data(spontidx);
    lfp.timestamps = lfp.timestamps(spontidx);
    
    lfp.samplingRate = lfp.samplingRate./downsamplefactor;
    lfp.data = downsample(double(lfp.data),downsamplefactor);
    lfp.timestamps = downsample(lfp.timestamps,downsamplefactor);
    
    %% Deconstruction
    sig = zeros(movingwin(1),nwin);
    timestamp = zeros(nwin,1);
    for i = 1 : nwin
        idx = [ceil((i-1)*movingwin(2))+1 : ceil((i-1)*movingwin(2))+movingwin(1)];
        sig(:,i) = lfp.data(idx);
        %figure out timestamp associated with window i
        timestamp(i) = mean(lfp.timestamps(idx));
    end
    clear lfp
    
    Frac = amri_sig_fractal_gpu(sig,srate,'detrend',1);
    Frac.timestamps = timestamp;
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
    
    deltaidx = find(Frac.freq >= 0.5 & Frac.freq <= 3);
    [dcorr.EMG(cc),dcorr.EMG_p(cc)] =...
        corr(log10(temp_EMG)',(nanmean(Frac.osci(deltaidx,:),1))',...
        'type','spearman','rows','complete');
    [dcorr.Pup(cc),dcorr.Pup_p(cc)] =...
        corr(log10(temp_Pup)',(nanmean(Frac.osci(deltaidx,:),1))',...
        'type','spearman','rows','complete');
    
    thetaidx = find(Frac.freq >= 4 & Frac.freq <= 10);
    [tcorr.EMG(cc),tcorr.EMG_p(cc)] =...
        corr(log10(temp_EMG)',(nanmean(Frac.osci(thetaidx,:),1))',...
        'type','spearman','rows','complete');
    [tcorr.Pup(cc),tcorr.Pup_p(cc)] =...
        corr(log10(temp_Pup)',(nanmean(Frac.osci(thetaidx,:),1))',...
        'type','spearman','rows','complete');
    
    gammaidx = find(Frac.freq >= 40 & Frac.freq <= 120);
    [gcorr.EMG(cc),gcorr.EMG_p(cc)] =...
        corr(log10(temp_EMG)',(nanmean(Frac.osci(gammaidx,:),1))',...
        'type','spearman','rows','complete');
    [gcorr.Pup(cc),gcorr.Pup_p(cc)] =...
        corr(log10(temp_Pup)',(nanmean(Frac.osci(gammaidx,:),1))',...
        'type','spearman','rows','complete');
    
    for x = 1:size(Frac.osci,1)
        [Oscicorr.EMG(x,cc),Oscicorr.EMG_p(x,cc)] =...
            corr(log10(temp_EMG)',Frac.osci(x,:)',...
            'type','spearman','rows','complete');
        [Oscicorr.Pup(x,cc),Oscicorr.Pup_p(x,cc)] =...
            corr(log10(temp_Pup)',Frac.osci(x,:)',...
            'type','spearman','rows','complete');
    end
    
end

%% FIGURE
figure;
subplot(2,3,1);
imagesc(PSSstatsdepth.bins,1:length(usechannels),PSSstatsdepth.dist)
ColorbarWithAxis([min(min(PSSstatsdepth.dist)) max(max(PSSstatsdepth.dist))],['counts'])
caxis([min(min(PSSstatsdepth.dist)) max(max(PSSstatsdepth.dist))])
xlabel('PSS distributions (au)')
xlim([-4 -1]);
ylabel('channel no. (depth-aligned)')
colormap(gca,'jet')
axis tight

subplot(2,3,2);
plot(PSScorr.EMG,1:length(usechannels),'k','linewidth',2)
hold on;
plot(PSScorr.Pup,1:length(usechannels),'r','linewidth',2)
legend('PSS-EMG corr','PSS-pupil diameter corr','location','eastoutside')
xlabel('PSS-Behavior correlation');
set(gca,'ydir','reverse')
axis tight

subplot(2,3,4);
imagesc(log10(Frac.freq),1:length(usechannels),Oscicorr.EMG')
LogScale('x',10);
colormap(gca,'jet')
ColorbarWithAxis([min(min(Oscicorr.EMG)) max(max(Oscicorr.EMG))],['Spearman corr'])
caxis([min(min(Oscicorr.EMG)) max(max(Oscicorr.EMG))])
xlabel('f (Hz)'); ylabel('channel no. (depth-aligned)');
title('Oscillatory LFP-EMG correlation');

subplot(2,3,5);
imagesc(log10(Frac.freq),1:length(usechannels),Oscicorr.Pup')
LogScale('x',10);
colormap(gca,'jet')
ColorbarWithAxis([min(min(Oscicorr.Pup)) max(max(Oscicorr.Pup))],['Spearman corr'])
caxis([min(min(Oscicorr.Pup)) max(max(Oscicorr.Pup))])
xlabel('f (Hz)'); ylabel('channel no. (depth-aligned)');
title('Oscillatory LFP-Pupil diameter correlation');

subplot(2,6,5);
plot(dcorr.EMG,1:length(usechannels),'r','linewidth',2)
hold on;
plot(tcorr.EMG,1:length(usechannels),'r','linewidth',2)
plot(gcorr.EMG,1:length(usechannels),'r','linewidth',2)
legend('d','t','g','location','southeast')
xlabel('Oscillatory-EMG correlation');
set(gca,'ydir','reverse')
axis tight

subplot(2,6,6);
plot(dcorr.Pup,1:length(usechannels),'r','linewidth',2)
hold on;
plot(tcorr.Pup,1:length(usechannels),'r','linewidth',2)
plot(gcorr.Pup,1:length(usechannels),'r','linewidth',2)
legend('d','t','g','location','southeast')
xlabel('Oscillatory-Pupil diameter correlation');
set(gca,'ydir','reverse')
axis tight

NiceSave('fPSS_Behavior_CorrbyDepth',figfolder,baseName)

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

%% FIGURE
figure;

subplot(1,2,1);
imagesc(1:length(usechannels),1:length(usechannels),PSSxcorr);
colormap(gca,'jet')
axis square
%set(gca,'ydir','reverse')
ColorbarWithAxis([min(min(PSSxcorr)) max(max(PSSxcorr))],['Spearman corr'])
caxis([min(min(PSSxcorr)) max(max(PSSxcorr))])
xlabel('channel no. (depth-aligned)');ylabel('channel no. (depth-aligned)');
title('PSS-PSS xcorr by depth');

subplot(1,2,2);
imagesc(1:length(usechannels),1:length(usechannels),Oscixcorr);
colormap(gca,'jet')
axis square
%set(gca,'ydir','reverse')
ColorbarWithAxis([min(min(Oscixcorr)) max(max(Oscixcorr))],['Spearman corr'])
caxis([min(min(Oscixcorr)) max(max(Oscixcorr))])
xlabel('channel no. (depth-aligned)');ylabel('channel no. (depth-aligned)');
title('Oscillatory-Oscillatory LFP xcorr by depth');

NiceSave('fPSS_Osci_Xcorrbydepth',figfolder,baseName);

%% Again... but slow!
lfp = bz_GetLFP(channel,'basepath',basePath,'noPrompts',true);
lfp.data = lfp.data(spontidx);
lfp.timestamps = lfp.timestamps(spontidx);

lfp.samplingRate = lfp.samplingRate./downsamplefactor;
lfp.data = downsample(double(lfp.data),downsamplefactor);
lfp.timestamps = downsample(lfp.timestamps,downsamplefactor);

movingwin = round([15 3.75].*srate);

nwin = floor((length(lfp.data) - movingwin(1))/movingwin(2));
sig = zeros(movingwin(1),nwin);
timestamp = zeros(nwin,1);
for i = 1 : nwin
    idx = [ceil((i-1)*movingwin(2))+1 : ceil((i-1)*movingwin(2))+movingwin(1)];
    sig(:,i) = lfp.data(idx);
    %figure out timestamp associated with window i
    timestamp(i) = mean(lfp.timestamps(idx));
end
clear lfp

Frac = amri_sig_fractal_gpu(sig,srate,'detrend',1);
Frac.timestamps = timestamp;
Frac = amri_sig_plawfit(Frac,Frange);

% Selected from prior analysis
Frange = [2.5, 100]; % define frequency range for power-law fitting

% Now calculating corrs...
PSS = []; Osci = [];

PSScorr.EMG = zeros(length(usechannels),1);
PSScorr.EMG_p = zeros(length(usechannels),1);
PSScorr.Pup = zeros(length(usechannels),1);
PSScorr.Pup_p = zeros(length(usechannels),1);

Oscicorr.EMG = zeros(size(Frac.osci,1),length(usechannels));
Oscicorr.EMG_p = zeros(size(Frac.osci,1),length(usechannels));
Oscicorr.Pup = zeros(size(Frac.osci,1),length(usechannels));
Oscicorr.Pup_p = zeros(size(Frac.osci,1),length(usechannels));

nbins = 100;
PSSstatsdepth.bins = linspace(-5,0,nbins);
PSSstatsdepth.dist = zeros(length(usechannels),nbins);

for cc = 1:length(usechannels)
    cc
    channum = usechannels(cc);
    lfp = bz_GetLFP(channum,'basepath',basePath,'noPrompts',true);
    
    lfp.data = lfp.data(spontidx);
    lfp.timestamps = lfp.timestamps(spontidx);
    
    lfp.samplingRate = lfp.samplingRate./downsamplefactor;
    lfp.data = downsample(double(lfp.data),downsamplefactor);
    lfp.timestamps = downsample(lfp.timestamps,downsamplefactor);
    
    %% Deconstruction
    sig = zeros(movingwin(1),nwin);
    timestamp = zeros(nwin,1);
    for i = 1 : nwin
        idx = [ceil((i-1)*movingwin(2))+1 : ceil((i-1)*movingwin(2))+movingwin(1)];
        sig(:,i) = lfp.data(idx);
        %figure out timestamp associated with window i
        timestamp(i) = mean(lfp.timestamps(idx));
    end
    clear lfp
    
    Frac = amri_sig_fractal_gpu(sig,srate,'detrend',1);
    Frac.timestamps = timestamp;
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
            corr(log10(temp_EMG)',Frac.osci(x,:)',...
            'type','spearman','rows','complete');
        [Oscicorr.Pup(x,cc),Oscicorr.Pup_p(x,cc)] =...
            corr(log10(temp_Pup)',Frac.osci(x,:)',...
            'type','spearman','rows','complete');
    end
end

%% FIGURE
figure;
subplot(2,2,1);
imagesc(PSSstatsdepth.bins,1:length(usechannels),PSSstatsdepth.dist)
ColorbarWithAxis([min(min(PSSstatsdepth.dist)) max(max(PSSstatsdepth.dist))],['counts'])
caxis([min(min(PSSstatsdepth.dist)) max(max(PSSstatsdepth.dist))])
xlabel('PSS distributions (au)')
xlim([-4 -1]);
ylabel('channel no. (depth-aligned)')
colormap(gca,'jet')
axis tight

subplot(2,2,2);
plot(PSScorr.EMG,1:length(usechannels),'k','linewidth',2)
hold on;
plot(PSScorr.Pup,1:length(usechannels),'r','linewidth',2)
legend('PSS-EMG corr','PSS-pupil diameter corr','location','eastoutside')
xlabel('PSS-Behavior correlation');
set(gca,'ydir','reverse')
axis tight

subplot(2,2,3);
imagesc(log10(Frac.freq),1:length(usechannels),Oscicorr.EMG')
LogScale('x',10);
colormap(gca,'jet')
ColorbarWithAxis([min(min(Oscicorr.EMG)) max(max(Oscicorr.EMG))],['Spearman corr'])
caxis([min(min(Oscicorr.EMG)) max(max(Oscicorr.EMG))])
xlabel('f (Hz)'); ylabel('channel no. (depth-aligned)');
title('Oscillatory LFP-EMG correlation');

subplot(2,2,4);
imagesc(log10(Frac.freq),1:length(usechannels),Oscicorr.Pup')
LogScale('x',10);
colormap(gca,'jet')
ColorbarWithAxis([min(min(Oscicorr.Pup)) max(max(Oscicorr.Pup))],['Spearman corr'])
caxis([min(min(Oscicorr.Pup)) max(max(Oscicorr.Pup))])
xlabel('f (Hz)'); ylabel('channel no. (depth-aligned)');
title('Oscillatory LFP-Pupil diameter correlation');

NiceSave('sPSS_Behavior_CorrbyDepth',figfolder,baseName)

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

%% FIGURE
figure;

subplot(1,2,1);
imagesc(1:length(usechannels),1:length(usechannels),PSSxcorr);
colormap(gca,'jet')
axis square
%set(gca,'ydir','reverse')
ColorbarWithAxis([min(min(PSSxcorr)) max(max(PSSxcorr))],['Spearman corr'])
caxis([min(min(PSSxcorr)) max(max(PSSxcorr))])
xlabel('channel no. (depth-aligned)');ylabel('channel no. (depth-aligned)');
title('PSS-PSS xcorr by depth');

subplot(1,2,2);
imagesc(1:length(usechannels),1:length(usechannels),Oscixcorr);
colormap(gca,'jet')
axis square
%set(gca,'ydir','reverse')
ColorbarWithAxis([min(min(Oscixcorr)) max(max(Oscixcorr))],['Spearman corr'])
caxis([min(min(Oscixcorr)) max(max(Oscixcorr))])
xlabel('channel no. (depth-aligned)');ylabel('channel no. (depth-aligned)');
title('Oscillatory-Oscillatory LFP xcorr by depth');

NiceSave('sPSS_Osci_Xcorrbydepth',figfolder,baseName);
