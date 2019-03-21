basePath = pwd;
baseName = bz_BasenameFromBasepath(basePath);
sessionInfo = bz_getSessionInfo(basePath,'noPrompts',true);

% badchannels = sessionInfo.badchannels;
% usechannels = sessionInfo.AnatGrps.Channels;
% usechannels(ismember(usechannels,badchannels))=[];
channels = sessionInfo.channels;

% Obtaining normalized layer info
rescaled = rescaleCx(basePath);
usechannels = rescaled.channels(~isnan(rescaled.ndepth));
normdepth = rescaled.ndepth(~isnan(rescaled.ndepth));

L1idx = find(rescaled.ndepth >= 0 & rescaled.ndepth <= 0.1);
L1idx = rescaled.channels(L1idx);

L23idx = find(rescaled.ndepth >= 0.1 & rescaled.ndepth <= 0.35);
L23idx = rescaled.channels(L23idx);

L4idx = find(rescaled.ndepth >= 0.35 & rescaled.ndepth <= 0.5);
L4idx = rescaled.channels(L4idx);

L5aidx = find(rescaled.ndepth >= 0.5 & rescaled.ndepth <= 0.6);
L5aidx = rescaled.channels(L5aidx);

L56idx = find(rescaled.ndepth >= 0.6 & rescaled.ndepth <= 0.9);
L56idx = rescaled.channels(L56idx);

L6idx = find(rescaled.ndepth >= 0.9 & rescaled.ndepth <= 1);
L6idx = rescaled.channels(L6idx);

% Specifying folders
figfolder = fullfile(basePath,'AnalysisFigures');
savefile = fullfile(basePath,[baseName,'.LaminarSpectralAnalysis.lfp.mat']);
savefolder = fullfile(basePath,'WaveSpec');

%Pending: better layer boundary detection and exclusion of bad channels

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

smoothwin = 2; %s
pupildilation.dpdt = diff(smooth(pupildilation.data,smoothwin.*pupildilation.samplingRate,'moving')).*pupildilation.samplingRate;
pupildilation.dpdt = smooth(pupildilation.dpdt,smoothwin.*pupildilation.samplingRate,'moving');
pupildilation.timestamps = pupildilation.timestamps(~nantimes);

% EMG
EMGwhisk = bz_LoadBehavior(basePath,'EMGwhisk');

% Specifying SPONT whisking
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

% Touch
% Piezotouch = bz_LoadBehavior(basePath,'Piezotouch');

% Pupil phase
lowfilter = [0.01 0.1];
pupil4filter = pupildilation;
lowpupildata = bz_Filter(pupil4filter,'passband',lowfilter,'filter' ,'fir1','order',3);
x = lowpupildata.timestamps;
x = repmat(x,[1 5]);
y = repmat((0:4) * (1/(lowpupildata.samplingRate*5)) ,[length(x) 1]);
x = x + y;
y = reshape(x',[length(x)*5 1]);
lowpupildata.timestamps = y;
lowpupildata.amp = resample(lowpupildata.amp,5,1);

pupthresh = nanmedian(log10(lowpupildata.amp));
highpup = log10(lowpupildata.amp)>pupthresh;

% Getting intervals in spec times
load(fullfile(savefolder,[baseName,'.',num2str(0),'.WaveSpec.lfp.mat']));

eventshipupil = interp1(wavespec.timestamps,...
    wavespec.timestamps,...
    lowpupildata.timestamps(highpup),'nearest').*wavespec.samplingRate;

eventslopupil = interp1(wavespec.timestamps,...
    wavespec.timestamps,...
    lowpupildata.timestamps(~highpup),'nearest').*wavespec.samplingRate;

events_Wh = interp1(wavespec.timestamps,wavespec.timestamps,...
    EMGwhisk.ints.Wh,'nearest');
allidx_Wh = [];
for e = 1:size(events_Wh,1)
    allidx_Wh = cat(1,allidx_Wh,wavespec.timestamps(find([wavespec.timestamps >= events_Wh(e,1)...
        & wavespec.timestamps <= events_Wh(e,2)])));
end
allidx_Wh = allidx_Wh.*wavespec.samplingRate;

events_NWh = interp1(wavespec.timestamps,wavespec.timestamps,...
    EMGwhisk.ints.NWh,'nearest');
allidx_NWh = [];
for e = 1:size(events_NWh,1)
    allidx_NWh = cat(1,allidx_NWh,wavespec.timestamps(find([wavespec.timestamps >= events_NWh(e,1)...
        & wavespec.timestamps <= events_NWh(e,2)])));
end
allidx_NWh = allidx_NWh.*wavespec.samplingRate;

% if ~isempty(Piezotouch)
%     events_T = interp1(wavespec.timestamps,wavespec.timestamps,...
%         Piezotouch.ints.Touch,'nearest');
%     allidx_T = [];
%     for e = 1:size(events_T,1)
%         allidx_T = cat(1,allidx_T,wavespec.timestamps(find([wavespec.timestamps >= events_T(e,1)...
%             & wavespec.timestamps <= events_T(e,2)])));
%     end
%     allidx_T = allidx_T.*wavespec.samplingRate;
% else
% end

%% Laminar Power spectra by layer, by state
cLayerSpec_all = NaN(size(wavespec.data,2),length(channels)); 
cLayerSpec_Wh = NaN(size(wavespec.data,2),length(channels)); 
cLayerSpec_NWh = NaN(size(wavespec.data,2),length(channels));  
cLayerSpec_loP = NaN(size(wavespec.data,2),length(channels)); 
cLayerSpec_hiP = NaN(size(wavespec.data,2),length(channels)); 

cLayerSpec_all_n = NaN(size(wavespec.data,2),length(channels)); 
cLayerSpec_Wh_n = NaN(size(wavespec.data,2),length(channels)); 
cLayerSpec_NWh_n = NaN(size(wavespec.data,2),length(channels));  
cLayerSpec_loP_n = NaN(size(wavespec.data,2),length(channels)); 
cLayerSpec_hiP_n = NaN(size(wavespec.data,2),length(channels)); 

cLayerSpec_all_z = NaN(size(wavespec.data,2),length(channels)); 
cLayerSpec_Wh_z = NaN(size(wavespec.data,2),length(channels)); 
cLayerSpec_NWh_z = NaN(size(wavespec.data,2),length(channels));  
cLayerSpec_loP_z = NaN(size(wavespec.data,2),length(channels)); 
cLayerSpec_hiP_z = NaN(size(wavespec.data,2),length(channels)); 

for i = 1:length(channels)
    i
    % Loading spectrograms
    load(fullfile(savefolder,[baseName,'.',num2str(channels(i)),'.WaveSpec.lfp.mat']));
    
    wavespec.dataz = NormToInt(log10(abs(wavespec.data)),'modZ');
    wavespec.datan = log10(abs(wavespec.data))./nanmedian(log10(abs(wavespec.data)),1);
    %wavespec.data = NormToInt(log10(abs(wavespec.data)),'Z',...
    %    [events_Wh; events_NWh],wavespec.samplingRate);
    
    cLayerSpec_all(:,i) = nanmedian(wavespec.data,1);
    cLayerSpec_Wh(:,i) = nanmedian(wavespec.data(round(allidx_Wh),:),1);
    cLayerSpec_NWh(:,i) = nanmedian(wavespec.data(round(allidx_NWh),:),1);
    cLayerSpec_hiP(:,i) = nanmedian(wavespec.data(round(eventshipupil),:),1);
    cLayerSpec_loP(:,i) = nanmedian(wavespec.data(round(eventslopupil),:),1);
    
    cLayerSpec_all_n(:,i) = nanmedian(wavespec.datan,1);
    cLayerSpec_Wh_n(:,i) = nanmedian(wavespec.datan(round(allidx_Wh),:),1);
    cLayerSpec_NWh_n(:,i) = nanmedian(wavespec.datan(round(allidx_NWh),:),1);
    cLayerSpec_hiP_n(:,i) = nanmedian(wavespec.datan(round(eventshipupil),:),1);
    cLayerSpec_loP_n(:,i) = nanmedian(wavespec.datan(round(eventslopupil),:),1);
    
    cLayerSpec_all_z(:,i) = nanmedian(wavespec.dataz,1);
    cLayerSpec_Wh_z(:,i) = nanmedian(wavespec.dataz(round(allidx_Wh),:),1);
    cLayerSpec_NWh_z(:,i) = nanmedian(wavespec.dataz(round(allidx_NWh),:),1);
    cLayerSpec_hiP_z(:,i) = nanmedian(wavespec.dataz(round(eventshipupil),:),1);
    cLayerSpec_loP_z(:,i) = nanmedian(wavespec.dataz(round(eventslopupil),:),1);

%     if ~isempty(Piezotouch)
%         cLayerSpec_T(:,i) = nanmedian(wavespec.data(round(allidx_T),:),1);
%     else
%     end
end

% Saving to struct
% LayerSpectral.cLayerSpec_all = cLayerSpec_all;
% LayerSpectral.cLayerSpec_Wh = cLayerSpec_Wh;
% LayerSpectral.cLayerSpec_NWh = cLayerSpec_NWh;
% LayerSpectral.cLayerSpec_loP = cLayerSpec_loP;
% LayerSpectral.cLayerSpec_hiP = cLayerSpec_hiP;

%% Laminar MUA power (0.15 - 2 kHz)
load(fullfile(basePath,[baseName,'.MUA.lfp.mat']));

cLayerMUA_all = nanmean(MUA.data,1);

% for Wh
events = interp1(MUA.timestamps,MUA.timestamps,...
    EMGwhisk.ints.Wh,'nearest');
allidx = [];
for e = 1:size(events,1)
    allidx = cat(1,allidx,wavespec.timestamps(find([wavespec.timestamps >= events(e,1)...
        & wavespec.timestamps <= events(e,2)])));
end
allidx = allidx.*MUA.samplingRate;
cLayerMUA_Wh = nanmean(MUA.data(round(allidx),:),1);
        
% for NWh
events = interp1(MUA.timestamps,MUA.timestamps,...
    EMGwhisk.ints.NWh,'nearest');
allidx = [];
for e = 1:size(events,1)
    allidx = cat(1,allidx,wavespec.timestamps(find([wavespec.timestamps >= events(e,1)...
        & wavespec.timestamps <= events(e,2)])));
end
allidx = allidx.*MUA.samplingRate;
cLayerMUA_NWh = nanmean(MUA.data(round(allidx),:),1);
        
% if ~isempty(Piezotouch)
%     % for Touch
%     events = interp1(MUA.timestamps,MUA.timestamps,...
%         Piezotouch.ints.Touch,'nearest');
%     allidx = [];
%     for e = 1:size(events,1)
%         allidx = cat(1,allidx,wavespec.timestamps(find([wavespec.timestamps >= events(e,1)...
%             & wavespec.timestamps <= events(e,2)])));
%     end
%     allidx = allidx.*MUA.samplingRate;
%     cLayerMUA_T = nanmean(MUA.data(round(allidx),:),1);
% else
% end

events = interp1(MUA.timestamps,MUA.timestamps,...
    lowpupildata.timestamps(highpup),'nearest').*MUA.samplingRate;
cLayerMUA_hiP = nanmean(MUA.data(round(events)),1);

events = interp1(MUA.timestamps,MUA.timestamps,...
    lowpupildata.timestamps(~highpup),'nearest').*MUA.samplingRate;
cLayerMUA_loP = nanmean(MUA.data(round(events)),1);

% event MUA
eventMUA_Wh = eventMUA(MUA,EMGwhisk.ints.Wh(:,1),...
    'channels',channels,'twin',[0.75 0.75],'spat_sm',0,'saveMat',false);

% Saving to struct
LayerSpectral.cLayerMUA_all = cLayerMUA_all;
LayerSpectral.cLayerMUA_Wh = cLayerMUA_Wh;
LayerSpectral.cLayerMUA_NWh = cLayerMUA_NWh;
LayerSpectral.cLayerMUA_hiP = cLayerMUA_hiP;
LayerSpectral.cLayerMUA_loP = cLayerMUA_loP;
LayerSpectral.eventMUA_Wh = eventMUA_Wh;

%% FIGURE: 
figure;
% subplot(1,7,1);
% stairs(cLayerMUA_all(usechannels+1),normdepth,'Color','k'); hold on;
% stairs(cLayerMUA_NWh(usechannels+1),normdepth,'Color','b'); hold on;
% stairs(cLayerMUA_Wh(usechannels+1),normdepth,'Color','r'); hold on;
% set(gca,'YDir','reverse');
% ylim([normdepth(1) normdepth(end)]);
% xlabel('MUA power (0.15-2 KHz)'); 
% set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
% set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
% set(gca,'YGrid','on', 'GridColor','k','GridAlpha',0.45);
% box on;
% legend({'ALL','NWh','Wh'},'location','northeast');

cmin = min(min(cLayerSpec_Wh(:,usechannels+1)-cLayerSpec_NWh(:,usechannels+1)));
cmax = max(max(cLayerSpec_Wh(:,usechannels+1)-cLayerSpec_NWh(:,usechannels+1)));

subplot(1,7,2:3);
imagesc(log10(wavespec.freqs),normdepth,...
    (cLayerSpec_Wh(:,usechannels+1)-cLayerSpec_NWh(:,usechannels+1))');
axis tight
LogScale('x',10)
colormap('jet'); 
%ColorbarWithAxis([cmin cmax],['power (dB)'])
caxis([cmax*-1 cmax]);
%set(gca,'YDir','reverse');
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('Power spectra Wh-NWh diff');

cmin = min([min(min(cLayerSpec_Wh(:,usechannels+1)))...
    min(min(cLayerSpec_NWh(:,usechannels+1)))]);
cmax = max([max(max(cLayerSpec_Wh(:,usechannels+1)))...
    max(max(cLayerSpec_NWh(:,usechannels+1)))]);

subplot(1,7,4:5);
imagesc(log10(wavespec.freqs),normdepth,cLayerSpec_NWh(:,usechannels+1)');
axis tight
LogScale('x',10)
colormap('jet'); 
%ColorbarWithAxis([cmin cmax],['power (dB)'])
caxis([cmax*-1 cmax]);
%set(gca,'YDir','reverse');
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('Power spectra NWh');

subplot(1,7,6:7);
imagesc(log10(wavespec.freqs),normdepth,cLayerSpec_Wh(:,usechannels+1)');
axis tight
LogScale('x',10)
colormap('jet'); 
%ColorbarWithAxis([cmin cmax],['power (dB)'])
caxis([cmax*-1 cmax]);
%set(gca,'YDir','reverse');
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('Power spectra Wh');

NiceSave('LaminarPspec_Wh_NWh',figfolder,baseName)

%% FIGURE:
figure;
% subplot(1,7,1);
% stairs(cLayerMUA_all(usechannels+1),normdepth,'Color','k'); hold on;
% stairs(cLayerMUA_loP(usechannels+1),normdepth,'Color','b'); hold on;
% stairs(cLayerMUA_hiP(usechannels+1),normdepth,'Color','r'); hold on;
% set(gca,'YDir','reverse');
% ylim([normdepth(1) normdepth(end)]);
% xlabel('MUA power (0.15-2 KHz)'); 
% set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
% set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
% set(gca,'YGrid','on', 'GridColor','k','GridAlpha',0.45);
% box on;
% legend({'ALL','loP','hiP'},'location','northeast');

cmin = min(min(cLayerSpec_hiP(:,usechannels+1)-cLayerSpec_loP(:,usechannels+1)));
cmax = max(max(cLayerSpec_hiP(:,usechannels+1)-cLayerSpec_loP(:,usechannels+1)));

subplot(1,7,2:3);
imagesc(log10(wavespec.freqs),normdepth,...
    (cLayerSpec_hiP(:,usechannels+1)-cLayerSpec_loP(:,usechannels+1))');
axis tight
LogScale('x',10)
colormap('jet'); 
%ColorbarWithAxis([cmin cmax],['power (dB)'])
caxis([cmax*-1 cmax]);
%set(gca,'YDir','reverse');
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('Power spectra hi-lo Pupil diff');

cmin = min([min(min(cLayerSpec_loP(:,usechannels+1)))...
    min(min(cLayerSpec_hiP(:,usechannels+1)))]);
cmax = max([max(max(cLayerSpec_loP(:,usechannels+1)))...
    max(max(cLayerSpec_hiP(:,usechannels+1)))]);

subplot(1,7,4:5);
imagesc(log10(wavespec.freqs),normdepth,cLayerSpec_loP(:,usechannels+1)');
axis tight
LogScale('x',10)
colormap('jet'); 
%ColorbarWithAxis([cmin cmax],['power (dB)'])
caxis([cmax*-1 cmax]);
%set(gca,'YDir','reverse');
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('Power spectra <median Pupil');

subplot(1,7,6:7);
imagesc(log10(wavespec.freqs),normdepth,cLayerSpec_hiP(:,usechannels+1)');
axis tight
LogScale('x',10)
colormap('jet'); 
%ColorbarWithAxis([cmin cmax],['power (dB)'])
caxis([cmax*-1 cmax]);
%set(gca,'YDir','reverse');
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('Power spectra >median Pupil');

NiceSave('LaminarPspec_lo_hiPup',figfolder,baseName)

%% State LFP-LFP cross-corr
SpecXcorr_all = zeros(size(cLayerSpec_all,2),size(cLayerSpec_all,2));
SpecXcorr_all_p = zeros(size(cLayerSpec_all,2),size(cLayerSpec_all,2));
for x = 1:size(cLayerSpec_all,2)
    for y = 1:size(cLayerSpec_all,2)
        [SpecXcorr_all(x,y),SpecXcorr_all_p(x,y)] = corr(cLayerSpec_all(:,x),cLayerSpec_all(:,y),...
            'type','spearman','rows','complete');
    end
end

SpecXcorr_NWh = zeros(size(cLayerSpec_NWh,2),size(cLayerSpec_NWh,2));
SpecXcorr_NWh_p = zeros(size(cLayerSpec_NWh,2),size(cLayerSpec_NWh,2));
for x = 1:size(cLayerSpec_NWh,2)
    for y = 1:size(cLayerSpec_NWh,2)
        [SpecXcorr_NWh(x,y),SpecXcorr_NWh_p(x,y)] = corr(cLayerSpec_NWh(:,x),cLayerSpec_NWh(:,y),...
            'type','spearman','rows','complete');
    end
end

SpecXcorr_Wh = zeros(size(cLayerSpec_Wh,2),size(cLayerSpec_Wh,2));
SpecXcorr_Wh_p = zeros(size(cLayerSpec_Wh,2),size(cLayerSpec_Wh,2));
for x = 1:size(cLayerSpec_Wh,2)
    for y = 1:size(cLayerSpec_Wh,2)
        [SpecXcorr_Wh(x,y),SpecXcorr_Wh_p(x,y)] = corr(cLayerSpec_Wh(:,x),cLayerSpec_Wh(:,y),...
            'type','spearman','rows','complete');
    end
end

SpecXcorr_loP = zeros(size(cLayerSpec_loP,2),size(cLayerSpec_loP,2));
SpecXcorr_loP_p = zeros(size(cLayerSpec_loP,2),size(cLayerSpec_loP,2));
for x = 1:size(cLayerSpec_loP,2)
    for y = 1:size(cLayerSpec_loP,2)
        [SpecXcorr_loP(x,y),SpecXcorr_loP_p(x,y)] = corr(cLayerSpec_loP(:,x),cLayerSpec_loP(:,y),...
            'type','spearman','rows','complete');
    end
end

SpecXcorr_hiP = zeros(size(cLayerSpec_hiP,2),size(cLayerSpec_hiP,2));
SpecXcorr_hiP_p = zeros(size(cLayerSpec_hiP,2),size(cLayerSpec_hiP,2));
for x = 1:size(cLayerSpec_hiP,2)
    for y = 1:size(cLayerSpec_hiP,2)
        [SpecXcorr_hiP(x,y),SpecXcorr_hiP_p(x,y)] = corr(cLayerSpec_hiP(:,x),cLayerSpec_hiP(:,y),...
            'type','spearman','rows','complete');
    end
end

% Saving to struct
LayerSpectral.SpecXcorr_all = SpecXcorr_all;
LayerSpectral.SpecXcorr_all_p = SpecXcorr_all_p;
LayerSpectral.SpecXcorr_NWh = SpecXcorr_NWh;
LayerSpectral.SpecXcorr_NWh_p = SpecXcorr_NWh_p;
LayerSpectral.SpecXcorr_Wh = SpecXcorr_Wh;
LayerSpectral.SpecXcorr_Wh_p = SpecXcorr_Wh_p;
LayerSpectral.SpecXcorr_loP = SpecXcorr_loP;
LayerSpectral.SpecXcorr_loP_p = SpecXcorr_loP_p;
LayerSpectral.SpecXcorr_hiP = SpecXcorr_hiP;
LayerSpectral.SpecXcorr_hiP_p = SpecXcorr_hiP_p;

%% FIGURE:
cmin = min([min(min(SpecXcorr_all(usechannels+1,usechannels+1)))...
    min(min(SpecXcorr_NWh(usechannels+1,usechannels+1)))...
    min(min(SpecXcorr_Wh(usechannels+1,usechannels+1)))]);
cmax = max([max(max(SpecXcorr_all(usechannels+1,usechannels+1)))...
    max(max(SpecXcorr_NWh(usechannels+1,usechannels+1)))...
    max(max(SpecXcorr_Wh(usechannels+1,usechannels+1)))]);

figure;
subplot(1,3,1);
imagesc(normdepth,normdepth,SpecXcorr_all(usechannels+1,usechannels+1));
colormap(gca,'jet')
axis square
ColorbarWithAxis([cmin cmax],['Spearman corr'])
caxis([cmin cmax])
xlabel('normalized depth');ylabel('normalized depth');
title('LFPxcorr ALL');

subplot(1,3,2);
imagesc(normdepth,normdepth,SpecXcorr_NWh(usechannels+1,usechannels+1));
colormap(gca,'jet')
axis square
ColorbarWithAxis([cmin cmax],['Spearman corr'])
caxis([cmin cmax])
xlabel('normalized depth');ylabel('normalized depth');
title('LFPxcorr NWh');

subplot(1,3,3);
imagesc(normdepth,normdepth,SpecXcorr_Wh(usechannels+1,usechannels+1));
colormap(gca,'jet')
axis square
ColorbarWithAxis([cmin cmax],['Spearman corr'])
caxis([cmin cmax])
xlabel('normalized depth');ylabel('normalized depth');
title('LFPxcorr Wh');

NiceSave('LaminarLFPXcorr_Wh_NWh',figfolder,baseName)

%% FIGURE:
cmin = min([min(min(SpecXcorr_all(usechannels+1,usechannels+1)))...
    min(min(SpecXcorr_loP(usechannels+1,usechannels+1)))...
    min(min(SpecXcorr_hiP(usechannels+1,usechannels+1)))]);
cmax = max([max(max(SpecXcorr_all(usechannels+1,usechannels+1)))...
    max(max(SpecXcorr_loP(usechannels+1,usechannels+1)))...
    max(max(SpecXcorr_hiP(usechannels+1,usechannels+1)))]);

figure;
subplot(1,3,1);
imagesc(normdepth,normdepth,SpecXcorr_all(usechannels+1,usechannels+1));
colormap(gca,'jet')
axis square
ColorbarWithAxis([cmin cmax],['Spearman corr'])
caxis([cmin cmax])
xlabel('normalized depth');ylabel('normalized depth');
title('LFPxcorr ALL');

subplot(1,3,2);
imagesc(normdepth,normdepth,SpecXcorr_loP(usechannels+1,usechannels+1));
colormap(gca,'jet')
axis square
ColorbarWithAxis([cmin cmax],['Spearman corr'])
caxis([cmin cmax])
xlabel('normalized depth');ylabel('normalized depth');
title('LFPxcorr <median Pupil area');

subplot(1,3,3);
imagesc(normdepth,normdepth,SpecXcorr_hiP(usechannels+1,usechannels+1));
colormap(gca,'jet')
axis square
ColorbarWithAxis([cmin cmax],['Spearman corr'])
caxis([cmin cmax])
xlabel('normalized depth');ylabel('normalized depth');
title('LFPxcorr >median Pupil area');

NiceSave('LaminarLFPXcorr_lo_hiPup',figfolder,baseName)

%% Laminar eventCSD/MUA
% sort thourgh high/low Wh/pupil/SW amplitude events
% adjust by phase of SlowWaves
% easyyy lines but specifying better intervals, trial selection
% establish and z score to baseline

% eventCSD = eventCSD (lfp, events, varargin);
% eventMUA = eventMUA (mua, events, varargin);
% eventSpec = eventSpec (spec, events, varargin);

%% Laminar CSD-CSD/MUA-MUA cross-corr
% thin about this...

%% Layer-averaged PSpecs, eventSpec, and comodulograms, by state
twin = [0.75 0.75].*wavespec.samplingRate;

dLayerSpec_all = NaN(size(wavespec.data,1),size(wavespec.data,2),6);
dLayerSpec_NWh = NaN(size(wavespec.data(round(allidx_NWh),:),1),...
    size(wavespec.data(round(allidx_NWh),:),2),6); 
dLayerSpec_Wh = NaN(size(wavespec.data(round(allidx_Wh),:),1),...
    size(wavespec.data(round(allidx_Wh),:),2),6);
% dLayerSpec_T = NaN(size(wavespec.data,1),size(wavespec.data,2),6);
dLayerSpec_hiP = NaN(size(wavespec.data(round(eventshipupil),:),1),...
    size(wavespec.data(round(eventshipupil),:),2),6);
dLayerSpec_loP = NaN(size(wavespec.data(round(eventslopupil),:),1),...
    size(wavespec.data(round(eventslopupil),:),2),6);

% Layer 1
L1eventSpec_Wh = NaN(twin(1)+twin(2)+1,size(wavespec.data,2),length(L1idx));
for x = 1: length(L1idx)
    load(fullfile(savefolder,[baseName,'.',num2str(L1idx(x)),'.WaveSpec.lfp.mat']));
    wavespec.data = (log10(abs(wavespec.data))-nanmean(log10(abs(wavespec.data)),1))...
        ./nanstd(log10(abs(wavespec.data)),0,1);
    %wavespec.data = NormToInt(log10(abs(wavespec.data)),'modZ');
    
    tempSpec = cat(3,dLayerSpec_all(:,:,1),wavespec.data);
    dLayerSpec_all(:,:,1) = nansum(tempSpec,3);
    
    tempSpec = cat(3,dLayerSpec_NWh(:,:,1),wavespec.data(round(allidx_NWh),:));
    dLayerSpec_NWh(:,:,1) = nansum(tempSpec,3);
    
    tempSpec = cat(3,dLayerSpec_Wh(:,:,1),wavespec.data(round(allidx_Wh),:));
    dLayerSpec_Wh(:,:,1) = nansum(tempSpec,3);
    
    tempSpec = cat(3,dLayerSpec_loP(:,:,1),wavespec.data(round(eventslopupil),:));
    dLayerSpec_loP(:,:,1) = nansum(tempSpec,3);
    
    tempSpec = cat(3,dLayerSpec_hiP(:,:,1),wavespec.data(round(eventshipupil),:));
    dLayerSpec_hiP(:,:,1) = nansum(tempSpec,3);
    
    % eventSpec
    events = round(events_Wh.*wavespec.samplingRate);
    spec_temp = NaN(twin(1)+twin(2)+1,size(wavespec.data,2),size(events,1));
    for e = 1:size(events,1)
        if events(e,1)-twin(1) > 0 && events(e,1)+twin(2) < size(wavespec.data,1)
            spec_temp(:,:,e) = wavespec.data(events(e,1)-twin(1):events(e,1)+twin(2),:);
        else
        end
    end
    L1eventSpec_Wh(:,:,x) = nanmean(spec_temp,3);
end
L1eventSpec_Wh = squeeze(nanmean(L1eventSpec_Wh,3));
dLayerSpec_all(:,:,1) = dLayerSpec_all(:,:,1)./length(L1idx);
dLayerSpec_Wh(:,:,1) = dLayerSpec_Wh(:,:,1)./length(L1idx);
dLayerSpec_NWh(:,:,1) = dLayerSpec_NWh(:,:,1)./length(L1idx);
dLayerSpec_loP(:,:,1) = dLayerSpec_loP(:,:,1)./length(L1idx);
dLayerSpec_hiP(:,:,1) = dLayerSpec_hiP(:,:,1)./length(L1idx);

% Layer 2/3
L23eventSpec_Wh = NaN(twin(1)+twin(2)+1,size(wavespec.data,2),length(L23idx));
for x = 1: length(L23idx)
    load(fullfile(savefolder,[baseName,'.',num2str(L23idx(x)),'.WaveSpec.lfp.mat']));
    wavespec.data = (log10(abs(wavespec.data))-nanmean(log10(abs(wavespec.data)),1))...
        ./nanstd(log10(abs(wavespec.data)),0,1);
    %wavespec.data = NormToInt(log10(abs(wavespec.data)),'modZ');

    tempSpec = cat(3,dLayerSpec_all(:,:,2),wavespec.data);
    dLayerSpec_all(:,:,2) = nansum(tempSpec,3);
    
    tempSpec = cat(3,dLayerSpec_NWh(:,:,2),wavespec.data(round(allidx_NWh),:));
    dLayerSpec_NWh(:,:,2) = nansum(tempSpec,3);
    
    tempSpec = cat(3,dLayerSpec_Wh(:,:,2),wavespec.data(round(allidx_Wh),:));
    dLayerSpec_Wh(:,:,2) = nansum(tempSpec,3);
    
    tempSpec = cat(3,dLayerSpec_loP(:,:,2),wavespec.data(round(eventslopupil),:));
    dLayerSpec_loP(:,:,2) = nansum(tempSpec,3);
    
    tempSpec = cat(3,dLayerSpec_hiP(:,:,2),wavespec.data(round(eventshipupil),:));
    dLayerSpec_hiP(:,:,2) = nansum(tempSpec,3);
    
    % eventSpec
    events = round(events_Wh.*wavespec.samplingRate);
    spec_temp = NaN(twin(1)+twin(2)+1,size(wavespec.data,2),size(events,1));
    for e = 1:size(events,1)
        if events(e,1)-twin(1) > 0 && events(e,1)+twin(2) < size(wavespec.data,1)
            spec_temp(:,:,e) = wavespec.data(events(e,1)-twin(1):events(e,1)+twin(2),:);
        else
        end
    end
    L23eventSpec_Wh(:,:,x) = nanmean(spec_temp,3);
end
L23eventSpec_Wh = squeeze(nanmean(L23eventSpec_Wh,3));
dLayerSpec_all(:,:,2) = dLayerSpec_all(:,:,2)./length(L23idx);
dLayerSpec_Wh(:,:,2) = dLayerSpec_Wh(:,:,2)./length(L23idx);
dLayerSpec_NWh(:,:,2) = dLayerSpec_NWh(:,:,2)./length(L23idx);
dLayerSpec_loP(:,:,2) = dLayerSpec_loP(:,:,2)./length(L23idx);
dLayerSpec_hiP(:,:,2) = dLayerSpec_hiP(:,:,2)./length(L23idx);

% Layer 4
L4eventSpec_Wh = NaN(twin(1)+twin(2)+1,size(wavespec.data,2),length(L4idx));
for x = 1: length(L4idx)
    load(fullfile(savefolder,[baseName,'.',num2str(L4idx(x)),'.WaveSpec.lfp.mat']));
    wavespec.data = (log10(abs(wavespec.data))-nanmean(log10(abs(wavespec.data)),1))...
        ./nanstd(log10(abs(wavespec.data)),0,1);
    %wavespec.data = NormToInt(log10(abs(wavespec.data)),'modZ');

    tempSpec = cat(3,dLayerSpec_all(:,:,3),wavespec.data);
    dLayerSpec_all(:,:,3) = nansum(tempSpec,3);
    
    tempSpec = cat(3,dLayerSpec_NWh(:,:,3),wavespec.data(round(allidx_NWh),:));
    dLayerSpec_NWh(:,:,3) = nansum(tempSpec,3);
    
    tempSpec = cat(3,dLayerSpec_Wh(:,:,3),wavespec.data(round(allidx_Wh),:));
    dLayerSpec_Wh(:,:,3) = nansum(tempSpec,3);
    
    tempSpec = cat(3,dLayerSpec_loP(:,:,3),wavespec.data(round(eventslopupil),:));
    dLayerSpec_loP(:,:,3) = nansum(tempSpec,3);
    
    tempSpec = cat(3,dLayerSpec_hiP(:,:,3),wavespec.data(round(eventshipupil),:));
    dLayerSpec_hiP(:,:,3) = nansum(tempSpec,3);
    
     % eventSpec
    events = round(events_Wh.*wavespec.samplingRate);
    spec_temp = NaN(twin(1)+twin(2)+1,size(wavespec.data,2),size(events,1));
    for e = 1:size(events,1)
        if events(e,1)-twin(1) > 0 && events(e,1)+twin(2) < size(wavespec.data,1)
            spec_temp(:,:,e) = wavespec.data(events(e,1)-twin(1):events(e,1)+twin(2),:);
        else
        end
    end
    L4eventSpec_Wh(:,:,x) = nanmean(spec_temp,3);
end
L4eventSpec_Wh = squeeze(nanmean(L4eventSpec_Wh,3));
dLayerSpec_all(:,:,3) = dLayerSpec_all(:,:,3)./length(L4idx);
dLayerSpec_Wh(:,:,3) = dLayerSpec_Wh(:,:,3)./length(L4idx);
dLayerSpec_NWh(:,:,3) = dLayerSpec_NWh(:,:,3)./length(L4idx);
dLayerSpec_loP(:,:,3) = dLayerSpec_loP(:,:,3)./length(L4idx);
dLayerSpec_hiP(:,:,3) = dLayerSpec_hiP(:,:,3)./length(L4idx);

% Layer 5a
L5aeventSpec_Wh = NaN(twin(1)+twin(2)+1,size(wavespec.data,2),length(L5aidx));
for x = 1: length(L5aidx)
    load(fullfile(savefolder,[baseName,'.',num2str(L5aidx(x)),'.WaveSpec.lfp.mat']));
    wavespec.data = (log10(abs(wavespec.data))-nanmean(log10(abs(wavespec.data)),1))...
        ./nanstd(log10(abs(wavespec.data)),0,1);
    %wavespec.data = NormToInt(log10(abs(wavespec.data)),'modZ');

    tempSpec = cat(3,dLayerSpec_all(:,:,4),wavespec.data);
    dLayerSpec_all(:,:,4) = nansum(tempSpec,3);
    
    tempSpec = cat(3,dLayerSpec_NWh(:,:,4),wavespec.data(round(allidx_NWh),:));
    dLayerSpec_NWh(:,:,4) = nansum(tempSpec,3);
    
    tempSpec = cat(3,dLayerSpec_Wh(:,:,4),wavespec.data(round(allidx_Wh),:));
    dLayerSpec_Wh(:,:,4) = nansum(tempSpec,3);
    
    tempSpec = cat(3,dLayerSpec_loP(:,:,4),wavespec.data(round(eventslopupil),:));
    dLayerSpec_loP(:,:,4) = nansum(tempSpec,3);
    
    tempSpec = cat(3,dLayerSpec_hiP(:,:,4),wavespec.data(round(eventshipupil),:));
    dLayerSpec_hiP(:,:,4) = nansum(tempSpec,3);
    
    % eventSpec
    events = round(events_Wh.*wavespec.samplingRate);
    spec_temp = NaN(twin(1)+twin(2)+1,size(wavespec.data,2),size(events,1));
    for e = 1:size(events,1)
        if events(e,1)-twin(1) > 0 && events(e,1)+twin(2) < size(wavespec.data,1)
            spec_temp(:,:,e) = wavespec.data(events(e,1)-twin(1):events(e,1)+twin(2),:);
        else
        end
    end
    L5aeventSpec_Wh(:,:,x) = nanmean(spec_temp,3);
end
L5aeventSpec_Wh = squeeze(nanmean(L5aeventSpec_Wh,3));
dLayerSpec_all(:,:,4) = dLayerSpec_all(:,:,4)./length(L5aidx);
dLayerSpec_Wh(:,:,4) = dLayerSpec_Wh(:,:,4)./length(L5aidx);
dLayerSpec_NWh(:,:,4) = dLayerSpec_NWh(:,:,4)./length(L5aidx);
dLayerSpec_loP(:,:,4) = dLayerSpec_loP(:,:,4)./length(L5aidx);
dLayerSpec_hiP(:,:,4) = dLayerSpec_hiP(:,:,4)./length(L5aidx);

% Layer 5b/6
L56eventSpec_Wh = NaN(twin(1)+twin(2)+1,size(wavespec.data,2),length(L56idx));
for x = 1: length(L56idx)
    load(fullfile(savefolder,[baseName,'.',num2str(L56idx(x)),'.WaveSpec.lfp.mat']));
    wavespec.data = (log10(abs(wavespec.data))-nanmean(log10(abs(wavespec.data)),1))...
        ./nanstd(log10(abs(wavespec.data)),0,1);
    %wavespec.data = NormToInt(log10(abs(wavespec.data)),'modZ');

    tempSpec = cat(3,dLayerSpec_all(:,:,5),wavespec.data);
    dLayerSpec_all(:,:,5) = nansum(tempSpec,3);
    
    tempSpec = cat(3,dLayerSpec_NWh(:,:,5),wavespec.data(round(allidx_NWh),:));
    dLayerSpec_NWh(:,:,5) = nansum(tempSpec,3);
    
    tempSpec = cat(3,dLayerSpec_Wh(:,:,5),wavespec.data(round(allidx_Wh),:));
    dLayerSpec_Wh(:,:,5) = nansum(tempSpec,3);
    
    tempSpec = cat(3,dLayerSpec_loP(:,:,5),wavespec.data(round(eventslopupil),:));
    dLayerSpec_loP(:,:,5) = nansum(tempSpec,3);
    
    tempSpec = cat(3,dLayerSpec_hiP(:,:,5),wavespec.data(round(eventshipupil),:));
    dLayerSpec_hiP(:,:,5) = nansum(tempSpec,3);
    
    % eventSpec
    events = round(events_Wh.*wavespec.samplingRate);
    spec_temp = NaN(twin(1)+twin(2)+1,size(wavespec.data,2),size(events,1));
    for e = 1:size(events,1)
        if events(e,1)-twin(1) > 0 && events(e,1)+twin(2) < size(wavespec.data,1)
            spec_temp(:,:,e) = wavespec.data(events(e,1)-twin(1):events(e,1)+twin(2),:);
        else
        end
    end
    L56eventSpec_Wh(:,:,x) = nanmean(spec_temp,3);
end
L56eventSpec_Wh = squeeze(nanmean(L56eventSpec_Wh,3));
dLayerSpec_all(:,:,5) = dLayerSpec_all(:,:,5)./length(L56idx);
dLayerSpec_Wh(:,:,5) = dLayerSpec_Wh(:,:,5)./length(L56idx);
dLayerSpec_NWh(:,:,5) = dLayerSpec_NWh(:,:,5)./length(L56idx);
dLayerSpec_loP(:,:,5) = dLayerSpec_loP(:,:,5)./length(L56idx);
dLayerSpec_hiP(:,:,5) = dLayerSpec_hiP(:,:,5)./length(L56idx);

% Layer 6
L6eventSpec_Wh = NaN(twin(1)+twin(2)+1,size(wavespec.data,2),length(L6idx));
for x = 1: length(L6idx)
    load(fullfile(savefolder,[baseName,'.',num2str(L6idx(x)),'.WaveSpec.lfp.mat']));
    wavespec.data = (log10(abs(wavespec.data))-nanmean(log10(abs(wavespec.data)),1))...
        ./nanstd(log10(abs(wavespec.data)),0,1);
    %wavespec.data = NormToInt(log10(abs(wavespec.data)),'modZ');
    
    tempSpec = cat(3,dLayerSpec_all(:,:,6),wavespec.data);
    dLayerSpec_all(:,:,6) = nansum(tempSpec,3);
    
    tempSpec = cat(3,dLayerSpec_NWh(:,:,6),wavespec.data(round(allidx_NWh),:));
    dLayerSpec_NWh(:,:,6) = nansum(tempSpec,3);
    
    tempSpec = cat(3,dLayerSpec_Wh(:,:,6),wavespec.data(round(allidx_Wh),:));
    dLayerSpec_Wh(:,:,6) = nansum(tempSpec,3);
    
    tempSpec = cat(3,dLayerSpec_loP(:,:,6),wavespec.data(round(eventslopupil),:));
    dLayerSpec_loP(:,:,6) = nansum(tempSpec,3);
    
    tempSpec = cat(3,dLayerSpec_hiP(:,:,6),wavespec.data(round(eventshipupil),:));
    dLayerSpec_hiP(:,:,6) = nansum(tempSpec,3);
    
    % eventSpec
    events = round(events_Wh.*wavespec.samplingRate);
    spec_temp = NaN(twin(1)+twin(2)+1,size(wavespec.data,2),size(events,1));
    for e = 1:size(events,1)
        if events(e,1)-twin(1) > 0 && events(e,1)+twin(2) < size(wavespec.data,1)
            spec_temp(:,:,e) = wavespec.data(events(e,1)-twin(1):events(e,1)+twin(2),:);
        else
        end
    end
    L6eventSpec_Wh(:,:,x) = nanmean(spec_temp,3);
end
L6eventSpec_Wh = squeeze(nanmean(L6eventSpec_Wh,3));
dLayerSpec_all(:,:,6) = dLayerSpec_all(:,:,6)./length(L6idx);
dLayerSpec_Wh(:,:,6) = dLayerSpec_Wh(:,:,6)./length(L6idx);
dLayerSpec_NWh(:,:,6) = dLayerSpec_NWh(:,:,6)./length(L6idx);
dLayerSpec_loP(:,:,6) = dLayerSpec_loP(:,:,6)./length(L6idx);
dLayerSpec_hiP(:,:,6) = dLayerSpec_hiP(:,:,6)./length(L6idx);

% Saving to struct
LayerSpectral.L1chan = length(L1idx);
LayerSpectral.L1eventSpec_Wh = L1eventSpec_Wh;
LayerSpectral.L1Spec_all = squeeze(nanmean(dLayerSpec_all(:,:,1),1)); 
LayerSpectral.L1Spec_Wh = squeeze(nanmean(dLayerSpec_Wh(:,:,1),1)); 
LayerSpectral.L1Spec_NWh = squeeze(nanmean(dLayerSpec_NWh(:,:,1),1)); 
LayerSpectral.L1Spec_loP = squeeze(nanmean(dLayerSpec_loP(:,:,1),1));
LayerSpectral.L1Spec_hiP = squeeze(nanmean(dLayerSpec_hiP(:,:,1),1));

LayerSpectral.L23chan = length(L23idx);
LayerSpectral.L23eventSpec_Wh = L23eventSpec_Wh;
LayerSpectral.L23Spec_all = squeeze(nanmean(dLayerSpec_all(:,:,2),1)); 
LayerSpectral.L23Spec_Wh = squeeze(nanmean(dLayerSpec_Wh(:,:,2),1)); 
LayerSpectral.L23Spec_NWh = squeeze(nanmean(dLayerSpec_NWh(:,:,2),1)); 
LayerSpectral.L23Spec_loP = squeeze(nanmean(dLayerSpec_loP(:,:,2),1));
LayerSpectral.L23Spec_hiP = squeeze(nanmean(dLayerSpec_hiP(:,:,2),1));

LayerSpectral.L4chan = length(L4idx);
LayerSpectral.L4eventSpec_Wh = L4eventSpec_Wh;
LayerSpectral.L4Spec_all = squeeze(nanmean(dLayerSpec_all(:,:,3),1)); 
LayerSpectral.L4Spec_Wh = squeeze(nanmean(dLayerSpec_Wh(:,:,3),1)); 
LayerSpectral.L4Spec_NWh = squeeze(nanmean(dLayerSpec_NWh(:,:,3),1)); 
LayerSpectral.L4Spec_loP = squeeze(nanmean(dLayerSpec_loP(:,:,3),1));
LayerSpectral.L4Spec_hiP = squeeze(nanmean(dLayerSpec_hiP(:,:,3),1));

LayerSpectral.L5achan = length(L5aidx);
LayerSpectral.L5aeventSpec_Wh = L5aeventSpec_Wh;
LayerSpectral.L5aSpec_all = squeeze(nanmean(dLayerSpec_all(:,:,4),1)); 
LayerSpectral.L5aSpec_Wh = squeeze(nanmean(dLayerSpec_Wh(:,:,4),1)); 
LayerSpectral.L5aSpec_NWh = squeeze(nanmean(dLayerSpec_NWh(:,:,4),1)); 
LayerSpectral.L5aSpec_loP = squeeze(nanmean(dLayerSpec_loP(:,:,4),1));
LayerSpectral.L5aSpec_hiP = squeeze(nanmean(dLayerSpec_hiP(:,:,4),1));

LayerSpectral.L56chan = length(L56idx);
LayerSpectral.L56eventSpec_Wh = L56eventSpec_Wh;
LayerSpectral.L56Spec_all = squeeze(nanmean(dLayerSpec_all(:,:,5),1)); 
LayerSpectral.L56Spec_Wh = squeeze(nanmean(dLayerSpec_Wh(:,:,5),1)); 
LayerSpectral.L56Spec_NWh = squeeze(nanmean(dLayerSpec_NWh(:,:,5),1)); 
LayerSpectral.L56Spec_loP = squeeze(nanmean(dLayerSpec_loP(:,:,5),1));
LayerSpectral.L56Spec_hiP = squeeze(nanmean(dLayerSpec_hiP(:,:,5),1));

LayerSpectral.L6chan = length(L6idx);
LayerSpectral.L6eventSpec_Wh = L6eventSpec_Wh;
LayerSpectral.L6Spec_all = squeeze(nanmean(dLayerSpec_all(:,:,6),1)); 
LayerSpectral.L6Spec_Wh = squeeze(nanmean(dLayerSpec_Wh(:,:,6),1)); 
LayerSpectral.L6Spec_NWh = squeeze(nanmean(dLayerSpec_NWh(:,:,6),1)); 
LayerSpectral.L6Spec_loP = squeeze(nanmean(dLayerSpec_loP(:,:,6),1));
LayerSpectral.L6Spec_hiP = squeeze(nanmean(dLayerSpec_hiP(:,:,6),1));

%% FIGURE:
taxis = (-(twin(1)/wavespec.samplingRate):(1/wavespec.samplingRate):(twin(2)/wavespec.samplingRate))*1e3;
cmax = max([max(max(L1eventSpec_Wh)) max(max(L23eventSpec_Wh))...
    max(max(L4eventSpec_Wh)) max(max(L5aeventSpec_Wh))...
    max(max(L56eventSpec_Wh)) max(max(L6eventSpec_Wh))]);

figure;
subplot(3,2,1);
imagesc(taxis,log10(wavespec.freqs),L1eventSpec_Wh');hold on;
colormap jet;
LogScale('y',10);
caxis([-cmax cmax]);
axis xy
xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.freqs(1)) log10(wavespec.freqs(end))],'--k');hold on;
title('L1');
    
subplot(3,2,3);
imagesc(taxis,log10(wavespec.freqs),L23eventSpec_Wh');hold on;
colormap jet;
LogScale('y',10);
caxis([-cmax cmax]);
axis xy
xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.freqs(1)) log10(wavespec.freqs(end))],'--k');hold on;
title('L2/3');
    
subplot(3,2,5);
imagesc(taxis,log10(wavespec.freqs),L4eventSpec_Wh');hold on;
colormap jet;
LogScale('y',10);
caxis([-cmax cmax]);
axis xy
xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.freqs(1)) log10(wavespec.freqs(end))],'--k');hold on;
title('L4')

subplot(3,2,2);
imagesc(taxis,log10(wavespec.freqs),L5aeventSpec_Wh');hold on;
colormap jet;
LogScale('y',10);
caxis([-cmax cmax]);
axis xy
xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.freqs(1)) log10(wavespec.freqs(end))],'--k');hold on;
title('L5a')
    
subplot(3,2,4);
imagesc(taxis,log10(wavespec.freqs),L56eventSpec_Wh');hold on;
colormap jet;
LogScale('y',10);
caxis([-cmax cmax]);
axis xy
xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.freqs(1)) log10(wavespec.freqs(end))],'--k');hold on;
title('L5/6')
    
subplot(3,2,6);
imagesc(taxis,log10(wavespec.freqs),L6eventSpec_Wh');hold on;
colormap jet;
LogScale('y',10);
caxis([-cmax cmax]);
axis xy
xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.freqs(1)) log10(wavespec.freqs(end))],'--k');hold on;
title('L6')
  
NiceSave('Laminar_eventSpec_Wh',figfolder,baseName)
    
%% State comodulograms
i = 1;
L1comodcorrs_all = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
L1comodcorrs_NWh = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
L1comodcorrs_Wh = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
L1comodcorrs_loP = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
L1comodcorrs_hiP = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
for ii = 1:size(dLayerSpec_all,3)
    L1comodcorrs_all(:,:,ii) = corr(dLayerSpec_all(:,:,i),...
        dLayerSpec_all(:,:,ii),'type','spearman',...
        'rows','complete');
    L1comodcorrs_NWh(:,:,ii) = corr(dLayerSpec_NWh(:,:,i),...
        dLayerSpec_NWh(:,:,ii),'type','spearman',...
        'rows','complete');
    L1comodcorrs_Wh(:,:,ii) = corr(dLayerSpec_Wh(:,:,i),...
        dLayerSpec_Wh(:,:,ii),'type','spearman',...
        'rows','complete');
    L1comodcorrs_loP(:,:,ii) = corr(dLayerSpec_loP(:,:,i),...
        dLayerSpec_loP(:,:,ii),'type','spearman',...
        'rows','complete');
    L1comodcorrs_hiP(:,:,ii) = corr(dLayerSpec_hiP(:,:,i),...
        dLayerSpec_hiP(:,:,ii),'type','spearman',...
        'rows','complete');
end

i = 2;
L23comodcorrs_all = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
L23comodcorrs_NWh = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
L23comodcorrs_Wh = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
L23comodcorrs_loP = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
L23comodcorrs_hiP = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
for ii = 1:size(dLayerSpec_all,3)
    L23comodcorrs_all(:,:,ii) = corr(dLayerSpec_all(:,:,i),...
        dLayerSpec_all(:,:,ii),'type','spearman',...
        'rows','complete');
    L23comodcorrs_NWh(:,:,ii) = corr(dLayerSpec_NWh(:,:,i),...
        dLayerSpec_NWh(:,:,ii),'type','spearman',...
        'rows','complete');
    L23comodcorrs_Wh(:,:,ii) = corr(dLayerSpec_Wh(:,:,i),...
        dLayerSpec_Wh(:,:,ii),'type','spearman',...
        'rows','complete');
    L23comodcorrs_loP(:,:,ii) = corr(dLayerSpec_loP(:,:,i),...
        dLayerSpec_loP(:,:,ii),'type','spearman',...
        'rows','complete');
    L23comodcorrs_hiP(:,:,ii) = corr(dLayerSpec_hiP(:,:,i),...
        dLayerSpec_hiP(:,:,ii),'type','spearman',...
        'rows','complete');
end

i = 3;
L4comodcorrs_all = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
L4comodcorrs_NWh = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
L4comodcorrs_Wh = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
L4comodcorrs_loP = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
L4comodcorrs_hiP = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
for ii = 1:size(dLayerSpec_all,3)
    L4comodcorrs_all(:,:,ii) = corr(dLayerSpec_all(:,:,i),...
        dLayerSpec_all(:,:,ii),'type','spearman',...
        'rows','complete');
    L4comodcorrs_NWh(:,:,ii) = corr(dLayerSpec_NWh(:,:,i),...
        dLayerSpec_NWh(:,:,ii),'type','spearman',...
        'rows','complete');
    L4comodcorrs_Wh(:,:,ii) = corr(dLayerSpec_Wh(:,:,i),...
        dLayerSpec_Wh(:,:,ii),'type','spearman',...
        'rows','complete');
    L4comodcorrs_loP(:,:,ii) = corr(dLayerSpec_loP(:,:,i),...
        dLayerSpec_loP(:,:,ii),'type','spearman',...
        'rows','complete');
    L4comodcorrs_hiP(:,:,ii) = corr(dLayerSpec_hiP(:,:,i),...
        dLayerSpec_hiP(:,:,ii),'type','spearman',...
        'rows','complete');
end

i = 4;
L5acomodcorrs_all = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
L5acomodcorrs_NWh = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
L5acomodcorrs_Wh = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
L5acomodcorrs_loP = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
L5acomodcorrs_hiP = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
for ii = 1:size(dLayerSpec_all,3)
    L5acomodcorrs_all(:,:,ii) = corr(dLayerSpec_all(:,:,i),...
        dLayerSpec_all(:,:,ii),'type','spearman',...
        'rows','complete');
    L5acomodcorrs_NWh(:,:,ii) = corr(dLayerSpec_NWh(:,:,i),...
        dLayerSpec_NWh(:,:,ii),'type','spearman',...
        'rows','complete');
    L5acomodcorrs_Wh(:,:,ii) = corr(dLayerSpec_Wh(:,:,i),...
        dLayerSpec_Wh(:,:,ii),'type','spearman',...
        'rows','complete');
    L5acomodcorrs_loP(:,:,ii) = corr(dLayerSpec_loP(:,:,i),...
        dLayerSpec_loP(:,:,ii),'type','spearman',...
        'rows','complete');
    L5acomodcorrs_hiP(:,:,ii) = corr(dLayerSpec_hiP(:,:,i),...
        dLayerSpec_hiP(:,:,ii),'type','spearman',...
        'rows','complete');
end

i = 5;
L56comodcorrs_all = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
L56comodcorrs_NWh = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
L56comodcorrs_Wh = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
L56comodcorrs_loP = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
L56comodcorrs_hiP = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
for ii = 1:size(dLayerSpec_all,3)
   L56comodcorrs_all(:,:,ii) = corr(dLayerSpec_all(:,:,i),...
        dLayerSpec_all(:,:,ii),'type','spearman',...
        'rows','complete');
    L56comodcorrs_NWh(:,:,ii) = corr(dLayerSpec_NWh(:,:,i),...
        dLayerSpec_NWh(:,:,ii),'type','spearman',...
        'rows','complete');
    L56comodcorrs_Wh(:,:,ii) = corr(dLayerSpec_Wh(:,:,i),...
        dLayerSpec_Wh(:,:,ii),'type','spearman',...
        'rows','complete');
    L56comodcorrs_loP(:,:,ii) = corr(dLayerSpec_loP(:,:,i),...
        dLayerSpec_loP(:,:,ii),'type','spearman',...
        'rows','complete');
    L56comodcorrs_hiP(:,:,ii) = corr(dLayerSpec_hiP(:,:,i),...
        dLayerSpec_hiP(:,:,ii),'type','spearman',...
        'rows','complete');
end

i = 6;
L6comodcorrs_all = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
L6comodcorrs_NWh = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
L6comodcorrs_Wh = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
L6comodcorrs_loP = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
L6comodcorrs_hiP = NaN(size(dLayerSpec_all,2),size(dLayerSpec_all,2),6);
for ii = 1:size(dLayerSpec_all,3)
    L6comodcorrs_all(:,:,ii) = corr(dLayerSpec_all(:,:,i),...
        dLayerSpec_all(:,:,ii),'type','spearman',...
        'rows','complete');
    L6comodcorrs_NWh(:,:,ii) = corr(dLayerSpec_NWh(:,:,i),...
        dLayerSpec_NWh(:,:,ii),'type','spearman',...
        'rows','complete');
    L6comodcorrs_Wh(:,:,ii) = corr(dLayerSpec_Wh(:,:,i),...
        dLayerSpec_Wh(:,:,ii),'type','spearman',...
        'rows','complete');
    L6comodcorrs_loP(:,:,ii) = corr(dLayerSpec_loP(:,:,i),...
        dLayerSpec_loP(:,:,ii),'type','spearman',...
        'rows','complete');
    L6comodcorrs_hiP(:,:,ii) = corr(dLayerSpec_hiP(:,:,i),...
        dLayerSpec_hiP(:,:,ii),'type','spearman',...
        'rows','complete');
end

% Saving to struct
LayerSpectral.L1comodcorrs_all.L1 = squeeze(L1comodcorrs_all(:,:,1));
LayerSpectral.L1comodcorrs_all.L23 = squeeze(L1comodcorrs_all(:,:,2));
LayerSpectral.L1comodcorrs_all.L4 = squeeze(L1comodcorrs_all(:,:,3));
LayerSpectral.L1comodcorrs_all.L5a = squeeze(L1comodcorrs_all(:,:,4));
LayerSpectral.L1comodcorrs_all.L56 = squeeze(L1comodcorrs_all(:,:,5));
LayerSpectral.L1comodcorrs_all.L6 = squeeze(L1comodcorrs_all(:,:,6));

LayerSpectral.L1comodcorrs_NWh.L1 = squeeze(L1comodcorrs_NWh(:,:,1));
LayerSpectral.L1comodcorrs_NWh.L23 = squeeze(L1comodcorrs_NWh(:,:,2));
LayerSpectral.L1comodcorrs_NWh.L4 = squeeze(L1comodcorrs_NWh(:,:,3));
LayerSpectral.L1comodcorrs_NWh.L5a = squeeze(L1comodcorrs_NWh(:,:,4));
LayerSpectral.L1comodcorrs_NWh.L56 = squeeze(L1comodcorrs_NWh(:,:,5));
LayerSpectral.L1comodcorrs_NWh.L6 = squeeze(L1comodcorrs_NWh(:,:,6));

LayerSpectral.L1comodcorrs_Wh.L1 = squeeze(L1comodcorrs_Wh(:,:,1));
LayerSpectral.L1comodcorrs_Wh.L23 = squeeze(L1comodcorrs_Wh(:,:,2));
LayerSpectral.L1comodcorrs_Wh.L4 = squeeze(L1comodcorrs_Wh(:,:,3));
LayerSpectral.L1comodcorrs_Wh.L5a = squeeze(L1comodcorrs_Wh(:,:,4));
LayerSpectral.L1comodcorrs_Wh.L56 = squeeze(L1comodcorrs_Wh(:,:,5));
LayerSpectral.L1comodcorrs_Wh.L6 = squeeze(L1comodcorrs_Wh(:,:,6));

LayerSpectral.L1comodcorrs_loP.L1 = squeeze(L1comodcorrs_loP(:,:,1));
LayerSpectral.L1comodcorrs_loP.L23 = squeeze(L1comodcorrs_loP(:,:,2));
LayerSpectral.L1comodcorrs_loP.L4 = squeeze(L1comodcorrs_loP(:,:,3));
LayerSpectral.L1comodcorrs_loP.L5a = squeeze(L1comodcorrs_loP(:,:,4));
LayerSpectral.L1comodcorrs_loP.L56 = squeeze(L1comodcorrs_loP(:,:,5));
LayerSpectral.L1comodcorrs_loP.L6 = squeeze(L1comodcorrs_loP(:,:,6));

LayerSpectral.L1comodcorrs_hiP.L1 = squeeze(L1comodcorrs_hiP(:,:,1));
LayerSpectral.L1comodcorrs_hiP.L23 = squeeze(L1comodcorrs_hiP(:,:,2));
LayerSpectral.L1comodcorrs_hiP.L4 = squeeze(L1comodcorrs_hiP(:,:,3));
LayerSpectral.L1comodcorrs_hiP.L5a = squeeze(L1comodcorrs_hiP(:,:,4));
LayerSpectral.L1comodcorrs_hiP.L56 = squeeze(L1comodcorrs_hiP(:,:,5));
LayerSpectral.L1comodcorrs_hiP.L6 = squeeze(L1comodcorrs_hiP(:,:,6));

%
LayerSpectral.L23comodcorrs_all.L1 = squeeze(L23comodcorrs_all(:,:,1));
LayerSpectral.L23comodcorrs_all.L23 = squeeze(L23comodcorrs_all(:,:,2));
LayerSpectral.L23comodcorrs_all.L4 = squeeze(L23comodcorrs_all(:,:,3));
LayerSpectral.L23comodcorrs_all.L5a = squeeze(L23comodcorrs_all(:,:,4));
LayerSpectral.L23comodcorrs_all.L56 = squeeze(L23comodcorrs_all(:,:,5));
LayerSpectral.L23comodcorrs_all.L6 = squeeze(L23comodcorrs_all(:,:,6));

LayerSpectral.L23comodcorrs_NWh.L1 = squeeze(L23comodcorrs_NWh(:,:,1));
LayerSpectral.L23comodcorrs_NWh.L23 = squeeze(L23comodcorrs_NWh(:,:,2));
LayerSpectral.L23comodcorrs_NWh.L4 = squeeze(L23comodcorrs_NWh(:,:,3));
LayerSpectral.L23comodcorrs_NWh.L5a = squeeze(L23comodcorrs_NWh(:,:,4));
LayerSpectral.L23comodcorrs_NWh.L56 = squeeze(L23comodcorrs_NWh(:,:,5));
LayerSpectral.L23comodcorrs_NWh.L6 = squeeze(L23comodcorrs_NWh(:,:,6));

LayerSpectral.L23comodcorrs_Wh.L1 = squeeze(L23comodcorrs_Wh(:,:,1));
LayerSpectral.L23comodcorrs_Wh.L23 = squeeze(L23comodcorrs_Wh(:,:,2));
LayerSpectral.L23comodcorrs_Wh.L4 = squeeze(L23comodcorrs_Wh(:,:,3));
LayerSpectral.L23comodcorrs_Wh.L5a = squeeze(L23comodcorrs_Wh(:,:,4));
LayerSpectral.L23comodcorrs_Wh.L56 = squeeze(L23comodcorrs_Wh(:,:,5));
LayerSpectral.L23comodcorrs_Wh.L6 = squeeze(L23comodcorrs_Wh(:,:,6));

LayerSpectral.L23comodcorrs_loP.L1 = squeeze(L23comodcorrs_loP(:,:,1));
LayerSpectral.L23comodcorrs_loP.L23 = squeeze(L23comodcorrs_loP(:,:,2));
LayerSpectral.L23comodcorrs_loP.L4 = squeeze(L23comodcorrs_loP(:,:,3));
LayerSpectral.L23comodcorrs_loP.L5a = squeeze(L23comodcorrs_loP(:,:,4));
LayerSpectral.L23comodcorrs_loP.L56 = squeeze(L23comodcorrs_loP(:,:,5));
LayerSpectral.L23comodcorrs_loP.L6 = squeeze(L23comodcorrs_loP(:,:,6));

LayerSpectral.L23comodcorrs_hiP.L1 = squeeze(L23comodcorrs_hiP(:,:,1));
LayerSpectral.L23comodcorrs_hiP.L23 = squeeze(L23comodcorrs_hiP(:,:,2));
LayerSpectral.L23comodcorrs_hiP.L4 = squeeze(L23comodcorrs_hiP(:,:,3));
LayerSpectral.L23comodcorrs_hiP.L5a = squeeze(L23comodcorrs_hiP(:,:,4));
LayerSpectral.L23comodcorrs_hiP.L56 = squeeze(L23comodcorrs_hiP(:,:,5));
LayerSpectral.L23comodcorrs_hiP.L6 = squeeze(L23comodcorrs_hiP(:,:,6));

%
LayerSpectral.L4comodcorrs_all.L1 = squeeze(L4comodcorrs_all(:,:,1));
LayerSpectral.L4comodcorrs_all.L23 = squeeze(L4comodcorrs_all(:,:,2));
LayerSpectral.L4comodcorrs_all.L4 = squeeze(L4comodcorrs_all(:,:,3));
LayerSpectral.L4comodcorrs_all.L5a = squeeze(L4comodcorrs_all(:,:,4));
LayerSpectral.L4comodcorrs_all.L56 = squeeze(L4comodcorrs_all(:,:,5));
LayerSpectral.L4comodcorrs_all.L6 = squeeze(L4comodcorrs_all(:,:,6));

LayerSpectral.L4comodcorrs_NWh.L1 = squeeze(L4comodcorrs_NWh(:,:,1));
LayerSpectral.L4comodcorrs_NWh.L23 = squeeze(L4comodcorrs_NWh(:,:,2));
LayerSpectral.L4comodcorrs_NWh.L4 = squeeze(L4comodcorrs_NWh(:,:,3));
LayerSpectral.L4comodcorrs_NWh.L5a = squeeze(L4comodcorrs_NWh(:,:,4));
LayerSpectral.L4comodcorrs_NWh.L56 = squeeze(L4comodcorrs_NWh(:,:,5));
LayerSpectral.L4comodcorrs_NWh.L6 = squeeze(L4comodcorrs_NWh(:,:,6));

LayerSpectral.L4comodcorrs_Wh.L1 = squeeze(L4comodcorrs_Wh(:,:,1));
LayerSpectral.L4comodcorrs_Wh.L23 = squeeze(L4comodcorrs_Wh(:,:,2));
LayerSpectral.L4comodcorrs_Wh.L4 = squeeze(L4comodcorrs_Wh(:,:,3));
LayerSpectral.L4comodcorrs_Wh.L5a = squeeze(L4comodcorrs_Wh(:,:,4));
LayerSpectral.L4comodcorrs_Wh.L56 = squeeze(L4comodcorrs_Wh(:,:,5));
LayerSpectral.L4comodcorrs_Wh.L6 = squeeze(L4comodcorrs_Wh(:,:,6));

LayerSpectral.L4comodcorrs_loP.L1 = squeeze(L4comodcorrs_loP(:,:,1));
LayerSpectral.L4comodcorrs_loP.L23 = squeeze(L4comodcorrs_loP(:,:,2));
LayerSpectral.L4comodcorrs_loP.L4 = squeeze(L4comodcorrs_loP(:,:,3));
LayerSpectral.L4comodcorrs_loP.L5a = squeeze(L4comodcorrs_loP(:,:,4));
LayerSpectral.L4comodcorrs_loP.L56 = squeeze(L4comodcorrs_loP(:,:,5));
LayerSpectral.L4comodcorrs_loP.L6 = squeeze(L4comodcorrs_loP(:,:,6));

LayerSpectral.L4comodcorrs_hiP.L1 = squeeze(L4comodcorrs_hiP(:,:,1));
LayerSpectral.L4comodcorrs_hiP.L23 = squeeze(L4comodcorrs_hiP(:,:,2));
LayerSpectral.L4comodcorrs_hiP.L4 = squeeze(L4comodcorrs_hiP(:,:,3));
LayerSpectral.L4comodcorrs_hiP.L5a = squeeze(L4comodcorrs_hiP(:,:,4));
LayerSpectral.L4comodcorrs_hiP.L56 = squeeze(L4comodcorrs_hiP(:,:,5));
LayerSpectral.L4comodcorrs_hiP.L6 = squeeze(L4comodcorrs_hiP(:,:,6));

%
LayerSpectral.L5acomodcorrs_all.L1 = squeeze(L5acomodcorrs_all(:,:,1));
LayerSpectral.L5acomodcorrs_all.L23 = squeeze(L5acomodcorrs_all(:,:,2));
LayerSpectral.L5acomodcorrs_all.L4 = squeeze(L5acomodcorrs_all(:,:,3));
LayerSpectral.L5acomodcorrs_all.L5a = squeeze(L5acomodcorrs_all(:,:,4));
LayerSpectral.L5acomodcorrs_all.L56 = squeeze(L5acomodcorrs_all(:,:,5));
LayerSpectral.L5acomodcorrs_all.L6 = squeeze(L5acomodcorrs_all(:,:,6));

LayerSpectral.L5acomodcorrs_NWh.L1 = squeeze(L5acomodcorrs_NWh(:,:,1));
LayerSpectral.L5acomodcorrs_NWh.L23 = squeeze(L5acomodcorrs_NWh(:,:,2));
LayerSpectral.L5acomodcorrs_NWh.L4 = squeeze(L5acomodcorrs_NWh(:,:,3));
LayerSpectral.L5acomodcorrs_NWh.L5a = squeeze(L5acomodcorrs_NWh(:,:,4));
LayerSpectral.L5acomodcorrs_NWh.L56 = squeeze(L5acomodcorrs_NWh(:,:,5));
LayerSpectral.L5acomodcorrs_NWh.L6 = squeeze(L5acomodcorrs_NWh(:,:,6));

LayerSpectral.L5acomodcorrs_Wh.L1 = squeeze(L5acomodcorrs_Wh(:,:,1));
LayerSpectral.L5acomodcorrs_Wh.L23 = squeeze(L5acomodcorrs_Wh(:,:,2));
LayerSpectral.L5acomodcorrs_Wh.L4 = squeeze(L5acomodcorrs_Wh(:,:,3));
LayerSpectral.L5acomodcorrs_Wh.L5a = squeeze(L5acomodcorrs_Wh(:,:,4));
LayerSpectral.L5acomodcorrs_Wh.L56 = squeeze(L5acomodcorrs_Wh(:,:,5));
LayerSpectral.L5acomodcorrs_Wh.L6 = squeeze(L5acomodcorrs_Wh(:,:,6));

LayerSpectral.L5acomodcorrs_loP.L1 = squeeze(L5acomodcorrs_loP(:,:,1));
LayerSpectral.L5acomodcorrs_loP.L23 = squeeze(L5acomodcorrs_loP(:,:,2));
LayerSpectral.L5acomodcorrs_loP.L4 = squeeze(L5acomodcorrs_loP(:,:,3));
LayerSpectral.L5acomodcorrs_loP.L5a = squeeze(L5acomodcorrs_loP(:,:,4));
LayerSpectral.L5acomodcorrs_loP.L56 = squeeze(L5acomodcorrs_loP(:,:,5));
LayerSpectral.L5acomodcorrs_loP.L6 = squeeze(L5acomodcorrs_loP(:,:,6));

LayerSpectral.L5acomodcorrs_hiP.L1 = squeeze(L5acomodcorrs_hiP(:,:,1));
LayerSpectral.L5acomodcorrs_hiP.L23 = squeeze(L5acomodcorrs_hiP(:,:,2));
LayerSpectral.L5acomodcorrs_hiP.L4 = squeeze(L5acomodcorrs_hiP(:,:,3));
LayerSpectral.L5acomodcorrs_hiP.L5a = squeeze(L5acomodcorrs_hiP(:,:,4));
LayerSpectral.L5acomodcorrs_hiP.L56 = squeeze(L5acomodcorrs_hiP(:,:,5));
LayerSpectral.L5acomodcorrs_hiP.L6 = squeeze(L5acomodcorrs_hiP(:,:,6));

%
LayerSpectral.L56comodcorrs_all.L1 = squeeze(L56comodcorrs_all(:,:,1));
LayerSpectral.L56comodcorrs_all.L23 = squeeze(L56comodcorrs_all(:,:,2));
LayerSpectral.L56comodcorrs_all.L4 = squeeze(L56comodcorrs_all(:,:,3));
LayerSpectral.L56comodcorrs_all.L5a = squeeze(L56comodcorrs_all(:,:,4));
LayerSpectral.L56comodcorrs_all.L56 = squeeze(L56comodcorrs_all(:,:,5));
LayerSpectral.L56comodcorrs_all.L6 = squeeze(L56comodcorrs_all(:,:,6));

LayerSpectral.L56comodcorrs_NWh.L1 = squeeze(L56comodcorrs_NWh(:,:,1));
LayerSpectral.L56comodcorrs_NWh.L23 = squeeze(L56comodcorrs_NWh(:,:,2));
LayerSpectral.L56comodcorrs_NWh.L4 = squeeze(L56comodcorrs_NWh(:,:,3));
LayerSpectral.L56comodcorrs_NWh.L5a = squeeze(L56comodcorrs_NWh(:,:,4));
LayerSpectral.L56comodcorrs_NWh.L56 = squeeze(L56comodcorrs_NWh(:,:,5));
LayerSpectral.L56comodcorrs_NWh.L6 = squeeze(L56comodcorrs_NWh(:,:,6));

LayerSpectral.L56comodcorrs_Wh.L1 = squeeze(L56comodcorrs_Wh(:,:,1));
LayerSpectral.L56comodcorrs_Wh.L23 = squeeze(L56comodcorrs_Wh(:,:,2));
LayerSpectral.L56comodcorrs_Wh.L4 = squeeze(L56comodcorrs_Wh(:,:,3));
LayerSpectral.L56comodcorrs_Wh.L5a = squeeze(L56comodcorrs_Wh(:,:,4));
LayerSpectral.L56comodcorrs_Wh.L56 = squeeze(L56comodcorrs_Wh(:,:,5));
LayerSpectral.L56comodcorrs_Wh.L6 = squeeze(L56comodcorrs_Wh(:,:,6));

LayerSpectral.L56comodcorrs_loP.L1 = squeeze(L56comodcorrs_loP(:,:,1));
LayerSpectral.L56comodcorrs_loP.L23 = squeeze(L56comodcorrs_loP(:,:,2));
LayerSpectral.L56comodcorrs_loP.L4 = squeeze(L56comodcorrs_loP(:,:,3));
LayerSpectral.L56comodcorrs_loP.L5a = squeeze(L56comodcorrs_loP(:,:,4));
LayerSpectral.L56comodcorrs_loP.L56 = squeeze(L56comodcorrs_loP(:,:,5));
LayerSpectral.L56comodcorrs_loP.L6 = squeeze(L56comodcorrs_loP(:,:,6));

LayerSpectral.L56comodcorrs_hiP.L1 = squeeze(L56comodcorrs_hiP(:,:,1));
LayerSpectral.L56comodcorrs_hiP.L23 = squeeze(L56comodcorrs_hiP(:,:,2));
LayerSpectral.L56comodcorrs_hiP.L4 = squeeze(L56comodcorrs_hiP(:,:,3));
LayerSpectral.L56comodcorrs_hiP.L5a = squeeze(L56comodcorrs_hiP(:,:,4));
LayerSpectral.L56comodcorrs_hiP.L56 = squeeze(L56comodcorrs_hiP(:,:,5));
LayerSpectral.L56comodcorrs_hiP.L6 = squeeze(L56comodcorrs_hiP(:,:,6));

%
LayerSpectral.L6comodcorrs_all.L1 = squeeze(L6comodcorrs_all(:,:,1));
LayerSpectral.L6comodcorrs_all.L23 = squeeze(L6comodcorrs_all(:,:,2));
LayerSpectral.L6comodcorrs_all.L4 = squeeze(L6comodcorrs_all(:,:,3));
LayerSpectral.L6comodcorrs_all.L5a = squeeze(L6comodcorrs_all(:,:,4));
LayerSpectral.L6comodcorrs_all.L56 = squeeze(L6comodcorrs_all(:,:,5));
LayerSpectral.L6comodcorrs_all.L6 = squeeze(L6comodcorrs_all(:,:,6));

LayerSpectral.L6comodcorrs_NWh.L1 = squeeze(L6comodcorrs_NWh(:,:,1));
LayerSpectral.L6comodcorrs_NWh.L23 = squeeze(L6comodcorrs_NWh(:,:,2));
LayerSpectral.L6comodcorrs_NWh.L4 = squeeze(L6comodcorrs_NWh(:,:,3));
LayerSpectral.L6comodcorrs_NWh.L5a = squeeze(L6comodcorrs_NWh(:,:,4));
LayerSpectral.L6comodcorrs_NWh.L56 = squeeze(L6comodcorrs_NWh(:,:,5));
LayerSpectral.L6comodcorrs_NWh.L6 = squeeze(L6comodcorrs_NWh(:,:,6));

LayerSpectral.L6comodcorrs_Wh.L1 = squeeze(L6comodcorrs_Wh(:,:,1));
LayerSpectral.L6comodcorrs_Wh.L23 = squeeze(L6comodcorrs_Wh(:,:,2));
LayerSpectral.L6comodcorrs_Wh.L4 = squeeze(L6comodcorrs_Wh(:,:,3));
LayerSpectral.L6comodcorrs_Wh.L5a = squeeze(L6comodcorrs_Wh(:,:,4));
LayerSpectral.L6comodcorrs_Wh.L56 = squeeze(L6comodcorrs_Wh(:,:,5));
LayerSpectral.L6comodcorrs_Wh.L6 = squeeze(L6comodcorrs_Wh(:,:,6));

LayerSpectral.L6comodcorrs_loP.L1 = squeeze(L6comodcorrs_loP(:,:,1));
LayerSpectral.L6comodcorrs_loP.L23 = squeeze(L6comodcorrs_loP(:,:,2));
LayerSpectral.L6comodcorrs_loP.L4 = squeeze(L6comodcorrs_loP(:,:,3));
LayerSpectral.L6comodcorrs_loP.L5a = squeeze(L6comodcorrs_loP(:,:,4));
LayerSpectral.L6comodcorrs_loP.L56 = squeeze(L6comodcorrs_loP(:,:,5));
LayerSpectral.L6comodcorrs_loP.L6 = squeeze(L6comodcorrs_loP(:,:,6));

LayerSpectral.L6comodcorrs_hiP.L1 = squeeze(L6comodcorrs_hiP(:,:,1));
LayerSpectral.L6comodcorrs_hiP.L23 = squeeze(L6comodcorrs_hiP(:,:,2));
LayerSpectral.L6comodcorrs_hiP.L4 = squeeze(L6comodcorrs_hiP(:,:,3));
LayerSpectral.L6comodcorrs_hiP.L5a = squeeze(L6comodcorrs_hiP(:,:,4));
LayerSpectral.L6comodcorrs_hiP.L56 = squeeze(L6comodcorrs_hiP(:,:,5));
LayerSpectral.L6comodcorrs_hiP.L6 = squeeze(L6comodcorrs_hiP(:,:,6));

% Finally saving all...
save(savefile,'-v7.3','LayerSpectral');

%% FIGURES
% cmin = min([min(min(min(L1comodcorrs_all)))...
%     ];
% cmax = 1;

figure;
for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i*1);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L1comodcorrs_all(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i+6);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L23comodcorrs_all(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i+12);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L4comodcorrs_all(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i+18);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L5acomodcorrs_all(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i+24);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L56comodcorrs_all(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i+30);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L6comodcorrs_all(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

NiceSave('LaminarCoMOD_all',figfolder,baseName)

%% FIGURE:
figure;
for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i*1);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L1comodcorrs_NWh(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i+6);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L23comodcorrs_NWh(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i+12);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L4comodcorrs_NWh(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i+18);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L5acomodcorrs_NWh(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i+24);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L56comodcorrs_NWh(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i+30);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L6comodcorrs_NWh(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

NiceSave('LaminarCoMOD_NWh',figfolder,baseName)

%% FIGURE:
figure;
for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i*1);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L1comodcorrs_Wh(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i+6);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L23comodcorrs_Wh(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i+12);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L4comodcorrs_Wh(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i+18);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L5acomodcorrs_Wh(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i+24);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L56comodcorrs_Wh(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i+30);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L6comodcorrs_Wh(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

NiceSave('LaminarCoMOD_Wh',figfolder,baseName)

%% FIGURE:
figure;
for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i*1);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L1comodcorrs_loP(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i+6);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L23comodcorrs_loP(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i+12);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L4comodcorrs_loP(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i+18);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L5acomodcorrs_loP(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i+24);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L56comodcorrs_loP(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i+30);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L6comodcorrs_loP(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

NiceSave('LaminarCoMOD_loPup',figfolder,baseName)

%% FIGURE:
figure;
for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i*1);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L1comodcorrs_hiP(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i+6);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L23comodcorrs_hiP(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i+12);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L4comodcorrs_hiP(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i+18);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L5acomodcorrs_hiP(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i+24);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L56comodcorrs_hiP(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(dLayerSpec_all,3)
    subplot(6,6,i+30);
    imagesc(log10(wavespec.freqs),log10(wavespec.freqs),L6comodcorrs_hiP(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

NiceSave('LaminarCoMOD_hiP',figfolder,baseName)

%% Laminar Phase-Amplitude coupling, by state
% PhaseAmpCouplingByAmp? Plan according to first results
