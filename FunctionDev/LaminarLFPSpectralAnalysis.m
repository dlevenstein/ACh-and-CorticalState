basePath = pwd;
baseName = bz_BasenameFromBasepath(basePath);
sessionInfo = bz_getSessionInfo(basePath,'noPrompts',true);

badchannels = sessionInfo.badchannels;
usechannels = sessionInfo.AnatGrps.Channels;
usechannels(ismember(usechannels,badchannels))=[];
channels = sessionInfo.channels;

figfolder = fullfile(basePath,'AnalysisFigures');
savefile = fullfile(basePath,[baseName,'.LaminarSpectralAnalysis.lfp.mat']);

%Pending: this would need the ultimate layer boundary detection and
%total exclusion of bad channels

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

% Touch
Piezotouch = bz_LoadBehavior(basePath,'Piezotouch');

% Pupil phase
lowfilter = [0.01 0.1];
pupil4filter = pupildilation;
lowpupildata = bz_Filter(pupil4filter,'passband',lowfilter,'filter' ,'fir1','order',3);

pupthresh = nanmedian(log10(lowpupildata.amp));
highpup = log10(lowpupildata.amp)>pupthresh;

eventshipupil = interp1(lowpupildata.timestamps(highpup),...
    lowpupildata.timestamps(highpup),t_lo,'nearest');
eventshipupil(isnan(eventshipupil)) = [];

eventslopupil = interp1(lowpupildata.timestamps(~highpup),...
    lowpupildata.timestamps(~highpup),t_lo,'nearest');
eventslopupil(isnan(eventslopupil)) = [];

%% All ints in spec times

%% Laminar Power spectra by layer, by state
%Pending... allocate space...

cLayerloSpec_all = []; cLayerhiSpec_all = []; cLayerOsci_all = []; cLayerMUA_all = [];
cLayerloSpec_NWh = []; cLayerhiSpec_NWh = [];cLayerOsci_NWh = []; cLayerMUA_NWh = [];
cLayerloSpec_Wh = []; cLayerhiSpec_Wh = []; cLayerOsci_Wh = []; cLayerMUA_Wh = [];
cLayerloSpec_T = []; cLayerhiSpec_T = []; cLayerOsci_T = []; cLayerMUA_T = [];
cLayerloSpec_loP = []; cLayerhiSpec_loP = []; cLayerOsci_loP = []; cLayerMUA_loP = [];
cLayerloSpec_hiP = []; cLayerhiSpec_hiP = []; cLayerOsci_hiP = []; cLayerMUA_hiP = [];

for i = 1:length(channels)
    i
    % Loading spectrograms
    load(fullfile(basePath,[baseName,'.',num2str(channels(i)),'.LayerID.lfp.mat']));
    
    lof = LayerID.lof;
    lospec = LayerID.lospec;
    t_lo = LayerID.t_lo;
    hif = LayerID.hif;
    hispec = LayerID.hispec;
    t_hi = LayerID.t_hi;
    MUA = LayerID.MUA;
    
    cLayerloSpec_all = cat(2,cLayerloSpec_all,nanmean(lospec,2));
    cLayerhiSpec_all = cat(2,cLayerhiSpec_all,nanmean(hispec,2));
    
    % Sorting through events
    losampRate = 1/mean(diff(t_lo));
    hisampRate = 1/mean(diff(t_hi));
    
    events = EMGwhisk.ints.Wh;
    lospec_temp = NaN(size(lospec,1),size(events,1));
    hispec_temp = NaN(size(hispec,1),size(events,1));
    
    for e = 1:size(events,1)
        %loSpec
        tempevents = round(events.*losampRate);
        if tempevents(e,1) > 0 && tempevents(e,2) < size(lospec,2)
            lospec_temp(:,e) = nanmean(lospec(:,tempevents(e,1):tempevents(e,2)),2);
        else
        end
        %hiSpec
        tempevents = round(events.*hisampRate);
        if tempevents(e,1) > 0 && tempevents(e,2) < size(hispec,2)
            hispec_temp(:,e) = nanmean(hispec(:,tempevents(e,1):tempevents(e,2)),2);
        else
        end
    end
    
    cLayerloSpec_Wh = cat(2,cLayerloSpec_Wh,nanmean(lospec_temp,2)); 
    cLayerhiSpec_Wh = cat(2,cLayerhiSpec_Wh,nanmean(hispec_temp,2));
    
    % Repeat for NWh
    events = EMGwhisk.ints.NWh;
    lospec_temp = NaN(size(lospec,1),size(events,1));
    hispec_temp = NaN(size(hispec,1),size(events,1));
    
    for e = 1:size(events,1)
        %loSpec
        tempevents = round(events.*losampRate);
        if tempevents(e,1) > 0 && tempevents(e,2) < size(lospec,2)
            lospec_temp(:,e) = nanmean(lospec(:,tempevents(e,1):tempevents(e,2)),2);
        else
        end
        %hiSpec
        tempevents = round(events.*hisampRate);
        if tempevents(e,1) > 0 && tempevents(e,2) < size(hispec,2)
            hispec_temp(:,e) = nanmean(hispec(:,tempevents(e,1):tempevents(e,2)),2);
        else
        end
    end
    
    cLayerloSpec_NWh = cat(2,cLayerloSpec_NWh,nanmean(lospec_temp,2)); 
    cLayerhiSpec_NWh = cat(2,cLayerhiSpec_NWh,nanmean(hispec_temp,2));
    
    % Repeat for Touch
    if ~isempty(Piezotouch)
        
        events = Piezotouch.ints.Touch;
        lospec_temp = NaN(size(lospec,1),size(events,1));
        hispec_temp = NaN(size(hispec,1),size(events,1));
        
        for e = 1:size(events,1)
            %loSpec
            tempevents = round(events.*losampRate);
            if tempevents(e,1) > 0 && tempevents(e,2) < size(lospec,2)
                lospec_temp(:,e) = nanmean(lospec(:,tempevents(e,1):tempevents(e,2)),2);
            else
            end
            %hiSpec
            tempevents = round(events.*hisampRate);
            if tempevents(e,1) > 0 && tempevents(e,2) < size(hispec,2)
                hispec_temp(:,e) = nanmean(hispec(:,tempevents(e,1):tempevents(e,2)),2);
            else
            end
        end
        
        cLayerloSpec_NWh = cat(2,cLayerloSpec_NWh,nanmean(lospec_temp,2));
        cLayerhiSpec_NWh = cat(2,cLayerhiSpec_NWh,nanmean(hispec_temp,2));
        
    else
    end
    
    % Repeat for >median pupil area
    
    events = lowpupildata.timestamps(highpup);
    lospec_temp = NaN(size(lospec,1),size(events,1));
    hispec_temp = NaN(size(hispec,1),size(events,1));
    
    for e = 1:size(events,1)
        %loSpec
        tempevents = round(events.*losampRate);
        if tempevents(e) > 0 && tempevents(e) < size(lospec,2)
            lospec_temp(:,e) = lospec(:,tempevents(e));
        else
        end
        %hiSpec
        tempevents = round(events.*hisampRate);
        if tempevents(e) > 0 && tempevents(e) < size(hispec,2)
            hispec_temp(:,e) = hispec(:,tempevents(e));
        else
        end
    end
    
    cLayerloSpec_hiP = cat(2,cLayerloSpec_hiP,nanmean(lospec_temp,2));
    cLayerhiSpec_hiP = cat(2,cLayerhiSpec_hiP,nanmean(hispec_temp,2));
        
    % Repeat for <median pupil area
    
    events = lowpupildata.timestamps(~highpup);
    lospec_temp = NaN(size(lospec,1),size(events,1));
    hispec_temp = NaN(size(hispec,1),size(events,1));
    
    for e = 1:size(events,1)
        %loSpec
        tempevents = round(events.*losampRate);
        if tempevents(e) > 0 && tempevents(e) < size(lospec,2)
            lospec_temp(:,e) = lospec(:,tempevents(e));
        else
        end
        %hiSpec
        tempevents = round(events.*hisampRate);
        if tempevents(e) > 0 && tempevents(e) < size(hispec,2)
            hispec_temp(:,e) = hispec(:,tempevents(e));
        else
        end
    end
    
    cLayerloSpec_loP = cat(2,cLayerloSpec_loP,nanmean(lospec_temp,2));
    cLayerhiSpec_loP = cat(2,cLayerhiSpec_loP,nanmean(hispec_temp,2));
        
end

clear LayerID lospec t_lo hispec t_hi 

% Saving to struct


%% z scores within cortex!
zsmlopspec = zeros(size(smlopspec,1),size(smlopspec,2));
for n = 1:size(smlopspec,1)
    zsmlopspec(n,:) = (smlopspec(n,:)-mean(smlopspec(n,:)))./std(smlopspec(n,:));
end

zsmhipspec = zeros(size(smhipspec,1),size(smhipspec,2));
for n = 1:size(smhipspec,1)
    zsmhipspec(n,:) = (smhipspec(n,:)-mean(smhipspec(n,:)))./std(smhipspec(n,:));
end

%% FIGURES! 
figure;
subplot(2,3,1);
imagesc(log10(lof),1:length(usechannels),cLayerloSpec_all(:,usechannels+1)')
axis xy
LogScale('x',10)
colormap('jet'); caxis([min(min(cLayerloSpec_all(:,usechannels+1)))...
    max(max(cLayerloSpec_all(:,usechannels+1)))]);
c = colorbar;
c.Label.String = 'power (dB)';
set(gca,'YDir','reverse');
set(gca,'Xtick',[log10([0.5 1 5 10 25 50 100 250])]);
set(gca,'Xticklabel',{'0.5','1','5','10','25','50','100','250'});
xlim([log10([0.5 250])]);
xtickangle(45);
xlabel('frequency (Hz)');
% set(gca,'Ytick',[0.5:1:size(smlopspec,2)-0.5]);
% set(gca,'Yticklabel',{});
% set(gca,'YGrid','on','Layer','top','GridColor',[0 0 0]);
% ylim([0 size(smlopspec,2)]);

subplot(2,3,2);
imagesc(log10(lof),1:length(usechannels),cLayerloSpec_NWh(:,usechannels+1)')
axis xy
LogScale('x',10)
colormap('jet'); caxis([min(min(cLayerloSpec_NWh(:,usechannels+1)))...
    max(max(cLayerloSpec_NWh(:,usechannels+1)))]);
c = colorbar;
c.Label.String = 'power (dB)';
set(gca,'YDir','reverse');
set(gca,'Xtick',[log10([0.5 1 5 10 25 50 100 250])]);
set(gca,'Xticklabel',{'0.5','1','5','10','25','50','100','250'});
xlim([log10([0.5 250])]);
xtickangle(45);
xlabel('frequency (Hz)');
% set(gca,'Ytick',[0.5:1:size(smlopspec,2)-0.5]);
% set(gca,'Yticklabel',{});
% set(gca,'YGrid','on','Layer','top','GridColor',[0 0 0]);
% ylim([0 size(smlopspec,2)]);

subplot(2,3,3);
imagesc(log10(lof),1:length(usechannels),cLayerloSpec_Wh(:,usechannels+1)')
axis xy
LogScale('x',10)
colormap('jet'); caxis([min(min(cLayerloSpec_Wh(:,usechannels+1)))...
    max(max(cLayerloSpec_Wh(:,usechannels+1)))]);
c = colorbar;
c.Label.String = 'power (dB)';
set(gca,'YDir','reverse');
set(gca,'Xtick',[log10([0.5 1 5 10 25 50 100 250])]);
set(gca,'Xticklabel',{'0.5','1','5','10','25','50','100','250'});
xlim([log10([0.5 250])]);
xtickangle(45);
xlabel('frequency (Hz)');
% set(gca,'Ytick',[0.5:1:size(smlopspec,2)-0.5]);
% set(gca,'Yticklabel',{});
% set(gca,'YGrid','on','Layer','top','GridColor',[0 0 0]);
% ylim([0 size(smlopspec,2)]);

subplot(2,3,4);
imagesc(log10(lof),1:length(usechannels),cLayerhiSpec_all(:,usechannels+1)')
axis xy
LogScale('x',10)
colormap('jet'); caxis([min(min(cLayerhiSpec_all(:,usechannels+1)))...
    max(max(cLayerhiSpec_all(:,usechannels+1)))]);
c = colorbar;
c.Label.String = 'power (dB)';
set(gca,'YDir','reverse');
set(gca,'Xtick',[log10([0.5 1 5 10 25 50 100 250])]);
set(gca,'Xticklabel',{'0.5','1','5','10','25','50','100','250'});
xlim([log10([0.5 250])]);
xtickangle(45);
xlabel('frequency (Hz)');
% set(gca,'Ytick',[0.5:1:size(smlopspec,2)-0.5]);
% set(gca,'Yticklabel',{});
% set(gca,'YGrid','on','Layer','top','GridColor',[0 0 0]);
% ylim([0 size(smlopspec,2)]);

subplot(2,3,5);
imagesc(log10(lof),1:length(usechannels),cLayerhiSpec_NWh(:,usechannels+1)')
axis xy
LogScale('x',10)
colormap('jet'); caxis([min(min(cLayerhiSpec_NWh(:,usechannels+1)))...
    max(max(cLayerhiSpec_NWh(:,usechannels+1)))]);
c = colorbar;
c.Label.String = 'power (dB)';
set(gca,'YDir','reverse');
set(gca,'Xtick',[log10([0.5 1 5 10 25 50 100 250])]);
set(gca,'Xticklabel',{'0.5','1','5','10','25','50','100','250'});
xlim([log10([0.5 250])]);
xtickangle(45);
xlabel('frequency (Hz)');
% set(gca,'Ytick',[0.5:1:size(smlopspec,2)-0.5]);
% set(gca,'Yticklabel',{});
% set(gca,'YGrid','on','Layer','top','GridColor',[0 0 0]);
% ylim([0 size(smlopspec,2)]);

subplot(2,3,6);
imagesc(log10(lof),1:length(usechannels),cLayerhiSpec_Wh(:,usechannels+1)')
axis xy
LogScale('x',10)
colormap('jet'); caxis([min(min(cLayerhiSpec_Wh(:,usechannels+1)))...
    max(max(cLayerhiSpec_Wh(:,usechannels+1)))]);
c = colorbar;
c.Label.String = 'power (dB)';
set(gca,'YDir','reverse');
set(gca,'Xtick',[log10([0.5 1 5 10 25 50 100 250])]);
set(gca,'Xticklabel',{'0.5','1','5','10','25','50','100','250'});
xlim([log10([0.5 250])]);
xtickangle(45);
xlabel('frequency (Hz)');
% set(gca,'Ytick',[0.5:1:size(smlopspec,2)-0.5]);
% set(gca,'Yticklabel',{});
% set(gca,'YGrid','on','Layer','top','GridColor',[0 0 0]);
% ylim([0 size(smlopspec,2)]);

%%
figure;
subplot(2,6,1); hold on;
stairs((MUAdepth-MUAdepth(1))./max(MUAdepth-MUAdepth(1)),[1:size(MUAdepth,2)],'Color','k');
set(gca,'YDir','reverse');
xlim([0 1]); ylim([0 size(MUAdepth,2)]);
xlabel('norm. MUA power (0.15-2 KHz)'); ylabel('channel');
set(gca,'Ytick',[0.5:1:size(MUAdepth,2)-0.5]);
set(gca,'Yticklabels',usechannels,'FontSize',4);
set(gca,'YGrid','on','Layer','top','GridColor',[0 0 0]);
box on;

subplot(2,6,2); hold on;
tempidx = find(LayerID.hif > 100 & LayerID.hif < 200);
muapow1 = mean(smhipspec(tempidx,:),1);
tempidx = find(LayerID.hif > 300 & LayerID.hif < 500);
muapow2 = mean(smhipspec(tempidx,:),1);
tempidx = find(LayerID.hif > 600 & LayerID.hif < 900);
muapow3 = mean(smhipspec(tempidx,:),1);
tempidx = find(LayerID.hif > 1000 & LayerID.hif < 5000);
muapow4 = mean(smhipspec(tempidx,:),1);
tempidx = find(LayerID.hif > 5000 & LayerID.hif < 10000);
muapow5 = mean(smhipspec(tempidx,:),1);
stairs((muapow1-muapow1(1))./max(muapow1-muapow1(1)),[1:size(smhipspec,2)],'Color','r');
stairs((muapow2-muapow2(1))./max(muapow2-muapow2(1)),[1:size(smhipspec,2)],'Color','g');
stairs((muapow3-muapow3(1))./max(muapow3-muapow3(1)),[1:size(smhipspec,2)],'Color','b');
stairs((muapow4-muapow4(1))./max(muapow4-muapow4(1)),[1:size(smhipspec,2)],'Color','k');
stairs((muapow5-muapow5(1))./max(muapow5-muapow5(1)),[1:size(smhipspec,2)],'Color','m');

set(gca,'YDir','reverse');
xlim([0 1]); ylim([0 size(smhipspec,2)]);
legend({'100-200 Hz','300-500 Hz','600-900 Hz','1-5 KHz','5-10 KHz'},'Location','southwest');
xlabel('norm. power');
set(gca,'Ytick',[0.5:1:size(smhipspec,2)-0.5]);
set(gca,'Yticklabel',{});
set(gca,'YGrid','on','Layer','top','GridColor',[0 0 0]);
box on;

subplot(2,6,3:4);
imagesc(log10(LayerID.lof),[0.5:size(smlopspec,2)-0.5],smlopspec')
axis xy
LogScale('x',10)
colormap('jet'); caxis([min(min(smlopspec)) max(max(smlopspec))]);
c = colorbar;
c.Label.String = 'power (dB)';
set(gca,'YDir','reverse');
set(gca,'Xtick',[log10([0.5 1 5 10 25 50 100 250])]);
set(gca,'Xticklabel',{'0.5','1','5','10','25','50','100','250'});
xlim([log10([0.5 250])]);
xtickangle(45);
xlabel('frequency (Hz)');
set(gca,'Ytick',[0.5:1:size(smlopspec,2)-0.5]);
set(gca,'Yticklabel',{});
set(gca,'YGrid','on','Layer','top','GridColor',[0 0 0]);
ylim([0 size(smlopspec,2)]);

subplot(2,6,5:6);
imagesc(log10(LayerID.hif),[0.5:size(smhipspec,2)-0.5],smhipspec')
axis xy
LogScale('x',10)
colormap('jet'); caxis([min(min(smhipspec)) max(max(smhipspec))]);
c = colorbar;
c.Label.String = 'power (dB)';
set(gca,'YDir','reverse');
set(gca,'Xtick',[log10([100 250 500 1000 2500 5000 10000])]);
set(gca,'Xticklabel',{'100','250','500','1000','2500','5000','10000'});
xlim([log10([100 10000])]);
xtickangle(45);
xlabel('frequency (Hz)');
set(gca,'Ytick',[0.5:1:size(smhipspec,2)-0.5]);
set(gca,'Yticklabel',{});
set(gca,'YGrid','on','Layer','top','GridColor',[0 0 0]);
ylim([0 size(smhipspec,2)]);

subplot(2,6,8);hold on;
tempidx = find(LayerID.hif > 100 & LayerID.hif < 200);
muapow1 = mean(zsmhipspec(tempidx,:),1);
tempidx = find(LayerID.hif > 300 & LayerID.hif < 500);
muapow2 = mean(zsmhipspec(tempidx,:),1);
tempidx = find(LayerID.hif > 600 & LayerID.hif < 900);
muapow3 = mean(zsmhipspec(tempidx,:),1);
tempidx = find(LayerID.hif > 1000 & LayerID.hif < 5000);
muapow4 = mean(zsmhipspec(tempidx,:),1);
tempidx = find(LayerID.hif > 5000 & LayerID.hif < 10000);
muapow5 = mean(zsmhipspec(tempidx,:),1);
stairs((muapow1-muapow1(1))./max(muapow1-muapow1(1)),[1:size(zsmhipspec,2)],'Color','r');
stairs((muapow2-muapow2(1))./max(muapow2-muapow2(1)),[1:size(zsmhipspec,2)],'Color','g');
stairs((muapow3-muapow3(1))./max(muapow3-muapow3(1)),[1:size(zsmhipspec,2)],'Color','b');
stairs((muapow4-muapow4(1))./max(muapow4-muapow4(1)),[1:size(zsmhipspec,2)],'Color','k');
stairs((muapow5-muapow5(1))./max(muapow5-muapow5(1)),[1:size(zsmhipspec,2)],'Color','m');

set(gca,'YDir','reverse');
xlim([0 1]); ylim([0 size(zsmhipspec,2)]);
legend({'100-200 Hz','300-500 Hz','600-900 Hz','1-5 KHz','5-10 KHz'},'Location','southwest');
xlabel('norm. power');
set(gca,'Ytick',[0.5:1:size(zsmhipspec,2)-0.5]);
set(gca,'Yticklabel',{});
set(gca,'YGrid','on','Layer','top','GridColor',[0 0 0]);
box on;

subplot(2,6,9:10);
imagesc(log10(LayerID.lof),[0.5:size(zsmlopspec,2)-0.5],zsmlopspec')
axis xy
LogScale('x',10)
colormap('jet'); caxis([max(max(zsmlopspec))*-1 max(max(zsmlopspec))]);
c = colorbar;
c.Label.String = 'power (z score)';
set(gca,'YDir','reverse');
set(gca,'Xtick',[log10([0.5 1 5 10 25 50 100 250])]);
set(gca,'Xticklabel',{'0.5','1','5','10','25','50','100','250'});
xlim([log10([0.5 250])]);
xtickangle(45);
xlabel('frequency (Hz)');
set(gca,'Ytick',[0.5:1:size(zsmlopspec,2)-0.5]);
set(gca,'Yticklabel',{});
set(gca,'YGrid','on','Layer','top','GridColor',[0 0 0]);
ylim([0 size(zsmlopspec,2)]);

subplot(2,6,11:12);
imagesc(log10(LayerID.hif),[0.5:size(zsmhipspec,2)-0.5],zsmhipspec')
axis xy
LogScale('x',10)
colormap('jet'); caxis([max(max(zsmhipspec))*-1 max(max(zsmhipspec))]);
c = colorbar;
c.Label.String = 'power (z score)';
set(gca,'YDir','reverse');
set(gca,'Xtick',[log10([100 250 500 1000 2500 5000 10000])]);
set(gca,'Xticklabel',{'100','250','500','1000','2500','5000','10000'});
xlim([log10([100 10000])]);
xtickangle(45);
xlabel('frequency (Hz)');
set(gca,'Ytick',[0.5:1:size(zsmhipspec,2)-0.5]);
set(gca,'Yticklabel',{});
set(gca,'YGrid','on','Layer','top','GridColor',[0 0 0]);
ylim([0 size(zsmhipspec,2)]);

%NiceSave('LaminarID_Spont',figfolder,baseName);

%% XCORR everything here


%% Obtaining channels for layer boundaries


%% Layer stuff...eventSpectrograms
%Struct w/ layer channels for later quantification... just the sum!

dLayerloSpec_all = []; dLayerhiSpec_all = []; dLayerOsci_all = [];
dLayerloSpec_NWh = []; dLayerhiSpec_NWh = []; dLayerOsci_NWh = [];
dLayerloSpec_Wh = []; dLayerhiSpec_Wh = []; dLayerOsci_Wh = [];
dLayerloSpec_T = []; dLayerhiSpec_T = []; dLayerOsci_T = [];
dLayerloSpec_hiP = []; dLayerhiSpec_hiP = []; dLayerOsci_hiP = [];
dLayerloSpec_loP = []; dLayerhiSpec_loP = []; dLayerOsci_loP = [];


%% Laminar Comodulogram, by state
%each layer to itself and to others...
%per states

%once average layer spectrogram.. dSpecs... .I could obtain comodulograms
%freq/freq axes w/ corr on top...

comod.corrs = corr(spec{1}',spec{2}','type','spearman');

%% Laminar LFP power distribution, by state
% also needs discrete layer info!!!
numpowerbins = 200;
%minpower = min(intspec(:)); maxpower = max(intspec(:));
switch spectype
    case 'FFT'
        minpower = 0.5;maxpower = 2.85e3;
    case 'wavelet'
        minpower = -2;maxpower = 2;
end
powerbins = linspace(log10(minpower),log10(maxpower),numpowerbins);
[powerdist_mean] = hist(log10(intspec),powerbins);

%
figure
%subplot(2,2,1)
imagesc(log2(freqlist),powerbins,powerdist_mean)
axis xy
LogScale('x',2)

%
alldists = {'birnbaumsaunders','exponential','extreme value','gamma',...
    'generalized extreme value','generalized pareto','inverse gaussian',...
    'logistic','loglogistic','lognormal','nakagami','normal','rayleigh',...
    'rician','tlocationscale','weibull'};

%Distirbutions were removed that consistently showed bad fit to cortical
%LFP data during NREM
testdists = {'gamma','loglogistic','lognormal',...
    'rayleigh','weibull'};

showexamples = [];
D = {};PF = {};bestfit={};
for ff = 1:nfreqs
    ff
    if ismember(ff,showexamples)
        [distfits] = allfitdist(intspec(:,ff),'PDF');
    else
        [distfits] = allfitdist(intspec(:,ff));
    end
    
    %Keep only the distributions tested
    distnames = {distfits(:).DistName};
    keepdists = ismember(distnames,testdists);
    distfits = distfits(keepdists);
    
    %Check the best and worst-fitting distribution
    bestfit{ff} = distfits(1).DistName;
    worstfit{ff} = distfits(end).DistName;
    
    %Sort alphabetically
    distnames = {distfits(:).DistName};
    [distnames,sortdist] = sort(distnames);
    
    D{ff} = distfits(sortdist);
end

%
AICs = cellfun(@(X) [X(:).AIC],D,'UniformOutput',false);
AICs = cat(1,AICs{:});

ndists = length(distnames);
bestfitmat = zeros(nfreqs,ndists);
for nn = 1:ndists
    bestfitmat(:,nn) = nn.*strcmp(distnames{nn},bestfit);
end
bestfitmat(bestfitmat==0)=nan;

%
figure
subplot(2,2,1)
imagesc(log2(freqlist),powerbins,powerdist_mean)
axis xy
LogScale('x',2);LogScale('y',10)
xlabel('f (Hz)')
ylabel('Power (AU)')
title([spectype,' Power Distribution'])

subplot(4,2,5)
plot(log2(freqlist),log10(AICs),'LineWidth',1)

legend(distnames,'location','southwest')
axis tight
LogScale('x',2);%LogScale('y',10)
xlabel('f (Hz)')
box off
ylabel('log(AIC)')

subplot(6,2,11)
plot(log2(freqlist),bestfitmat,'.','markersize',20)
set(gca,'YTick',1:ndists)
set(gca,'YTickLabels',distnames)
LogScale('x',2);%LogScale('y',10)
xlabel('f (Hz)')
box off

%NiceSave('LFPPowerDist',figfolder,recname)

%% Laminar Phase-Amplitude coupling, by state

% PhaseAmpCouplingByAmp seems the best bet.. think about it...

%% Laminar eventCSD/MUA
%sort thourgh high/low Wh/touch/SW amplitude events
%adjust by phase of SlowWaves

%easyyy lines but specifying better intervals, trial selection

%% MISC

% lfp = bz_GetLFP('all','basepath',basePath,'noPrompts',true);
%
% specparms.frange = [0.5 120];
% specparms.nfreqs = 100;
% specparms.space = 'log';
%
% figparms.plotname = 'LFP';
% figparms.figfolder = figfolder;
% figparms.baseName = baseName;
%
% comod = bz_Comodulogram(lfp,specparms,figparms);
% save(savefile,'comod');
