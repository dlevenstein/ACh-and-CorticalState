basePath = pwd;
[baseFolder,baseName] = fileparts(basePath);

%% DEV
% baseName = '171209_WT_EM1M3';

%%
savefile = fullfile(basePath,[baseName,'.LaminarPowerSpec.lfp.mat']);
figfolder = fullfile(basePath,'AnalysisFigures');

%%
load(fullfile(basePath,[baseName,'.sessionInfo.mat']));
channels = sessionInfo.channels;
usechannels = sessionInfo.AnatGrps.Channels;
badchannels = sessionInfo.badchannels;
usechannels(ismember(usechannels,badchannels))=[];

%%
mergefile = fullfile(basePath,[baseName,'.MergePoints.events.mat']);
load(mergefile,'MergePoints');

spontidx = find(startsWith(MergePoints.foldernames,"Spont"));
sponttimes = [MergePoints.timestamps(spontidx(1),1) MergePoints.timestamps(spontidx(end),2)];
spontendsample = MergePoints.timestamps_samples(spontidx(end),2)+1;

%%
LOSPEC = []; HISPEC = []; MUAdepth = [];
tLOSPEC = []; tHISPEC = []; tMUAdepth = [];
for i = 1:length(usechannels)
    savfile = fullfile(basePath,[baseName,'.',num2str(usechannels(i)),'.LayerID.lfp.mat']);
    load(savfile,'LayerID');
    
    tempsidx = find(LayerID.t_lo < sponttimes(2));
    LOSPEC = cat(2,LOSPEC,mean(LayerID.lospec(:,tempsidx),2));
    tempsidx = find(LayerID.t_hi < sponttimes(2));
    HISPEC = cat(2,HISPEC,mean(LayerID.hispec(:,tempsidx),2));
    
    tempsidx = find(LayerID.MUA.timestamps < sponttimes(2));
    MUAdepth = cat(2,MUAdepth,mean(LayerID.MUA.data(tempsidx)));
    
    if spontendsample < MergePoints.timestamps_samples(end,end)+1
        tempsidx = find(LayerID.t_lo > sponttimes(2));
        tLOSPEC = cat(2,tLOSPEC,mean(LayerID.lospec(:,tempsidx),2));
        tempsidx = find(LayerID.t_hi > sponttimes(2));
        tHISPEC = cat(2,tHISPEC,mean(LayerID.hispec(:,tempsidx),2));
        
        tempsidx = find(LayerID.MUA.timestamps > sponttimes(2));
        tMUAdepth = cat(2,tMUAdepth,mean(LayerID.MUA.data(tempsidx)));
    end
end

%% Saving 
laminarpspec.LOdata = LOSPEC;
laminarpspec.HIdata = HISPEC;
laminarpspec.MUAdata = MUAdepth;
laminarpspec.tLOdata = tLOSPEC;
laminarpspec.tHIdata = tHISPEC;
laminarpspec.tMUAdata = tMUAdepth;
laminarpspec.LOfreqs = LayerID.lof;
laminarpspec.HIfreqs = LayerID.hif;

save(savefile,'laminarpspec');

%% Spatial/temporal smoothening
% smlopspec = zeros(size(LOSPEC,1),size(LOSPEC,2));
% temp_sm = 5; spat_sm = 5;
% for ch = 1:size(LOSPEC,2)
%     smlopspec(:,ch) = smooth(LOSPEC(:,ch),temp_sm,'sgolay');
% end
%
% for t = 1:size(LOSPEC,1)
%     smlopspec(t,:) = smooth(LOSPEC(t,:),spat_sm,'lowess');
% end
%
% smhipspec = zeros(size(HISPEC,1),size(HISPEC,2));
% temp_sm = 5; spat_sm = 5;
% for ch = 1:size(HISPEC,2)
%     smhipspec(:,ch) = smooth(HISPEC(:,ch),temp_sm,'sgolay');
% end
%
% for t = 1:size(HISPEC,1)
%     smhipspec(t,:) = smooth(HISPEC(t,:),spat_sm,'lowess');
% end

%% DEV
smlopspec = LOSPEC;
smhipspec = HISPEC;

%% z-scoring across freq
zsmlopspec = zeros(size(smlopspec,1),size(smlopspec,2));
for n = 1:size(smlopspec,1)
    zsmlopspec(n,:) = (smlopspec(n,:)-mean(smlopspec(n,:)))./std(smlopspec(n,:));
end

zsmhipspec = zeros(size(smhipspec,1),size(smhipspec,2));
for n = 1:size(smhipspec,1)
    zsmhipspec(n,:) = (smhipspec(n,:)-mean(smhipspec(n,:)))./std(smhipspec(n,:));
end

%% Saving SPONT figure
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

NiceSave('LaminarID_Spont',figfolder,baseName);
close all;

%%
if spontendsample < MergePoints.timestamps_samples(end,end)+1
    %% Spatial/temporal smoothening
    % smlopspec = zeros(size(tLOSPEC,1),size(tLOSPEC,2));
    % temp_sm = 5; spat_sm = 5;
    % for ch = 1:size(tLOSPEC,2)
    %     smlopspec(:,ch) = smooth(tLOSPEC(:,ch),temp_sm,'sgolay');
    % end
    %
    % for t = 1:size(tLOSPEC,1)
    %     smlopspec(t,:) = smooth(tLOSPEC(t,:),spat_sm,'lowess');
    % end
    %
    % smhipspec = zeros(size(tHISPEC,1),size(tHISPEC,2));
    % temp_sm = 5; spat_sm = 5;
    % for ch = 1:size(tHISPEC,2)
    %     smhipspec(:,ch) = smooth(tHISPEC(:,ch),temp_sm,'sgolay');
    % end
    %
    % for t = 1:size(tHISPEC,1)
    %     smhipspec(t,:) = smooth(tHISPEC(t,:),spat_sm,'lowess');
    % end
    
    %% DEV
    smlopspec = tLOSPEC;
    smhipspec = tHISPEC;
    
    %% z-scoring across freq
    zsmlopspec = zeros(size(smlopspec,1),size(smlopspec,2));
    for n = 1:size(smlopspec,1)
        zsmlopspec(n,:) = (smlopspec(n,:)-mean(smlopspec(n,:)))./std(smlopspec(n,:));
    end
    
    zsmhipspec = zeros(size(smhipspec,1),size(smhipspec,2));
    for n = 1:size(smhipspec,1)
        zsmhipspec(n,:) = (smhipspec(n,:)-mean(smhipspec(n,:)))./std(smhipspec(n,:));
    end
    
    %%
    figure;
    subplot(2,6,1); hold on;
    stairs((tMUAdepth-tMUAdepth(1))./max(tMUAdepth-tMUAdepth(1)),[1:size(tMUAdepth,2)],'Color','k');
    set(gca,'YDir','reverse');
    xlim([0 1]); ylim([0 size(tMUAdepth,2)]);
    xlabel('norm. MUA power (0.15-2 KHz)'); ylabel('channel');
    set(gca,'Ytick',[0.5:1:size(tMUAdepth,2)-0.5]);
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
    
    NiceSave('LaminarID_Touch',figfolder,baseName);
    close all;    
end