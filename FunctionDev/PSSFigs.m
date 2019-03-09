%%
basePath = pwd;
baseName = 'WT_EM1M3';
%baseName = 'KO_EM1M3';

figfolder = fullfile(basePath,'SummaryFigures');
load(fullfile(basePath,[baseName,'.PowerSpectrumSlope.lfp.mat']));

%% Mean/std/stats for group...
% Pending proper error propagation and stats...
% Fisher Z transform for corr...
% Bartlett type correction for points are not independent...

movingwin = groupPSS.Opti.movingwin(:,:,1);
lowerbound = groupPSS.Opti.lowerbound(:,:,1);

for i = 1:size(groupPSS.Opti.FracEMGrho,3)
    [nonsigidx, nonsigidy] = find(groupPSS.Opti.FracEMGp(:,:,i) > 0.05);
    groupPSS.Opti.FracEMGrho(nonsigidx,nonsigidy,i) = NaN;
    
    [nonsigidx, nonsigidy] = find(groupPSS.Opti.FracPupilp(:,:,i) > 0.05);
    groupPSS.Opti.FracPupilrho(nonsigidx,nonsigidy,i) = NaN;
    
    [nonsigidx, nonsigidy] = find(groupPSS.Opti.Fracrsqp(:,:,i) > 0.05);
    groupPSS.Opti.Fracrsqcorr(nonsigidx,nonsigidy,i) = NaN;
end

for f = 1:size(groupPSS.Opti.FracEMGrho,1)
    for ff = 1:size(groupPSS.Opti.FracEMGrho,2)
        groupPSS.Opti.FracEMGrho(f,ff,:) = pear_fisherz(groupPSS.Opti.FracEMGrho(f,ff,:));
        groupPSS.Opti.FracPupilrho(f,ff,:) = pear_fisherz(groupPSS.Opti.FracPupilrho(f,ff,:));
        groupPSS.Opti.Fracrsqcorr(f,ff,:) = pear_fisherz(groupPSS.Opti.Fracrsqcorr(f,ff,:));
    end
end

FracEMGrho = squeeze(nanmean(groupPSS.Opti.FracEMGrho,3));
FracPupilrho = squeeze(nanmean(groupPSS.Opti.FracPupilrho,3));
Fracmeanrsq = squeeze(nanmean(groupPSS.Opti.Fracmeanrsq,3));
Fracrsqcorr = squeeze(nanmean(groupPSS.Opti.Fracrsqcorr,3));

%% FIGURE 1: Exploring parameter space...
figure;
subplot(2,2,1); hold on;
winds = movingwin(:,1);
imagesc(log10(lowerbound),log10(winds),FracEMGrho')
colormap(gca,'jet')
LogScale('x',10); LogScale('y',10);
xticks(log10([1 2.5 5 10 20 40 80]));
xticklabels({'1','2.5','5','10','20','40','80'});
yticks(log10([0.25 0.5 1 2.5 5 15 30 60 90]));
yticklabels({'0.25','0.5','1','2.5','5','15','30','60','90'});
axis square
axis tight
ColorbarWithAxis([min(min(FracEMGrho)) max(max(FracEMGrho))],['Spearman corr (Z-transformed)'])
caxis([min(min(FracEMGrho)) max(max(FracEMGrho))])
xlabel('lower f bound (Hz)');ylabel('interval window (s)');
title(['Frac-EMG correlation']);

subplot(2,2,2); hold on;
imagesc(log10(lowerbound),log10(winds),FracPupilrho')
colormap(gca,'jet')
LogScale('x',10); LogScale('y',10);
xticks(log10([1 2.5 5 10 20 40 80]));
xticklabels({'1','2.5','5','10','20','40','80'});
yticks(log10([0.25 0.5 1 2.5 5 15 30 60 90]));
yticklabels({'0.25','0.5','1','2.5','5','15','30','60','90'});
axis square
axis tight
ColorbarWithAxis([min(min(FracPupilrho)) max(max(FracPupilrho))],['Spearman corr (Z-transformed)'])
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
colormap(gca,'jet')
LogScale('x',10); LogScale('y',10);
xticks(log10([1 2.5 5 10 20 40 80]));
xticklabels({'1','2.5','5','10','20','40','80'});
yticks(log10([0.25 0.5 1 2.5 5 15 30 60 90]));
yticklabels({'0.25','0.5','1','2.5','5','15','30','60','90'});
axis square
axis tight
ColorbarWithAxis([-1 1],['Spearman corr (Z-transformed)'])
caxis([-1 1])
xlabel('lower f bound (Hz)'); ylabel('interval window (s)');
title('Frac-RSQ correlation');

NiceSave('Frac_EMG_Pupil_RSQ_Corr_p_win_lof',figfolder,baseName);

%% Laminar rescaling/resampling for plotting

% PSSstatsdepth = groupPSS.Shortwin.PSSstatsdepth;
% PSSstatsdepth.bins = squeeze(PSSstatsdepth.bins(:,:,1));
% PSScorr =  groupPSS.Shortwin.PSScorr;
% PSSEMGtimeDiff = groupPSS.Shortwin.PSSEMGtimeDiff;
% PSSPuptimeDiff = groupPSS.Shortwin.PSSPuptimeDiff;

PSSstatsdepth = groupPSS.Longwin.PSSstatsdepth;
PSSstatsdepth.bins = squeeze(PSSstatsdepth.bins(:,:,1));
PSScorr =  groupPSS.Longwin.PSScorr;
PSSEMGtimeDiff = groupPSS.Longwin.PSSEMGtimeDiff;
PSSPuptimeDiff = groupPSS.Longwin.PSSPuptimeDiff;

%
numberchans = 30;
normcolumn = [0:1/numberchans:1-1/numberchans];

%
normPSSdist = NaN(numberchans,length(PSSstatsdepth.bins),size(PSSstatsdepth.dist,3));
normPSSEMG = NaN(numberchans,size(PSSstatsdepth.dist,3));
normPSSEMGlag = NaN(numberchans,size(PSSstatsdepth.dist,3));
normPSSPup = NaN(numberchans,size(PSSstatsdepth.dist,3));
normPSSPuplag = NaN(numberchans,size(PSSstatsdepth.dist,3));

for i = 1:size(PSSstatsdepth.dist,3)
    [bincounts,ind]=histc(squeeze(groupDepth.ndepth(:,:,i)),normcolumn);
    for ii = 1:max(ind)
        tempidx = find(ind == ii);
        
        normPSSdist(ii,:,i) = squeeze(nanmean(PSSstatsdepth.dist(tempidx,:,i),1));
        normPSSEMG(ii,i) = squeeze(nanmean(PSScorr.EMG(tempidx,:,i),1));
        normPSSEMGlag(ii,i) = squeeze(nanmean(PSSEMGtimeDiff(tempidx,:,i),1));
        normPSSPup(ii,i) = squeeze(nanmean(PSScorr.Pup(tempidx,:,i),1));
        normPSSPuplag(ii,i) = squeeze(nanmean(PSSPuptimeDiff(tempidx,:,i),1));
    end
end
% handle p values...

PSSstatsdepth.dist = squeeze(nanmean(normPSSdist,3));
PSSEMGtimeDiff = nanmean(normPSSEMGlag,2);
PSSPuptimeDiff = nanmean(normPSSPuplag,2);

% Fisher Z transformation
for f = 1:size(normPSSEMG,1)
    normPSSEMG(f,:) = pear_fisherz(normPSSEMG(f,:));
    normPSSPup(f,:) = pear_fisherz(normPSSPup(f,:));
end
PSScorr.EMG =  nanmean(normPSSEMG,2);
PSScorr.Pup =  nanmean(normPSSPup,2);

% Spatial smoothening
spat_sm = 5;
for t = 1:size(PSSstatsdepth.dist,2)
    PSSstatsdepth.dist(:,t) = smooth(PSSstatsdepth.dist(:,t),spat_sm,'lowess');
end

spat_sm = 5;
for t = 1:size(PSScorr.EMG,2)
    PSScorr.EMG(:,t) = smooth(PSScorr.EMG(:,t),spat_sm,'lowess');
    PSScorr.Pup(:,t) = smooth(PSScorr.Pup(:,t),spat_sm,'lowess');
    PSSEMGtimeDiff(:,t) = smooth(PSSEMGtimeDiff(:,t),spat_sm,'lowess');
    PSSPuptimeDiff(:,t) = smooth(PSSPuptimeDiff(:,t),spat_sm,'lowess');
end

% Normalize
PSSstatsdepth.dist = PSSstatsdepth.dist-min(min(PSSstatsdepth.dist));
PSSstatsdepth.dist = PSSstatsdepth.dist./max(max(PSSstatsdepth.dist));

%% FIGURE 2:
figure;
subplot(2,2,1);
imagesc(PSSstatsdepth.bins,normcolumn,PSSstatsdepth.dist)
ColorbarWithAxis([min(min(PSSstatsdepth.dist)) max(max(PSSstatsdepth.dist))],['normalized counts'])
caxis([min(min(PSSstatsdepth.dist)) max(max(PSSstatsdepth.dist))])
xlabel('PSS (au)')
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
colormap(gca,'jet')
axis tight;
xlim([-4 0]); %ylim([0 1]);

subplot(2,2,3);
plot(PSScorr.EMG,normcolumn,'k','linewidth',2); hold on;
h1 = plot(PSScorr.EMG-(nanstd(normPSSEMG,1,2)./sqrt(size(normPSSEMG,2))),...
    normcolumn,'k.','linewidth',2)
h2 = plot(PSScorr.EMG+(nanstd(normPSSEMG,1,2)./sqrt(size(normPSSEMG,2))),...
    normcolumn,'k.','linewidth',2)
plot(PSScorr.Pup,normcolumn,'r','linewidth',2)
h3 = plot(PSScorr.Pup-(nanstd(normPSSPup,1,2)./sqrt(size(normPSSPup,2))),...
    normcolumn,'r.','linewidth',2)
h4 = plot(PSScorr.Pup+(nanstd(normPSSPup,1,2)./sqrt(size(normPSSPup,2))),...
    normcolumn,'r.','linewidth',2)
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','k','GridAlpha',0.45);
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend('PSS-EMG','PSS-Pupil','location','eastoutside')
xlabel('Spearman corr (Z-transformed)');
set(gca,'ydir','reverse')
axis tight

subplot(2,2,4);
plot(PSSEMGtimeDiff,normcolumn,'k','linewidth',2);hold on;
h1 = plot(PSSEMGtimeDiff-(nanstd(normPSSEMGlag,1,2)./sqrt(size(normPSSEMGlag,2))),...
    normcolumn,'k.','linewidth',2)
h2 = plot(PSSEMGtimeDiff+(nanstd(normPSSEMGlag,1,2)./sqrt(size(normPSSEMGlag,2))),...
    normcolumn,'k.','linewidth',2)
plot(PSSPuptimeDiff,normcolumn,'r','linewidth',2)
h3 = plot(PSSPuptimeDiff-(nanstd(normPSSPuplag,1,2)./sqrt(size(normPSSPuplag,2))),...
    normcolumn,'r.','linewidth',2)
h4 = plot(PSSPuptimeDiff+(nanstd(normPSSPuplag,1,2)./sqrt(size(normPSSPuplag,2))),...
    normcolumn,'r.','linewidth',2)
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','k','GridAlpha',0.45);
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend('PSS-EMG','PSS-Pupil','location','eastoutside')
xlabel('Xcorr lag (s)');
set(gca,'ydir','reverse')
axis tight

%NiceSave('shortwinPSS_Behavior_CorrbyDepth',figfolder,baseName)
NiceSave('longwinPSS_Behavior_CorrbyDepth',figfolder,baseName)

%%
% Oscicorr = groupPSS.Shortwin.Oscicorr;
% dcorr = groupPSS.Shortwin.dcorr;
% tcorr = groupPSS.Shortwin.tcorr;
% gcorr = groupPSS.Shortwin.gcorr;

Oscicorr = groupPSS.Longwin.Oscicorr;
dcorr = groupPSS.Longwin.dcorr;
tcorr = groupPSS.Longwin.tcorr;
gcorr = groupPSS.Longwin.gcorr;

normOsciEMG = NaN(size(Oscicorr.EMG,1),numberchans,size(Oscicorr.EMG,3));
normOsciPup = NaN(size(Oscicorr.Pup,1),numberchans,size(Oscicorr.Pup,3));
normdEMG = NaN(numberchans,size(dcorr.EMG,3));
normdPup = NaN(numberchans,size(dcorr.Pup,3));
normtEMG = NaN(numberchans,size(tcorr.EMG,3));
normtPup = NaN(numberchans,size(tcorr.Pup,3));
normgEMG = NaN(numberchans,size(gcorr.EMG,3));
normgPup = NaN(numberchans,size(gcorr.Pup,3));

for i = 1:size(Oscicorr.EMG,3)
    [bincounts,ind]=histc(squeeze(groupDepth.ndepth(:,:,i)),normcolumn);
    for ii = 1:max(ind)
        tempidx = find(ind == ii);
        
        normOsciEMG(:,ii,i) = squeeze(nanmean(Oscicorr.EMG(:,tempidx,i),2));
        normOsciPup(:,ii,i) = squeeze(nanmean(Oscicorr.Pup(:,tempidx,i),2));
        
        normdEMG(ii,i) = squeeze(nanmean(dcorr.EMG(tempidx,i),1));
        normdPup(ii,i) = squeeze(nanmean(dcorr.Pup(tempidx,i),1));
        normtEMG(ii,i) = squeeze(nanmean(tcorr.EMG(tempidx,i),1));
        normtPup(ii,i) = squeeze(nanmean(tcorr.Pup(tempidx,i),1));
        normgEMG(ii,i) = squeeze(nanmean(gcorr.EMG(tempidx,i),1));
        normgPup(ii,i) = squeeze(nanmean(gcorr.Pup(tempidx,i),1));
    end
end
% handle p values...

% Fisher Z transformation
for f = 1:size(normOsciEMG,1)
    for ff = 1:size(normOsciEMG,2)
        normOsciEMG(f,ff,:) = pear_fisherz(normOsciEMG(f,ff,:));
        normOsciPup(f,ff,:) = pear_fisherz(normOsciPup(f,ff,:));
    end
end
Oscicorr.EMG = squeeze(nanmean(normOsciEMG,3));
Oscicorr.Pup = squeeze(nanmean(normOsciPup,3));

for f = 1:size(normdEMG,1)
    normdEMG(f,:) = pear_fisherz(normdEMG(f,:));
    normdPup(f,:) = pear_fisherz(normdPup(f,:));
    normtEMG(f,:) = pear_fisherz(normtEMG(f,:));
    normtPup(f,:) = pear_fisherz(normtPup(f,:));
    normgEMG(f,:) = pear_fisherz(normgEMG(f,:));
    normgPup(f,:) = pear_fisherz(normgPup(f,:));
end
dcorr.EMG =  nanmean(normdEMG,2);
dcorr.Pup =  nanmean(normdPup,2);
tcorr.EMG =  nanmean(normtEMG,2);
tcorr.Pup =  nanmean(normtPup,2);
gcorr.EMG =  nanmean(normgEMG,2);
gcorr.Pup =  nanmean(normgPup,2);

% Spatial smoothening
spat_sm = 5;
for t = 1:size(Oscicorr.EMG,1)
    Oscicorr.EMG(t,:) = smooth(Oscicorr.EMG(t,:),spat_sm,'lowess');
    Oscicorr.Pup(t,:) = smooth(Oscicorr.Pup(t,:),spat_sm,'lowess');
end

spat_sm = 5;
for t = 1:size(dcorr.EMG,2)
    dcorr.EMG(:,t) = smooth(dcorr.EMG(:,t),spat_sm,'lowess');
    dcorr.Pup(:,t) = smooth(dcorr.Pup(:,t),spat_sm,'lowess');
    tcorr.EMG(:,t) = smooth(tcorr.EMG(:,t),spat_sm,'lowess');
    tcorr.Pup(:,t) = smooth(tcorr.Pup(:,t),spat_sm,'lowess');
    gcorr.EMG(:,t) = smooth(gcorr.EMG(:,t),spat_sm,'lowess');
    gcorr.Pup(:,t) = smooth(gcorr.Pup(:,t),spat_sm,'lowess');
end

%% FIGURE 3:
%Eventually fuckin' save the freqsssss when running PSS analysis!!!
%freq = fastfreq;
freq = slowfreq;

mincaxis = min([min(min(Oscicorr.EMG)) min(min(Oscicorr.Pup))]);
maxcaxis = max([max(max(Oscicorr.EMG)) max(max(Oscicorr.Pup))]);

figure;
subplot(2,2,1);
imagesc(log10(freq),normcolumn,Oscicorr.EMG')
LogScale('x',10);
colormap(gca,'jet')
ColorbarWithAxis([mincaxis maxcaxis],['Spearman corr (Z-transformed)'])
caxis([mincaxis maxcaxis])
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.65);
xlabel('f (Hz)');
title('Oscillatory LFP-EMG correlation');
axis square

subplot(2,2,3);
imagesc(log10(freq),normcolumn,Oscicorr.Pup')
LogScale('x',10);
colormap(gca,'jet')
ColorbarWithAxis([mincaxis maxcaxis],['Spearman corr (Z-transformed)'])
caxis([mincaxis maxcaxis])
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.65);
xlabel('f (Hz)'); ylabel('channel no. (depth-aligned)');
title('Oscillatory LFP-Pupil diameter correlation');
axis square

subplot(2,2,2);
plot(dcorr.EMG,normcolumn,'r','linewidth',2); hold on;
h1 = plot(dcorr.EMG-(nanstd(normdEMG,1,2)./sqrt(size(normdEMG,2))),...
    normcolumn,'r.','linewidth',2)
h2 = plot(dcorr.EMG+(nanstd(normdEMG,1,2)./sqrt(size(normdEMG,2))),...
    normcolumn,'r.','linewidth',2)
plot(tcorr.EMG,normcolumn,'g','linewidth',2); hold on;
h3 = plot(tcorr.EMG-(nanstd(normtEMG,1,2)./sqrt(size(normtEMG,2))),...
    normcolumn,'g.','linewidth',2)
h4 = plot(tcorr.EMG+(nanstd(normtEMG,1,2)./sqrt(size(normtEMG,2))),...
    normcolumn,'g.','linewidth',2)
plot(gcorr.EMG,normcolumn,'b','linewidth',2); hold on;
h5 = plot(gcorr.EMG-(nanstd(normgEMG,1,2)./sqrt(size(normgEMG,2))),...
    normcolumn,'b.','linewidth',2)
h6 = plot(gcorr.EMG+(nanstd(normgEMG,1,2)./sqrt(size(normgEMG,2))),...
    normcolumn,'b.','linewidth',2)
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','k','GridAlpha',0.45);
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h5,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h6,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend('d','t','g','location','southeast')
xlabel('Oscillatory-EMG Spearman corr (Z-transformed)');
set(gca,'ydir','reverse')
axis tight
axis square

subplot(2,2,4);
plot(dcorr.Pup,normcolumn,'r','linewidth',2); hold on;
h1 = plot(dcorr.Pup-(nanstd(normdPup,1,2)./sqrt(size(normdPup,2))),...
    normcolumn,'r.','linewidth',2)
h2 = plot(dcorr.Pup+(nanstd(normdPup,1,2)./sqrt(size(normdPup,2))),...
    normcolumn,'r.','linewidth',2)
plot(tcorr.Pup,normcolumn,'g','linewidth',2); hold on;
h3 = plot(tcorr.Pup-(nanstd(normtPup,1,2)./sqrt(size(normtPup,2))),...
    normcolumn,'g.','linewidth',2)
h4 = plot(tcorr.Pup+(nanstd(normtPup,1,2)./sqrt(size(normtPup,2))),...
    normcolumn,'g.','linewidth',2)
plot(gcorr.Pup,normcolumn,'b','linewidth',2); hold on;
h5 = plot(gcorr.Pup-(nanstd(normgPup,1,2)./sqrt(size(normgPup,2))),...
    normcolumn,'b.','linewidth',2)
h6 = plot(gcorr.Pup+(nanstd(normgPup,1,2)./sqrt(size(normgPup,2))),...
    normcolumn,'b.','linewidth',2)
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','k','GridAlpha',0.45);
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h5,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h6,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend('d','t','g','location','southeast')
xlabel('Oscillatory-Pupil Spearman corr (Z-transformed)');
set(gca,'ydir','reverse')
axis tight
axis square

%NiceSave('shortwinOsci_Behavior_CorrbyDepth',figfolder,baseName)
NiceSave('longwinOsci_Behavior_CorrbyDepth',figfolder,baseName)

%%
% PSSxcorr = groupPSS.Shortwin.PSSxcorr;
% Oscixcorr = groupPSS.Shortwin.Oscixcorr;
% dxcorr = groupPSS.Shortwin.dxcorr;
% txcorr = groupPSS.Shortwin.txcorr;
% gxcorr = groupPSS.Shortwin.gxcorr;

PSSxcorr = groupPSS.Longwin.PSSxcorr;
Oscixcorr = groupPSS.Longwin.Oscixcorr;
dxcorr = groupPSS.Longwin.dxcorr;
txcorr = groupPSS.Longwin.txcorr;
gxcorr = groupPSS.Longwin.gxcorr;

%
numberchans = 30;
normcolumn = [0:1/numberchans:1-1/numberchans];

%
tempxPSS = NaN(numberchans,size(PSSxcorr,2),size(PSSxcorr,3));
tempxOsci = NaN(numberchans,size(PSSxcorr,2),size(PSSxcorr,3));
tempxd = NaN(numberchans,size(PSSxcorr,2),size(PSSxcorr,3));
tempxt = NaN(numberchans,size(PSSxcorr,2),size(PSSxcorr,3));
tempxg = NaN(numberchans,size(PSSxcorr,2),size(PSSxcorr,3));

for i = 1:size(PSSxcorr,3)
    [bincounts,ind]=histc(squeeze(groupDepth.ndepth(:,:,i)),normcolumn);
    for ii = 1:max(ind)
        tempidx = find(ind == ii);
        
        tempxPSS(ii,:,i) = squeeze(nanmean(PSSxcorr(tempidx,:,i),1));
        tempxOsci(ii,:,i) = squeeze(nanmean(Oscixcorr(tempidx,:,i),1));
        tempxd(ii,:,i) = squeeze(nanmean(dxcorr(tempidx,:,i),1));
        tempxt(ii,:,i) = squeeze(nanmean(txcorr(tempidx,:,i),1));
        tempxg(ii,:,i) = squeeze(nanmean(gxcorr(tempidx,:,i),1));
    end
end
% handle p values...

%
normxPSS = NaN(numberchans,numberchans,size(PSSxcorr,3));
normxOsci = NaN(numberchans,numberchans,size(PSSxcorr,3));
normxd = NaN(numberchans,numberchans,size(PSSxcorr,3));
normxt = NaN(numberchans,numberchans,size(PSSxcorr,3));
normxg = NaN(numberchans,numberchans,size(PSSxcorr,3));

for i = 1:size(PSSxcorr,3)
    [bincounts,ind]=histc(squeeze(groupDepth.ndepth(:,:,i)),normcolumn);
    for ii = 1:max(ind)
        tempidx = find(ind == ii);
        
        normxPSS(:,ii,i) = squeeze(nanmean(tempxPSS(:,tempidx,i),2));
        normxOsci(:,ii,i) = squeeze(nanmean(tempxOsci(:,tempidx,i),2));
        normxd(:,ii,i) = squeeze(nanmean(tempxd(:,tempidx,i),2));
        normxt(:,ii,i) = squeeze(nanmean(tempxt(:,tempidx,i),2));
        normxg(:,ii,i) = squeeze(nanmean(tempxg(:,tempidx,i),2));
    end
end
% handle p values...

% Fisher Z transformation
for f = 1:size(normxPSS,1)
    for ff = 1:size(normxPSS,2)
        normxPSS(f,ff,:) = pear_fisherz(normxPSS(f,ff,:));
        normxOsci(f,ff,:) = pear_fisherz(normxOsci(f,ff,:));
        normxd(f,ff,:) = pear_fisherz(normxd(f,ff,:));
        normxt(f,ff,:) = pear_fisherz(normxt(f,ff,:));
        normxg(f,ff,:) = pear_fisherz(normxg(f,ff,:));
    end
end
PSSxcorr = squeeze(nanmean(normxPSS,3));
Oscixcorr = squeeze(nanmean(normxOsci,3));
dxcorr = squeeze(nanmean(normxd,3));
txcorr = squeeze(nanmean(normxt,3));
gxcorr = squeeze(nanmean(normxg,3));

%% FIGURE 4:
figure;
subplot(2,3,1);
imagesc(normcolumn,normcolumn,PSSxcorr);
colormap(gca,'jet')
axis square
%set(gca,'ydir','reverse')
ColorbarWithAxis([min(min(PSSxcorr)) 2],['Spearman corr (Z-transformed)'])
caxis([min(min(PSSxcorr)) 2])
% ColorbarWithAxis([min(min(PSSxcorr)) max(max(PSSxcorr))],['Spearman corr (Z-transformed)'])
% caxis([min(min(PSSxcorr)) max(max(PSSxcorr))])
set(gca,'Xtick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Xticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'XGrid','on', 'GridColor','w','GridAlpha',0.65);
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.65);
title('PSS-PSS xcorr by depth');

subplot(2,3,2);
imagesc(normcolumn,normcolumn,Oscixcorr);
colormap(gca,'jet')
axis square
%set(gca,'ydir','reverse')
ColorbarWithAxis([min(min(Oscixcorr)) 2.5],['Spearman corr (Z-transformed)'])
caxis([min(min(Oscixcorr)) 2.5])
% ColorbarWithAxis([min(min(PSSxcorr)) max(max(PSSxcorr))],['Spearman corr (Z-transformed)'])
% caxis([min(min(PSSxcorr)) max(max(PSSxcorr))])
set(gca,'Xtick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Xticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'XGrid','on', 'GridColor','w','GridAlpha',0.65);
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.65);
title('Oscillatory-Oscillatory LFP xcorr by depth');

subplot(2,3,4);
imagesc(normcolumn,normcolumn,dxcorr);
colormap(gca,'jet')
axis square
%set(gca,'ydir','reverse')
ColorbarWithAxis([min(min(dxcorr)) 1.5],['Spearman corr (Z-transformed)'])
caxis([min(min(dxcorr)) 1.5])
% ColorbarWithAxis([min(min(PSSxcorr)) max(max(PSSxcorr))],['Spearman corr (Z-transformed)'])
% caxis([min(min(PSSxcorr)) max(max(PSSxcorr))])
set(gca,'Xtick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Xticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'XGrid','on', 'GridColor','w','GridAlpha',0.65);
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.65);
title('Delta-Delta LFP xcorr by depth');

subplot(2,3,5);
imagesc(normcolumn,normcolumn,txcorr);
colormap(gca,'jet')
axis square
%set(gca,'ydir','reverse')
ColorbarWithAxis([min(min(txcorr)) 1.5],['Spearman corr (Z-transformed)'])
caxis([min(min(txcorr)) 1.5])
% ColorbarWithAxis([min(min(PSSxcorr)) max(max(PSSxcorr))],['Spearman corr (Z-transformed)'])
% caxis([min(min(PSSxcorr)) max(max(PSSxcorr))])
set(gca,'Xtick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Xticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'XGrid','on', 'GridColor','w','GridAlpha',0.65);
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.65);
title('Theta-Theta LFP xcorr by depth');

subplot(2,3,6);
imagesc(normcolumn,normcolumn,gxcorr);
colormap(gca,'jet')
axis square
%set(gca,'ydir','reverse')
ColorbarWithAxis([min(min(gxcorr)) 1.5],['Spearman corr (Z-transformed)'])
caxis([min(min(gxcorr)) 1.5])
% ColorbarWithAxis([min(min(PSSxcorr)) max(max(PSSxcorr))],['Spearman corr (Z-transformed)'])
% caxis([min(min(PSSxcorr)) max(max(PSSxcorr))])
set(gca,'Xtick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Xticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'XGrid','on', 'GridColor','w','GridAlpha',0.65);
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.65);
title('Gamma-Gamma LFP xcorr by depth');

%NiceSave('shortwinPSS_Osci_Xcorrbydepth',figfolder,baseName);
NiceSave('longwinPSS_Osci_Xcorrbydepth',figfolder,baseName);