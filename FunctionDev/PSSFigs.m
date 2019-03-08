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
ColorbarWithAxis([min(min(FracEMGrho)) max(max(FracEMGrho))],['Spearman corr (Z-trasnformed)'])
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
ColorbarWithAxis([min(min(FracPupilrho)) max(max(FracPupilrho))],['Spearman corr (Z-trasnformed)'])
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
ColorbarWithAxis([-1 1],['Spearman corr (Z-trasnformed)'])
caxis([-1 1])
xlabel('lower f bound (Hz)'); ylabel('interval window (s)');
title('Frac-RSQ correlation');

NiceSave('Frac_EMG_Pupil_RSQ_Corr_p_win_lof',figfolder,baseName);

%% Laminar rescaling/resampling for plotting


% Fisher Z transformation
for f = 1:size(groupPSS.Opti.FracEMGrho,1)
    for ff = 1:size(groupPSS.Opti.FracEMGrho,2)
        groupPSS.Opti.FracEMGrho(f,ff,:) = pear_fisherz(groupPSS.Opti.FracEMGrho(f,ff,:));
        groupPSS.Opti.FracPupilrho(f,ff,:) = pear_fisherz(groupPSS.Opti.FracPupilrho(f,ff,:));
        groupPSS.Opti.Fracrsqcorr(f,ff,:) = pear_fisherz(groupPSS.Opti.Fracrsqcorr(f,ff,:));
    end
end

%% FIGURE 2: 
figure;
subplot(2,2,1);
imagesc(PSSstatsdepth.bins,1:length(usechannels),PSSstatsdepth.dist(usechannels+1,:))
ColorbarWithAxis([min(min(PSSstatsdepth.dist(usechannels+1,:))) max(max(PSSstatsdepth.dist(usechannels+1,:)))],['counts'])
caxis([min(min(PSSstatsdepth.dist(usechannels+1,:))) max(max(PSSstatsdepth.dist(usechannels+1,:)))])
xlabel('PSS distributions (au)')
xlim([-4 -1]);
ylabel('channel no. (depth-aligned)')
colormap(gca,'jet')
axis tight

subplot(4,2,2);
plot(PSSEMGlag./srate,mean(PSSEMGxcorr,2),'k','linewidth',0.5);
plot(PSSEMGlag./srate,mean(PSSEMGxcorr,2)+std(PSSEMGxcorr,0,2),'k.','linewidth',0.5);
plot(PSSEMGlag./srate,mean(PSSEMGxcorr,2)-std(PSSEMGxcorr,0,2),'k.','linewidth',0.5);
axis square
xlabel('PSS-EMG Lag (s)')
a3 = gca;
a3.XTick = [-5:2.5:5];
xlim([-5 5]);

subplot(4,2,4);
plot(PSSPuplag./srate,mean(PSSPupxcorr,2),'r','linewidth',0.5);
plot(PSSPuplag./srate,mean(PSSPupxcorr,2)+std(PSSPupxcorr,0,2),'r.','linewidth',0.5);
plot(PSSPuplag./srate,mean(PSSPupxcorr,2)-std(PSSPupxcorr,0,2),'r.','linewidth',0.5);
axis square
xlabel('PSS-Pupil Lag (s)')
a3 = gca;
a3.XTick = [-5:2.5:5];
xlim([-5 5]);

subplot(2,2,3);
plot(PSScorr.EMG(usechannels+1),1:length(usechannels),'k','linewidth',2)
hold on;
plot(PSScorr.Pup(usechannels+1),1:length(usechannels),'r','linewidth',2)
legend('PSS-EMG','PSS-Pupil','location','eastoutside')
xlabel('PSS-Behavior correlation');
set(gca,'ydir','reverse')
axis tight

subplot(2,2,4);
plot(PSSEMGtimeDiff(usechannels+1),1:length(usechannels),'k','linewidth',2)
hold on;
plot(PSSPuptimeDiff(usechannels+1),1:length(usechannels),'r','linewidth',2)
xlabel('PSS-Behavior lag (s)');
set(gca,'ydir','reverse')
axis tight
xlim([min([min(PSSEMGtimeDiff(usechannels+1)) min(PSSPuptimeDiff(usechannels+1))])-3 max([max(PSSEMGtimeDiff(usechannels+1)) max(PSSPuptimeDiff(usechannels+1))])+3]);
legend('PSS-EMG','PSS-Pupil','location','eastoutside')

NiceSave('fPSS_Behavior_CorrbyDepth',figfolder,baseName)