%%
baseGroup = ['/gpfs/data/rudylab/William/171006_WT_EM1M3';...
    '/gpfs/data/rudylab/William/171206_WT_EM1M3';...
    '/gpfs/data/rudylab/William/171209_WT_EM1M3';...
    '/gpfs/data/rudylab/William/180209_WT_EM1M3';...
    '/gpfs/data/rudylab/William/180211_WT_EM1M3';...
    '/gpfs/data/rudylab/William/180213_WT_EM1M3';...
    '/gpfs/data/rudylab/William/180214_WT_EM1M3';...
    '/gpfs/data/rudylab/William/180603_WT_EM1M3';...
    '/gpfs/data/rudylab/William/180605_WT_EM1M3';...
    '/gpfs/data/rudylab/William/180607_WT_EM1M3';...
    '/gpfs/data/rudylab/William/180703_WT_EM1M3';...
    '/gpfs/data/rudylab/William/180706_WT_EM1M3';...
    '/gpfs/data/rudylab/William/180708_WT_EM1M3';...
    '/gpfs/data/rudylab/William/180710_WT_EM1M3'];

baseGroup = ['/gpfs/data/rudylab/William/171211_KO_EM1M3';...
'/gpfs/data/rudylab/William/180208_KO_EM1M3';...
'/gpfs/data/rudylab/William/180210_KO_EM1M3';...
'/gpfs/data/rudylab/William/180212_KO_EM1M3';...
'/gpfs/data/rudylab/William/180525_KO_EM1M3';...
'/gpfs/data/rudylab/William/180530_KO_EM1M3';...
'/gpfs/data/rudylab/William/180531_KO_EM1M3';...
'/gpfs/data/rudylab/William/180601_KO_EM1M3';...
'/gpfs/data/rudylab/William/180602_KO_EM1M3';...
'/gpfs/data/rudylab/William/180608_KO_EM1M3';...
'/gpfs/data/rudylab/William/180704_KO_EM1M3';...
'/gpfs/data/rudylab/William/180705_KO_EM1M3';...
'/gpfs/data/rudylab/William/180709_KO_EM1M3';...
'/gpfs/data/rudylab/William/180711_KO_EM1M3'];

basePath = pwd;
%baseName = 'WT_EM1M3';
baseName = 'KO_EM1M3';
figfolder = fullfile(basePath,'SummaryFigures');
savefile = fullfile(basePath,[baseName,'.Optimization.PSS.lfp.mat']);

%% Rescale and collect data

groupPSS = [];
%groupDepth = [];
for i = 1:size(baseGroup,1)
    
    basename = bz_BasenameFromBasepath(baseGroup(i,:));
    
    %sessionInfo = bz_getSessionInfo(baseGroup(i,:),'noPrompts',true);
    %animaldepth = rescaleCx(baseGroup(i,:),'BADOUT',true);
    %groupDepth = bz_CollapseStruct([groupDepth animaldepth],3,'justcat',true);
    
    % For PSS-Behavior analysis
    load(fullfile(baseGroup(i,:),[basename,'.Optimization.PSS.lfp.mat']));
    PSpecSlope = PSpecSlope.Opti;
    groupPSS = bz_CollapseStruct([groupPSS PSpecSlope],3,'justcat',true);
    
end

% Saving group structs
save(savefile,'groupPSS');

%% Mean/std for group...
%Pending proper error propagation and stats...
movingwin = groupPSS.movingwin(:,:,1);
lowerbound = groupPSS.lowerbound(:,:,1);

for i = 1:size(groupPSS.FracEMGrho,3)
    [nonsigidx, nonsigidy] = find(groupPSS.FracEMGp(:,:,i) > 0.05);
    groupPSS.FracEMGrho(nonsigidx,nonsigidy,i) = NaN;
    
    [nonsigidx, nonsigidy] = find(groupPSS.FracPupilp(:,:,i) > 0.05);
    groupPSS.FracPupilrho(nonsigidx,nonsigidy,i) = NaN;
    
    [nonsigidx, nonsigidy] = find(groupPSS.Fracrsqp(:,:,i) > 0.05);
    groupPSS.Fracrsqcorr(nonsigidx,nonsigidy,i) = NaN;
end

FracEMGrho = squeeze(nanmean(groupPSS.FracEMGrho,3));
FracPupilrho = squeeze(nanmean(groupPSS.FracPupilrho,3));
Fracmeanrsq = squeeze(nanmean(groupPSS.Fracmeanrsq,3));
Fracrsqcorr = squeeze(nanmean(groupPSS.Fracrsqcorr,3));

%% FIGURES

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
ColorbarWithAxis([min(min(FracEMGrho)) max(max(FracEMGrho))],['Spearman corr'])
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

