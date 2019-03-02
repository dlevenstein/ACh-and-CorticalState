basePath = pwd;
baseName = bz_BasenameFromBasepath(basePath);
figfolder = fullfile(basePath,'PSSstructs');

%%
lowerbound = [1:20];
upperbound = [60:120];

%%
grandWhiskPSScorr.EMG = zeros(length(lowerbound),length(upperbound));
grandWhiskPSScorr.pup = zeros(length(lowerbound),length(upperbound));
grandWhiskPSScorr.dpdt = zeros(length(lowerbound),length(upperbound));
grandWhiskPSScorr.phasecoupling = zeros(length(lowerbound),length(upperbound));

for f = 1:length(lowerbound)
    for ff = 1:length(upperbound)
        load(fullfile(figfolder,[baseName,'.',num2str(f),'.',num2str(ff),'.WhiskPSScorr.lfp.mat']),'WhiskPSScorr');
        
        grandWhiskPSScorr.EMG(f,ff) = max(WhiskPSScorr.EMG);
        grandWhiskPSScorr.pup(f,ff) = max(WhiskPSScorr.pup);
        grandWhiskPSScorr.dpdt(f,ff) = max(WhiskPSScorr.dpdt);
        grandWhiskPSScorr.phasecoupling(f,ff) = max(WhiskPSScorr.phasecoupling);
        
    end
end

save(fullfile(basePath,[baseName,'win2.dt0.5.grandWhiskPSScorr.lfp.mat']),'grandWhiskPSScorr');

%% FIGURE
rwbcolormap = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);

figure;

subplot(2,2,1);
imagesc(lowerbound,upperbound,grandWhiskPSScorr.EMG)
colormap(gca,rwbcolormap)
axis xy
axis tight
ColorbarWithAxis([min(min(grandWhiskPSScorr.EMG)) max(max(grandWhiskPSScorr.EMG))],['Spearman corr'])
caxis([min(min(grandWhiskPSScorr.EMG)) max(max(grandWhiskPSScorr.EMG))])
xlabel('lower f bound (Hz)');ylabel('upper f bound (Hz)');
title('PSS-EMG');

subplot(2,2,2);
imagesc(lowerbound,upperbound,grandWhiskPSScorr.pup)
colormap(gca,rwbcolormap)
axis xy
axis tight
ColorbarWithAxis([min(min(grandWhiskPSScorr.pup)) max(max(grandWhiskPSScorr.pup))],['Spearman corr'])
caxis([min(min(grandWhiskPSScorr.pup)) max(max(grandWhiskPSScorr.pup))])
xlabel('lower f bound (Hz)');ylabel('upper f bound (Hz)');
title('PSS-Pupil diameter');

subplot(2,2,3);
imagesc(lowerbound,upperbound,grandWhiskPSScorr.dpdt)
colormap(gca,rwbcolormap)
axis xy
axis tight
ColorbarWithAxis([min(min(grandWhiskPSScorr.dpdt)) max(max(grandWhiskPSScorr.dpdt))],['Spearman corr'])
caxis([min(min(grandWhiskPSScorr.dpdt)) max(max(grandWhiskPSScorr.dpdt))])
xlabel('lower f bound (Hz)');ylabel('upper f bound (Hz)');
title('PSS-dPdt');

subplot(2,2,4);
imagesc(lowerbound,upperbound,grandWhiskPSScorr.phasecoupling)
colormap(gca,rwbcolormap)
axis xy
axis tight
ColorbarWithAxis([min(min(grandWhiskPSScorr.phasecoupling)) max(max(grandWhiskPSScorr.phasecoupling))],['Spearman corr'])
caxis([min(min(grandWhiskPSScorr.phasecoupling)) max(max(grandWhiskPSScorr.phasecoupling))])
xlabel('lower f bound (Hz)');ylabel('upper f bound (Hz)');
title('PSS-Pupil phase coupling');

NiceSave('PSScorrs',figfolder,baseName);