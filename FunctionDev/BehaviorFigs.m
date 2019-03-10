%%
basePath = pwd;
baseName = 'EM1M3';

figfolder = fullfile(basePath,'SummaryFigures');

load(fullfile(basePath,['WT_EM1M3.BehaviorAnalysis.mat']));
WTBehavior = groupBehavior;
load(fullfile(basePath,['KO_EM1M3.BehaviorAnalysis.mat']));
KOBehavior = groupBehavior;

%%

WTpuphist = WTBehavior.puphist;
WTpupdthist = WTBehavior.pupdthist;
for i = 1:size(WTpuphist.counts,3)
    WTpuphist.counts(:,:,i) = WTpuphist.counts(:,:,i)./max(WTpuphist.counts(:,:,i));
    WTpupdthist.counts(:,:,i) = WTpupdthist.counts(:,:,i)./max(max(WTpupdthist.counts(:,:,i)));
end
WTpupPSD = WTBehavior.pupPSD;

KOpuphist = KOBehavior.puphist;
KOpupdthist = KOBehavior.pupdthist;
for i = 1:size(KOpuphist.counts,3)
    KOpuphist.counts(:,:,i) = KOpuphist.counts(:,:,i)./max(KOpuphist.counts(:,:,i));
    KOpupdthist.counts(:,:,i) = KOpupdthist.counts(:,:,i)./max(max(KOpupdthist.counts(:,:,i)));
end
KOpupPSD = KOBehavior.pupPSD;

%% FIGURE 1:
figure;
subplot(2,2,1);
bar(WTpuphist.bins(:,:,1),...
    nanmean(WTpuphist.counts,3),...
    'facecolor','k','facealpha',0.85); hold on;
errorbar(WTpuphist.bins(:,:,1),...
    nanmean(WTpuphist.counts,3),...
    nanstd(WTpuphist.counts,0,3),'k.');
bar(KOpuphist.bins(:,:,1),...
    nanmean(KOpuphist.counts,3),...
    'facecolor','r','facealpha',0.65)
errorbar(KOpuphist.bins(:,:,1),...
    nanmean(KOpuphist.counts,3),...
    nanstd(KOpuphist.counts,0,3),'r.');
axis tight;
xlim([0 3]); ylim([0 1]);
xlabel('diameter (norm.)'); ylabel('counts (norm.)');

subplot(2,2,2);
shadederror(log10(WTpupPSD.freqs(:,:,1)),...
    nanmean(WTpupPSD.psd,3)',...
    nanstd(WTpupPSD.psd,0,3)','k','none',1,'k');
shadederror(log10(KOpupPSD.freqs(:,:,1)),...
    nanmean(KOpupPSD.psd,3)',...
    nanstd(KOpupPSD.psd,0,3)','r','none',1,'r');
LogScale('x',10)
xlabel('f (Hz)'); ylabel('power (dB)')
axis tight
ylim([-1.5 1.5]);
title('Pupil diameter power spectrum');

subplot(2,2,3);
imagesc(WTpupdthist.bins{:,1,1},WTpupdthist.bins{:,2,1},...
    nanmean(WTpupdthist.counts,3)'); hold on;
colormap(gca,[1 1 1; colormap('jet')])
plot(get(gca,'xlim'),[0 0],'r-')
ColorbarWithAxis([-0.1 1],['counts (au)'])
caxis([-0.1 1]);
xlabel('diameter (norm)'); ylabel('dP/dt')
LogScale('x',10);
axis square
ylim([-0.25 0.25]); xlim([-0.25 0.25]);
title('floxed M1M3');

subplot(2,2,4);
imagesc(KOpupdthist.bins{:,1,1},KOpupdthist.bins{:,2,1},...
    nanmean(KOpupdthist.counts,3)'); hold on;
colormap(gca,[1 1 1; colormap('jet')])
plot(get(gca,'xlim'),[0 0],'r-')
ColorbarWithAxis([-0.1 1],['counts (au)'])
caxis([-0.1 1]);
xlabel('diameter (norm)'); ylabel('dP/dt')
LogScale('x',10);
axis square
ylim([-0.25 0.25]); xlim([-0.25 0.25]);
title('Emx1;M1M3');

NiceSave('PupilStats',figfolder,baseName)

%%
WTEMGhist = WTBehavior.EMGhist;
WTWhdurhist = WTBehavior.Whdurhist;
for i = 1:size(WTWhdurhist.Whdurs,3)
    WTEMGhist.logcounts(:,:,i) = WTEMGhist.logcounts(:,:,i)./max(WTEMGhist.logcounts(:,:,i));
    WTWhdurhist.Whdurs(:,:,i) = WTWhdurhist.Whdurs(:,:,i)./max(WTWhdurhist.Whdurs(:,:,i));
    WTWhdurhist.InterWhdurs(:,:,i) = WTWhdurhist.InterWhdurs(:,:,i)./max(WTWhdurhist.InterWhdurs(:,:,i));
end
WTEMGPSD = WTBehavior.EMGPSD;

KOEMGhist = KOBehavior.EMGhist;
KOWhdurhist = KOBehavior.Whdurhist;
for i = 1:size(KOWhdurhist.Whdurs,3)
        KOEMGhist.logcounts(:,:,i) = KOEMGhist.logcounts(:,:,i)./max(KOEMGhist.logcounts(:,:,i));
    KOWhdurhist.Whdurs(:,:,i) = KOWhdurhist.Whdurs(:,:,i)./max(KOWhdurhist.Whdurs(:,:,i));
    KOWhdurhist.InterWhdurs(:,:,i) = KOWhdurhist.InterWhdurs(:,:,i)./max(KOWhdurhist.InterWhdurs(:,:,i));
end
KOEMGPSD = KOBehavior.EMGPSD;

%% FIGURE 2:
figure;
subplot(2,2,1);
bar(WTEMGhist.logbins(:,:,1),...
    nanmean(WTEMGhist.logcounts,3),...
    'facecolor','k','facealpha',0.85); hold on;
h1 = errorbar(WTEMGhist.logbins(:,:,1),...
    nanmean(WTEMGhist.logcounts,3),...
    nanstd(WTEMGhist.logcounts,0,3),'k.');
bar(KOEMGhist.logbins(:,:,1),...
    nanmean(KOEMGhist.logcounts,3),...
    'facecolor','r','facealpha',0.65); hold on;
h2 = errorbar(KOEMGhist.logbins(:,:,1),...
    nanmean(KOEMGhist.logcounts,3),...
    nanstd(KOEMGhist.logcounts,0,3),'r.');
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
LogScale('x',10);
axis tight
ylim([0 1]);
xlabel('EMG (norm.)'); ylabel('occupancy (norm.)');
legend({'floxed M1/M3','Emx1;M1/M3'},'location','northeast');

subplot(2,2,2);
shadederror(WTWhdurhist.bins(:,:,1),...
    nanmean(WTWhdurhist.Whdurs,3)',...
    nanstd(WTWhdurhist.Whdurs,0,3)','k','none',1,'k');
shadederror(KOWhdurhist.bins(:,:,1),...
    nanmean(KOWhdurhist.Whdurs,3)',...
    nanstd(KOWhdurhist.Whdurs,0,3)','r','none',1,'r');
LogScale('x',10)
%legend({'floxed M1/M3','Emx1;M1/M3'},'location','northeast');
xlabel('duration (s)'); ylabel('counts (norm.)')
title('Wh durations')

subplot(2,2,4);
shadederror(WTWhdurhist.bins(:,:,1),...
    nanmean(WTWhdurhist.InterWhdurs,3)',...
    nanstd(WTWhdurhist.InterWhdurs,0,3)','k','none',1,'k');
shadederror(KOWhdurhist.bins(:,:,1),...
    nanmean(KOWhdurhist.InterWhdurs,3)',...
    nanstd(KOWhdurhist.InterWhdurs,0,3)','r','none',1,'r');
LogScale('x',10)
%legend({'floxed M1/M3','Emx1;M1/M3'},'location','northeast');
xlabel('duration (s)'); ylabel('counts (norm.)')
title('NWh durations')

subplot(2,2,3);
shadederror(log10(WTEMGPSD.freqs(:,:,1)),nanmean(WTEMGPSD.psd,3)',...
    nanstd(WTEMGPSD.psd,0,3)','k','none',1,'k');
shadederror(log10(KOEMGPSD.freqs(:,:,1)),nanmean(KOEMGPSD.psd,3)',...
    nanstd(KOEMGPSD.psd,0,3)','r','none',1,'r');
LogScale('x',10)
xlabel('f (Hz)'); ylabel('power (dB)')
axis tight
ylim([-1 2.5]);
title('EMG power spectrum');

NiceSave('EMGStats',figfolder,baseName)

%%
WTpupildynamicsEMG = WTBehavior.pupildynamicsEMG;
WTpupilphaseEMG = WTBehavior.pupilphaseEMG;

KOpupildynamicsEMG = KOBehavior.pupildynamicsEMG;
KOpupilphaseEMG = KOBehavior.pupilphaseEMG;

%% FIGURE 3:
emgcolor = [1 1 1;makeColorMap([0.5 0.5 0.5],[0 0 0.8])];

figure;
subplot(4,4,1);
imagesc(WTpupildynamicsEMG.bins(:,:,1),WTpupildynamicsEMG.bins(:,:,1),...
    nanmean(WTpupildynamicsEMG.mean,3)'); hold on;
colormap(gca,emgcolor); colorbar; caxis([-0.75 0.5])
axis xy
LogScale('c',10);
xlim([-0.5 0.5]); ylim([-0.5 0.5])
ylabel('dP/dt')
title('WT EMG')
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color','r')

subplot(4,4,2);
imagesc(KOpupildynamicsEMG.bins(:,:,1),KOpupildynamicsEMG.bins(:,:,1),...
    nanmean(KOpupildynamicsEMG.mean,3)'); hold on;
colormap(gca,emgcolor); colorbar; caxis([-0.75 0.5])
axis xy
LogScale('c',10);
xlim([-0.5 0.5]); ylim([-0.5 0.5])
ylabel('dP/dt')
title('KO EMG')
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color','r')

subplot(4,4,5);
imagesc(WTpupildynamicsEMG.bins(:,:,1),WTpupildynamicsEMG.bins(:,:,1),...
    nanmean(WTpupildynamicsEMG.pWhisk,3)'); hold on;
colormap(gca,emgcolor)
axis xy
xlim([-0.5 0.5]);ylim([-0.5 0.5])
caxis([-0.01 1])
colorbar
ylabel('dP/dt')
title('p(Wh)')
hold on
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color','r')

subplot(4,4,6);
imagesc(KOpupildynamicsEMG.bins(:,:,1),KOpupildynamicsEMG.bins(:,:,1),...
    nanmean(KOpupildynamicsEMG.pWhisk,3)'); hold on;
colormap(gca,emgcolor)
axis xy
xlim([-0.5 0.5]);ylim([-0.5 0.5])
caxis([-0.01 1])
colorbar
ylabel('dP/dt')
title('p(Wh)')
hold on
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color','r')

subplot(4,4,9);
imagesc(WTpupildynamicsEMG.bins(:,:,1),WTpupildynamicsEMG.bins(:,:,1),...
    nanmean(WTpupildynamicsEMG.Whstartrate,3)'); hold on;
colormap(gca,emgcolor)
axis xy
xlim([-0.5 0.5]);ylim([-0.5 0.5])
colorbar
ylabel('dP/dt')
title('p(Wh onset)')
hold on
caxis([-0.1 1])
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color','r')

subplot(4,4,10);
imagesc(KOpupildynamicsEMG.bins(:,:,1),KOpupildynamicsEMG.bins(:,:,1),...
    nanmean(KOpupildynamicsEMG.Whstartrate,3)'); hold on;
colormap(gca,emgcolor)
axis xy
xlim([-0.5 0.5]);ylim([-0.5 0.5])
colorbar
ylabel('dP/dt')
title('p(Wh onset)')
hold on
caxis([-0.1 1])
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color','r')

subplot(4,4,13);
imagesc(WTpupildynamicsEMG.bins(:,:,1),WTpupildynamicsEMG.bins(:,:,1),...
    log10(nanmean(WTpupildynamicsEMG.meanWhdur,3))'); hold on;
axis xy
xlim([-0.5 0.5]);ylim([-0.5 0.5])
colorbar
xlabel('Pupil Diameter');ylabel('dP/dt')
title('Whisk Duration')
hold on
caxis([-1 0.5])
LogScale('c',10)
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color','r')

subplot(4,4,14);
imagesc(KOpupildynamicsEMG.bins(:,:,1),KOpupildynamicsEMG.bins(:,:,1),...
    log10(nanmean(KOpupildynamicsEMG.meanWhdur,3))'); hold on;
axis xy
xlim([-0.5 0.5]);ylim([-0.5 0.5])
colorbar
xlabel('Pupil Diameter');ylabel('dP/dt')
title('Whisk Duration')
hold on
caxis([-1 0.5])
LogScale('c',10)
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color','r')

% THE other side...
subplot(4,4,3);
imagesc(WTpupilphaseEMG.bins(:,:,1),WTpupilphaseEMG.bins(:,:,1),...
    nanmean(WTpupilphaseEMG.mean,3)); hold on;
imagesc(WTpupilphaseEMG.bins(:,:,1)+2*pi,WTpupilphaseEMG.bins(:,:,1),...
    nanmean(WTpupilphaseEMG.mean,3)); hold on;
colormap(gca,emgcolor)
axis xy
ylim([-2 0.25]);xlim([-pi 3*pi])
colorbar
ylabel('power (dB)')
title('WT EMG')
hold on
caxis([-0.75 0.5])
LogScale('c',10)
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color',0.5.*[1 1 1])

subplot(4,4,4);
imagesc(KOpupilphaseEMG.bins(:,:,1),KOpupilphaseEMG.bins(:,:,1),...
    nanmean(KOpupilphaseEMG.mean,3)); hold on;
imagesc(KOpupilphaseEMG.bins(:,:,1)+2*pi,KOpupilphaseEMG.bins(:,:,1),...
    nanmean(KOpupilphaseEMG.mean,3)); hold on;
colormap(gca,emgcolor)
axis xy
ylim([-2 0.25]);xlim([-pi 3*pi])
colorbar
ylabel('power (dB)')
title('KO EMG')
hold on
caxis([-0.75 0.5])
LogScale('c',10)
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color',0.5.*[1 1 1]);

subplot(4,4,7);
imagesc(WTpupilphaseEMG.bins(:,:,1),WTpupilphaseEMG.bins(:,:,1),...
    nanmean(WTpupilphaseEMG.pWhisk,3)); hold on;
imagesc(WTpupilphaseEMG.bins(:,:,1)+2*pi,WTpupilphaseEMG.bins(:,:,1),...
    nanmean(WTpupilphaseEMG.pWhisk,3)); hold on;
colormap(gca,emgcolor)
axis xy
ylim([-2 0.25]);xlim([-pi 3*pi])
colorbar
ylabel('power (dB)')
title('p(Wh)')
hold on
caxis([-0.01 1])
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color',0.5.*[1 1 1])

subplot(4,4,8);
imagesc(KOpupilphaseEMG.bins(:,:,1),KOpupilphaseEMG.bins(:,:,1),...
    nanmean(KOpupilphaseEMG.pWhisk,3)); hold on;
imagesc(KOpupilphaseEMG.bins(:,:,1)+2*pi,KOpupilphaseEMG.bins(:,:,1),...
    nanmean(KOpupilphaseEMG.pWhisk,3)); hold on;
colormap(gca,emgcolor)
axis xy
ylim([-2 0.25]);xlim([-pi 3*pi])
colorbar
ylabel('power (dB)')
title('p(Wh)')
hold on
caxis([-0.01 1])
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color',0.5.*[1 1 1])

subplot(4,4,11);
imagesc(WTpupilphaseEMG.bins(:,:,1),WTpupilphaseEMG.bins(:,:,1),...
    nanmean(WTpupilphaseEMG.Whstartrate,3)); hold on;
imagesc(WTpupilphaseEMG.bins(:,:,1)+2*pi,WTpupilphaseEMG.bins(:,:,1),...
    nanmean(WTpupilphaseEMG.Whstartrate,3)); hold on;
colormap(gca,emgcolor)
axis xy
ylim([-2 0.25]);xlim([-pi 3*pi])
colorbar
ylabel('power (dB)')
title('p(Wh onset)')
hold on
caxis([-0.1 1])
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color',0.5.*[1 1 1])

subplot(4,4,12);
imagesc(KOpupilphaseEMG.bins(:,:,1),KOpupilphaseEMG.bins(:,:,1),...
    nanmean(KOpupilphaseEMG.Whstartrate,3)); hold on;
imagesc(KOpupilphaseEMG.bins(:,:,1)+2*pi,KOpupilphaseEMG.bins(:,:,1),...
    nanmean(KOpupilphaseEMG.Whstartrate,3)); hold on;
colormap(gca,emgcolor)
axis xy
ylim([-2 0.25]);xlim([-pi 3*pi])
colorbar
ylabel('power (dB)')
title('p(Wh onset)')
hold on
caxis([-0.1 1])
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color',0.5.*[1 1 1])

subplot(4,4,15);
imagesc(WTpupilphaseEMG.bins(:,:,1),WTpupilphaseEMG.bins(:,:,1),...
    log10(nanmean(WTpupilphaseEMG.meanWhdur,3))); hold on;
imagesc(WTpupilphaseEMG.bins(:,:,1)+2*pi,WTpupilphaseEMG.bins(:,:,1),...
    log10(nanmean(WTpupilphaseEMG.meanWhdur,3))); hold on;
axis xy
ylim([-2 0.25]);xlim([-pi 3*pi])
colorbar
xlabel('pupil phase');ylabel('power (dB)')
title('Wh duration')
hold on
caxis([-1 0.5])
LogScale('c',10)
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color',0.5.*[1 1 1])

subplot(4,4,16);
imagesc(KOpupilphaseEMG.bins(:,:,1),KOpupilphaseEMG.bins(:,:,1),...
    log10(nanmean(KOpupilphaseEMG.meanWhdur,3))); hold on;
imagesc(KOpupilphaseEMG.bins(:,:,1)+2*pi,KOpupilphaseEMG.bins(:,:,1),...
    log10(nanmean(KOpupilphaseEMG.meanWhdur,3))); hold on;
axis xy
ylim([-2 0.25]);xlim([-pi 3*pi])
colorbar
xlabel('pupil phase');ylabel('power (dB)')
title('Wh duration')
hold on
caxis([-1 0.5])
LogScale('c',10)
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color',0.5.*[1 1 1])

NiceSave('EMG_Pupil_Space',figfolder,baseName);

%%

WTcoupling = WTBehavior.WPcoupling;
KOcoupling = KOBehavior.WPcoupling;

%% FIGURE 4:
figure;
shadederror(log10(WTcoupling.freqs(:,:,1)),...
    nanmean(WTcoupling.coupling,3),...
    nanstd(WTcoupling.coupling,0,3),'k','none',1,'k');
shadederror(log10(KOcoupling.freqs(:,:,1)),...
    nanmean(KOcoupling.coupling,3),...
    nanstd(KOcoupling.coupling,0,3),'r','none',1,'r');
LogScale('x',10);
xlabel('f (Hz)');
axis tight
title('Pupil-EMG phase coupling')

NiceSave('Pupil_EMG_Phase_coupling',figfolder,baseName)
