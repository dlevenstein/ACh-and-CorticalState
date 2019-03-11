%%
basePath = pwd;
baseName = 'EM1M3';

figfolder = fullfile(basePath,'SummaryFigures');

load(fullfile(basePath,['WT_EM1M3.PSSBehaviorAnalysis.mat']));
WTPSSBehavior = groupPSSBehavior;
load(fullfile(basePath,['KO_EM1M3.PSSBehaviorAnalysis.mat']));
KOPSSBehavior = groupPSSBehavior;

%%

WTPSShist = WTPSSBehavior.PSShist;
WTpupcyclePSS = WTPSSBehavior.pupcyclePSS;
for i = 1:size(WTPSShist.hist,3)
    WTPSShist.hist(:,:,i) = WTPSShist.hist(:,:,i)./max(WTPSShist.hist(:,:,i));
end

KOPSShist = KOPSSBehavior.PSShist;
KOpupcyclePSS = KOPSSBehavior.pupcyclePSS;
for i = 1:size(KOPSShist.hist,3)
    KOPSShist.hist(:,:,i) = KOPSShist.hist(:,:,i)./max(KOPSShist.hist(:,:,i));
end

%% FIGURE:
figure;
subplot(2,2,1)
bar(WTPSShist.bins(:,:,1),nanmean(WTPSShist.hist,3),...
    'facecolor','k','facealpha',0.65); hold on;
errorbar(WTPSShist.bins(:,:,1),...
    nanmean(WTPSShist.hist,3),...
    nanstd(WTPSShist.hist,0,3)./sqrt(size(WTPSShist.hist,3)),'k.');
bar(KOPSShist.bins(:,:,1),nanmean(KOPSShist.hist,3),...
        'facecolor','r','facealpha',0.65);
errorbar(KOPSShist.bins(:,:,1),...
    nanmean(KOPSShist.hist,3),...
    nanstd(KOPSShist.hist,0,3)./sqrt(size(KOPSShist.hist,3)),'r.');
xlabel('PSS (au)'); ylabel('occupancy (norm.)');
axis tight
xlim([-3.5 0]); ylim([0 1])

cmin = min([min(min(nanmean(WTpupcyclePSS.meanPSS,3)))...
    min(min(nanmean(KOpupcyclePSS.meanPSS,3)))]);
cmax = max([max(max(nanmean(WTpupcyclePSS.meanPSS,3)))...
    max(max(nanmean(KOpupcyclePSS.meanPSS,3)))]);

subplot(2,2,3)
h = imagesc(WTpupcyclePSS.bincenters(:,:,1),...
    WTpupcyclePSS.bincenters(:,:,1),nanmean(WTpupcyclePSS.meanPSS,3)');
set(h,'AlphaData',~isnan(nanmean(WTpupcyclePSS.meanPSS,3)'));
hold on
plot(WTpupcyclePSS.bincenters(:,[1 end],1),[0 0],'k--')
LogScale('x',10)
axis xy
axis square
colorbar; caxis([cmin cmax]);
xlabel('Pupil area (med^-^1)')
ylabel('dP/dt')
LogScale('x',10)
title('WT PSS')

subplot(2,2,4)
h = imagesc(KOpupcyclePSS.bincenters(:,:,1),...
    KOpupcyclePSS.bincenters(:,:,1),nanmean(KOpupcyclePSS.meanPSS,3)');
set(h,'AlphaData',~isnan(nanmean(KOpupcyclePSS.meanPSS,3)'));
hold on
plot(KOpupcyclePSS.bincenters(:,[1 end],1),[0 0],'k--')
LogScale('x',10)
axis xy
axis square
colorbar; caxis([cmin cmax]);
xlabel('Pupil area (med^-^1)')
ylabel('dP/dt')
LogScale('x',10)
title('KO PSS')

NiceSave('PSS_Phase_dPdt',figfolder,baseName)

%%
WTpupildist = WTPSSBehavior.pupildist;
KOpupildist = KOPSSBehavior.pupildist;

%% FIGURE:
cosx = linspace(-pi,3*pi,100);

cmin = min([min(min(nanmean(WTpupildist.joint,3)))...
    min(min(nanmean(KOpupildist.joint,3)))]);
cmax = max([max(max(nanmean(WTpupildist.joint,3)))...
    max(max(nanmean(KOpupildist.joint,3)))]);
figure;
subplot(2,2,1);
imagesc(WTpupildist.bins{:,1,1},WTpupildist.bins{:,2,1},nanmean(WTpupildist.joint,3)')
hold on
imagesc(WTpupildist.bins{:,1,1}+2*pi,WTpupildist.bins{:,2,1},nanmean(WTpupildist.joint,3)')
plot(cosx,cos(cosx)./3-1.1,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil Phase');ylabel('Pupil Power')
title('WT P(Pupil power,phase)')

subplot(2,2,2);
imagesc(KOpupildist.bins{:,1,1},KOpupildist.bins{:,2,1},nanmean(KOpupildist.joint,3)')
hold on
imagesc(KOpupildist.bins{:,1,1}+2*pi,KOpupildist.bins{:,2,1},nanmean(KOpupildist.joint,3)')
plot(cosx,cos(cosx)./3-1.1,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil Phase');ylabel('Pupil Power')
title('KO P(Pupil power,phase)')

cmin = min([min(min(nanmean(WTpupildist.conditional,3)))...
    min(min(nanmean(KOpupildist.conditional,3)))]);
cmax = max([max(max(nanmean(WTpupildist.conditional,3)))...
    max(max(nanmean(KOpupildist.conditional,3)))]);
subplot(2,2,3);
imagesc(WTpupildist.bins{:,1,1},WTpupildist.bins{:,2,1},nanmean(WTpupildist.conditional,3)')
hold on
imagesc(WTpupildist.bins{:,1,1}+2*pi,WTpupildist.bins{:,2,1},nanmean(WTpupildist.conditional,3)')
plot(cosx,cos(cosx)./3-1.1,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil Phase');ylabel('Pupil Power')
title('P(Pupil power|phase)')

subplot(2,2,4);
imagesc(KOpupildist.bins{:,1,1},KOpupildist.bins{:,2,1},nanmean(KOpupildist.conditional,3)')
hold on
imagesc(KOpupildist.bins{:,1,1}+2*pi,KOpupildist.bins{:,2,1},nanmean(KOpupildist.conditional,3)')
plot(cosx,cos(cosx)./3-1.1,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil Phase');ylabel('Pupil Power')
title('P(Pupil power|phase)')

NiceSave('PSS_PhasePow',figfolder,baseName)

%%
WTpupilPSSdist = WTPSSBehavior.pupilPSSdist;
KOpupilPSSdist = KOPSSBehavior.pupilPSSdist;

%% FIGURE:
cmin = min([min(min(nanmean(WTpupilPSSdist.joint,3)))...
    min(min(nanmean(KOpupilPSSdist.joint,3)))]);
cmax = max([max(max(nanmean(WTpupilPSSdist.joint,3)))...
    max(max(nanmean(KOpupilPSSdist.joint,3)))]);

figure;
subplot(2,2,1);
imagesc(WTpupilPSSdist.bins{:,1,1},WTpupilPSSdist.bins{:,2,1},nanmean(WTpupilPSSdist.joint,3)')
hold on
imagesc(WTpupilPSSdist.bins{:,1,1}+2*pi,WTpupilPSSdist.bins{:,2,1},nanmean(WTpupilPSSdist.joint,3)')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('PSS')
title('WT P(PSS,phase)')

subplot(2,2,2);
imagesc(KOpupilPSSdist.bins{:,1,1},KOpupilPSSdist.bins{:,2,1},nanmean(KOpupilPSSdist.joint,3)')
hold on
imagesc(KOpupilPSSdist.bins{:,1,1}+2*pi,KOpupilPSSdist.bins{:,2,1},nanmean(KOpupilPSSdist.joint,3)')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil Phase');ylabel('PSS')
title('KO P(PSS,phase)')

cmin = min([min(min(nanmean(WTpupilPSSdist.conditional,3)))...
    min(min(nanmean(KOpupilPSSdist.conditional,3)))]);
cmax = max([max(max(nanmean(WTpupilPSSdist.conditional,3)))...
    max(max(nanmean(KOpupilPSSdist.conditional,3)))]);

subplot(2,2,3);
imagesc(WTpupilPSSdist.bins{:,1,1},WTpupilPSSdist.bins{:,2,1},...
    nanmean(WTpupilPSSdist.conditional,3)')
hold on
imagesc(WTpupilPSSdist.bins{:,1,1}+2*pi,WTpupilPSSdist.bins{:,2,1},...
    nanmean(WTpupilPSSdist.conditional,3)')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('PSS')
title('P(PSS|phase)')

subplot(2,2,4);
imagesc(KOpupilPSSdist.bins{:,1,1},KOpupilPSSdist.bins{:,2,1},...
    nanmean(KOpupilPSSdist.conditional,3)')
hold on
imagesc(KOpupilPSSdist.bins{:,1,1}+2*pi,KOpupilPSSdist.bins{:,2,1},...
    nanmean(KOpupilPSSdist.conditional,3)')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('PSS')
title('P(PSS|phase)')

NiceSave('PSSbyPupilPhase',figfolder,baseName)

%% FIGURE:
cmin = min([min(min(nanmean(WTpupilPSSdist.counts_high,3)))...
    min(min(nanmean(KOpupilPSSdist.counts_high,3)))]);
cmax = max([max(max(nanmean(WTpupilPSSdist.counts_high,3)))...
    max(max(nanmean(KOpupilPSSdist.counts_high,3)))]);

figure;
subplot(2,2,1);
imagesc(WTpupilPSSdist.bins{:,1,1},WTpupilPSSdist.bins{:,2,1},...
    nanmean(WTpupilPSSdist.counts_high,3)'); hold on;
imagesc(WTpupilPSSdist.bins{:,1,1}+2*pi,WTpupilPSSdist.bins{:,2,1},...
    nanmean(WTpupilPSSdist.counts_high,3)')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('PSS')
title('WT >median pupil');

subplot(2,2,2);
imagesc(KOpupilPSSdist.bins{:,1,1},KOpupilPSSdist.bins{:,2,1},...
    nanmean(KOpupilPSSdist.counts_high,3)'); hold on;
imagesc(KOpupilPSSdist.bins{:,1,1}+2*pi,KOpupilPSSdist.bins{:,2,1},...
    nanmean(KOpupilPSSdist.counts_high,3)')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('PSS')
title('KO >median pupil');

cmin = min([min(min(nanmean(WTpupilPSSdist.counts_low,3)))...
    min(min(nanmean(KOpupilPSSdist.counts_low,3)))]);
cmax = max([max(max(nanmean(WTpupilPSSdist.counts_low,3)))...
    max(max(nanmean(KOpupilPSSdist.counts_low,3)))]);

subplot(2,2,3);
imagesc(WTpupilPSSdist.bins{:,1,1},WTpupilPSSdist.bins{:,2,1},...
    nanmean(WTpupilPSSdist.counts_low,3)')
hold on
imagesc(WTpupilPSSdist.bins{:,1,1}+2*pi,WTpupilPSSdist.bins{:,2,1},...
    nanmean(WTpupilPSSdist.counts_low,3)')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil Phase');ylabel('PSS')
title('<median pupil');

subplot(2,2,4);
imagesc(KOpupilPSSdist.bins{:,1,1},KOpupilPSSdist.bins{:,2,1},...
    nanmean(KOpupilPSSdist.counts_low,3)')
hold on
imagesc(KOpupilPSSdist.bins{:,1,1}+2*pi,KOpupilPSSdist.bins{:,2,1},...
    nanmean(KOpupilPSSdist.counts_low,3)')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil Phase');ylabel('PSS')
title('<median pupil');

NiceSave('PSSlohiPupilPhase',figfolder,baseName)

%%
WTpupilEMGdist = WTPSSBehavior.pupilEMGdist;
KOpupilEMGdist = KOPSSBehavior.pupilEMGdist;

%% FIGURE:
cosx = linspace(-pi,3*pi,100);
cmin = min([min(min(nanmean(WTpupilEMGdist.joint,3)))...
    min(min(nanmean(KOpupilEMGdist.joint,3)))]);
cmax = max([max(max(nanmean(WTpupilEMGdist.joint,3)))...
    max(max(nanmean(KOpupilEMGdist.joint,3)))]);

figure;
subplot(2,2,1);
imagesc(WTpupilEMGdist.bins{:,1,1},WTpupilEMGdist.bins{:,2,1},...
    nanmean(WTpupilEMGdist.joint,3)'); hold on;
imagesc(WTpupilEMGdist.bins{:,1,1}+2*pi,WTpupilEMGdist.bins{:,2,1},...
    nanmean(WTpupilEMGdist.joint,3)')
plot(cosx,cos(cosx)./3-1.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('EMG')
title('WT P(EMG,phase)')

subplot(2,2,2);
imagesc(KOpupilEMGdist.bins{:,1,1},KOpupilEMGdist.bins{:,2,1},...
    nanmean(KOpupilEMGdist.joint,3)'); hold on;
imagesc(KOpupilEMGdist.bins{:,1,1}+2*pi,KOpupilEMGdist.bins{:,2,1},...
    nanmean(KOpupilEMGdist.joint,3)')
plot(cosx,cos(cosx)./3-1.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('EMG')
title('KO P(EMG,phase)')

cmin = min([min(min(nanmean(WTpupilEMGdist.conditional,3)))...
    min(min(nanmean(KOpupilEMGdist.conditional,3)))]);
cmax = max([max(max(nanmean(WTpupilEMGdist.conditional,3)))...
    max(max(nanmean(KOpupilEMGdist.conditional,3)))]);

subplot(2,2,3);
imagesc(WTpupilEMGdist.bins{:,1,1},WTpupilEMGdist.bins{:,2,1},...
    nanmean(WTpupilEMGdist.conditional,3)'); hold on;
imagesc(WTpupilEMGdist.bins{:,1,1}+2*pi,WTpupilEMGdist.bins{:,2,1},...
    nanmean(WTpupilEMGdist.conditional,3)')
plot(cosx,cos(cosx)./3-1.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('EMG')
title('P(EMG|phase)')

subplot(2,2,4);
imagesc(KOpupilEMGdist.bins{:,1,1},KOpupilEMGdist.bins{:,2,1},...
    nanmean(KOpupilEMGdist.conditional,3)'); hold on;
imagesc(KOpupilEMGdist.bins{:,1,1}+2*pi,KOpupilEMGdist.bins{:,2,1},...
    nanmean(KOpupilEMGdist.conditional,3)')
plot(cosx,cos(cosx)./3-1.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('EMG')
title('P(EMG|phase)')

NiceSave('EMGbyPupilPhase',figfolder,baseName)

%% FIGURE:
cmin = min([min(min(nanmean(WTpupilEMGdist.counts_high,3)))...
    min(min(nanmean(KOpupilEMGdist.counts_high,3)))]);
cmax = max([max(max(nanmean(WTpupilEMGdist.counts_high,3)))...
    max(max(nanmean(KOpupilEMGdist.counts_high,3)))]);

figure;
subplot(2,2,1);
imagesc(WTpupilEMGdist.bins{:,1,1},WTpupilEMGdist.bins{:,2,1},...
    nanmean(WTpupilEMGdist.counts_high,3)'); hold on;
imagesc(WTpupilEMGdist.bins{:,1,1}+2*pi,WTpupilEMGdist.bins{:,2,1},...
    nanmean(WTpupilEMGdist.counts_high,3)')
plot(cosx,cos(cosx)./3-1.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('EMG')
title('WT >median pupil');

subplot(2,2,2);
imagesc(KOpupilEMGdist.bins{:,1,1},KOpupilEMGdist.bins{:,2,1},...
    nanmean(KOpupilEMGdist.counts_high,3)'); hold on;
imagesc(KOpupilEMGdist.bins{:,1,1}+2*pi,KOpupilEMGdist.bins{:,2,1},...
    nanmean(KOpupilEMGdist.counts_high,3)')
plot(cosx,cos(cosx)./3-1.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('EMG')
title('KO >median pupil');

cmin = min([min(min(nanmean(WTpupilEMGdist.counts_low,3)))...
    min(min(nanmean(KOpupilEMGdist.counts_low,3)))]);
cmax = max([max(max(nanmean(WTpupilEMGdist.counts_low,3)))...
    max(max(nanmean(KOpupilEMGdist.counts_low,3)))]);

subplot(2,2,3);
imagesc(WTpupilEMGdist.bins{:,1,1},WTpupilEMGdist.bins{:,2,1},...
    nanmean(WTpupilEMGdist.counts_low,3)'); hold on;
imagesc(WTpupilEMGdist.bins{:,1,1}+2*pi,WTpupilEMGdist.bins{:,2,1},...
    nanmean(WTpupilEMGdist.counts_low,3)')
plot(cosx,cos(cosx)./3-1.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('EMG')
title('<median pupil');

subplot(2,2,4);
imagesc(KOpupilEMGdist.bins{:,1,1},KOpupilEMGdist.bins{:,2,1},...
    nanmean(KOpupilEMGdist.counts_low,3)'); hold on;
imagesc(KOpupilEMGdist.bins{:,1,1}+2*pi,KOpupilEMGdist.bins{:,2,1},...
    nanmean(KOpupilEMGdist.counts_low,3)')
plot(cosx,cos(cosx)./3-1.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('EMG')
title('<median pupil');

NiceSave('EMGlohiPupilPhase',figfolder,baseName)

%%
WTPSSEMGdist = WTPSSBehavior.PSSEMGdist;
KOPSSEMGdist = KOPSSBehavior.PSSEMGdist;

%% FIGURE:
cmin = min([min(min(nanmean(WTPSSEMGdist.joint,3)))...
    min(min(nanmean(KOPSSEMGdist.joint,3)))]);
cmax = max([max(max(nanmean(WTPSSEMGdist.joint,3)))...
    max(max(nanmean(KOPSSEMGdist.joint,3)))]);

figure;
subplot(2,2,1);
imagesc(WTPSSEMGdist.edges{:,1,1},WTPSSEMGdist.edges{:,2,1},...
    nanmean(WTPSSEMGdist.joint,3)'); hold on;
axis xy
ylim([-3 -0.5])
caxis([cmin cmax]);
xlabel('EMG');ylabel('PSS')
title('WT P(EMG,PSS)')

subplot(2,2,2);
imagesc(KOPSSEMGdist.edges{:,1,1},KOPSSEMGdist.edges{:,2,1},...
    nanmean(KOPSSEMGdist.joint,3)'); hold on;
axis xy
ylim([-3 -0.5])
caxis([cmin cmax]);
xlabel('EMG');ylabel('PSS')
title('KO P(EMG,PSS)')

cmin = min([min(min(nanmean(WTPSSEMGdist.conditional,3)))...
    min(min(nanmean(KOPSSEMGdist.conditional,3)))]);
cmax = max([max(max(nanmean(WTPSSEMGdist.conditional,3)))...
    max(max(nanmean(KOPSSEMGdist.conditional,3)))]);

subplot(2,2,3);
imagesc(WTPSSEMGdist.edges{:,1,1},WTPSSEMGdist.edges{:,2,1},...
    nanmean(WTPSSEMGdist.conditional,3)'); hold on;
axis xy
ylim([-3 -0.5])
caxis([cmin cmax]);
xlabel('EMG');ylabel('PSS')
title('P(EMG|PSS)')

subplot(2,2,4);
imagesc(KOPSSEMGdist.edges{:,1,1},KOPSSEMGdist.edges{:,2,1},...
    nanmean(KOPSSEMGdist.conditional,3)'); hold on;
axis xy
ylim([-3 -0.5])
caxis([cmin cmax]);
xlabel('EMG');ylabel('PSS')
title('P(EMG|PSS)')

NiceSave('EMGPSS',figfolder,baseName)

%%
WTphasePETH = WTPSSBehavior.phasePETH;
KOphasePETH = KOPSSBehavior.phasePETH;

%% FIGURE:
cmin = min([min(min(nanmean(WTphasePETH.high.mean,3)))...
    min(min(nanmean(KOphasePETH.high.mean,3)))]);
cmax = max([max(max(nanmean(WTphasePETH.high.mean,3)))...
    max(max(nanmean(KOphasePETH.high.mean,3)))]);

figure;
subplot(2,2,1);
imagesc(WTphasePETH.high.bincenters(:,:,1),WTphasePETH.high.bincenters(:,:,1),...
    nanmean(WTphasePETH.high.mean,3)'); hold on;
imagesc(WTphasePETH.high.bincenters(:,:,1),WTphasePETH.high.bincenters(:,:,1)+2*pi,...
    nanmean(WTphasePETH.high.mean,3)')
plot(cos(cosx),cosx,'w','linewidth',2)
plot([0 0],[-pi 3*pi],'r')
colorbar
axis xy
caxis([cmin cmax])
xlim([-1 4]);ylim([-pi 3*pi])
xlabel('t (s, aligned to Wh Onset)');ylabel('Pupil phase')
title('WT >median pupil')

subplot(2,2,2);
imagesc(KOphasePETH.high.bincenters(:,:,1),KOphasePETH.high.bincenters(:,:,1),...
    nanmean(KOphasePETH.high.mean,3)'); hold on;
imagesc(KOphasePETH.high.bincenters(:,:,1),KOphasePETH.high.bincenters(:,:,1)+2*pi,...
    nanmean(KOphasePETH.high.mean,3)')
plot(cos(cosx),cosx,'w','linewidth',2)
plot([0 0],[-pi 3*pi],'r')
colorbar
axis xy
caxis([cmin cmax])
xlim([-1 4]);ylim([-pi 3*pi])
xlabel('t (s, aligned to Wh Onset)');ylabel('Pupil phase')
title('KO >median pupil')

cmin = min([min(min(nanmean(WTphasePETH.low.mean,3)))...
    min(min(nanmean(KOphasePETH.low.mean,3)))]);
cmax = max([max(max(nanmean(WTphasePETH.low.mean,3)))...
    max(max(nanmean(KOphasePETH.low.mean,3)))]);

subplot(2,2,3);
imagesc(WTphasePETH.low.bincenters(:,:,1),WTphasePETH.low.bincenters(:,:,1),...
    nanmean(WTphasePETH.low.mean,3)'); hold on
imagesc(WTphasePETH.low.bincenters(:,:,1),WTphasePETH.low.bincenters(:,:,1)+2*pi,...
    nanmean(WTphasePETH.low.mean,3)');
plot(cos(cosx),cosx,'w','linewidth',2)
plot([0 0],[-pi 3*pi],'r')
colorbar
axis xy
caxis([cmin cmax])
xlim([-1 4]);ylim([-pi 3*pi])
xlabel('t (s, aligned to Wh Onset)');ylabel('Pupil phase')
title('<median pupil')

subplot(2,2,4);
imagesc(KOphasePETH.low.bincenters(:,:,1),KOphasePETH.low.bincenters(:,:,1),...
    nanmean(KOphasePETH.low.mean,3)'); hold on
imagesc(KOphasePETH.low.bincenters(:,:,1),KOphasePETH.low.bincenters(:,:,1)+2*pi,...
    nanmean(KOphasePETH.low.mean,3)');
plot(cos(cosx),cosx,'w','linewidth',2)
plot([0 0],[-pi 3*pi],'r')
colorbar
axis xy
caxis([cmin cmax])
xlim([-1 4]);ylim([-pi 3*pi])
xlabel('t (s, aligned to Wh Onset)');ylabel('Pupil phase')
title('<median pupil')

NiceSave('PETHbyPhase_Whaligned',figfolder,baseName)
