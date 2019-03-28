%%
basePath = pwd;
baseName = 'EM1M3';

figfolder = fullfile(basePath,'SummaryFigures');
if (~exist(figfolder,'dir'))
    mkdir(figfolder)
end

load(fullfile(basePath,['WT_EM1M3.PSSBehaviorAnalysis2.mat']));
WTPSSBehavior = groupPSSBehavior;
WTpupsorts = pupsorts;
WTpuplockedPSS = puplockedPSS;
WTwhisksorts = whisksorts;
WTtimelockedPSS = timelockedPSS;

load(fullfile(basePath,['KO_EM1M3.PSSBehaviorAnalysis2.mat']));
KOPSSBehavior = groupPSSBehavior;
KOpupsorts = pupsorts;
KOpuplockedPSS = puplockedPSS;
KOwhisksorts = whisksorts;
KOtimelockedPSS = timelockedPSS;

%% FIGURE:
upcolor = [1 1 1;makeColorMap([1 1 1],[0.8 0 0])];
downcolor = [1 1 1;makeColorMap([1 1 1],[0 0 0.8])];

figure;
subplot(2,5,1);
barh(WTPSSBehavior.UPdur.bins(:,:,1),nanmean(WTPSSBehavior.UPdur.hist,3),'r'); hold on;
barh(WTPSSBehavior.DOWNdur.bins(:,:,1),nanmean(WTPSSBehavior.DOWNdur.hist,3),'b'); hold on;
ylabel('UP/DOWN duration (s)');
set(gca,'Xdir','reverse');
LogScale('y',10);
xlabel('norm counts');
title('WT');

subplot(2,5,6);
barh(KOPSSBehavior.UPdur.bins(:,:,1),nanmean(KOPSSBehavior.UPdur.hist,3),'r'); hold on;
barh(KOPSSBehavior.DOWNdur.bins(:,:,1),nanmean(KOPSSBehavior.DOWNdur.hist,3),'b'); hold on;
ylabel('UP/DOWN duration (s)');
set(gca,'Xdir','reverse');
LogScale('y',10);
xlabel('norm counts');
title('KO');

subplot(2,5,2:3);
imagesc(WTPSSBehavior.PSSUPDOWNhist.bins{:,1,1},WTPSSBehavior.PSSUPDOWNhist.bins{:,2,1},...
    nanmean(WTPSSBehavior.PSSUPDOWNhist.UPcounts,3)'); hold on;
colormap(gca,upcolor)
%ColorbarWithAxis([min(min(PSSUPDOWNhist.UPcounts)) max(max(PSSUPDOWNhist.UPcounts))/4],['counts (au)'])
caxis([min(min(nanmean(WTPSSBehavior.PSSUPDOWNhist.UPcounts,3)))...
    max(max(nanmean(WTPSSBehavior.PSSUPDOWNhist.UPcounts,3)))/4]);
xlabel('PSS (au)'); %ylabel('UP duration (s)')
LogScale('y',10);
axis square
axis xy

subplot(2,5,4:5);
imagesc(WTPSSBehavior.PSSUPDOWNhist.bins{:,1,1},WTPSSBehavior.PSSUPDOWNhist.bins{:,2,1},...
    nanmean(WTPSSBehavior.PSSUPDOWNhist.DOWNcounts,3)'); hold on;
colormap(gca,downcolor)
%ColorbarWithAxis([min(min(PSSUPDOWNhist.DOWNcounts)) max(max(PSSUPDOWNhist.DOWNcounts))/4],['counts (au)'])
caxis([min(min(nanmean(WTPSSBehavior.PSSUPDOWNhist.DOWNcounts,3)))...
    max(max(nanmean(WTPSSBehavior.PSSUPDOWNhist.DOWNcounts,3)))/4]);
xlabel('PSS (au)'); %ylabel('DOWN duration (s)')
LogScale('y',10);
axis square
axis xy

subplot(2,5,7:8);
imagesc(KOPSSBehavior.PSSUPDOWNhist.bins{:,1,1},KOPSSBehavior.PSSUPDOWNhist.bins{:,2,1},...
    nanmean(KOPSSBehavior.PSSUPDOWNhist.UPcounts,3)'); hold on;
colormap(gca,upcolor)
%ColorbarWithAxis([min(min(PSSUPDOWNhist.UPcounts)) max(max(PSSUPDOWNhist.UPcounts))/4],['counts (au)'])
caxis([min(min(nanmean(KOPSSBehavior.PSSUPDOWNhist.UPcounts,3)))...
    max(max(nanmean(KOPSSBehavior.PSSUPDOWNhist.UPcounts,3)))/4]);
xlabel('PSS (au)'); %ylabel('UP duration (s)')
LogScale('y',10);
axis square
axis xy

subplot(2,5,9:10);
imagesc(KOPSSBehavior.PSSUPDOWNhist.bins{:,1,1},KOPSSBehavior.PSSUPDOWNhist.bins{:,2,1},...
    nanmean(KOPSSBehavior.PSSUPDOWNhist.DOWNcounts,3)'); hold on;
colormap(gca,downcolor)
%ColorbarWithAxis([min(min(PSSUPDOWNhist.DOWNcounts)) max(max(PSSUPDOWNhist.DOWNcounts))/4],['counts (au)'])
caxis([min(min(nanmean(KOPSSBehavior.PSSUPDOWNhist.DOWNcounts,3)))...
    max(max(nanmean(KOPSSBehavior.PSSUPDOWNhist.DOWNcounts,3)))/4]);
xlabel('PSS (au)'); %ylabel('DOWN duration (s)')
LogScale('y',10);
axis square
axis xy

NiceSave('PSS_SlowWaves',figfolder,baseName)

%% FIGURE:
figure;
subplot(2,4,1)
bar(WTPSSBehavior.PSShist.bins(:,:,1),nanmean(WTPSSBehavior.PSShist.hist,3),...
    'facecolor','k','facealpha',0.65); hold on;
errorbar(WTPSSBehavior.PSShist.bins(:,:,1),...
    nanmean(WTPSSBehavior.PSShist.hist,3),...
    nanstd(WTPSSBehavior.PSShist.hist,0,3)./sqrt(size(WTPSSBehavior.PSShist.hist,3)),'k.');
bar(KOPSSBehavior.PSShist.bins(:,:,1),nanmean(KOPSSBehavior.PSShist.hist,3),...
        'facecolor','r','facealpha',0.65);
errorbar(KOPSSBehavior.PSShist.bins(:,:,1),...
    nanmean(KOPSSBehavior.PSShist.hist,3),...
    nanstd(KOPSSBehavior.PSShist.hist,0,3)./sqrt(size(KOPSSBehavior.PSShist.hist,3)),'r.');
xlabel('PSS (au)'); ylabel('counts (norm.)');
axis tight
xlim([-3.5 0]); %ylim([0 1])

cmin = min([min(min(nanmean(WTPSSBehavior.pupcyclePSS.meanZ,3)))...
    min(min(nanmean(KOPSSBehavior.pupcyclePSS.meanZ,3)))]);
cmax = max([max(max(nanmean(WTPSSBehavior.pupcyclePSS.meanZ,3)))...
    max(max(nanmean(KOPSSBehavior.pupcyclePSS.meanZ,3)))]);

subplot(2,4,2)    
h = imagesc(WTPSSBehavior.pupcyclePSS.Xbins(:,:,1),WTPSSBehavior.pupcyclePSS.Ybins(:,:,1),...
    nanmean(WTPSSBehavior.pupcyclePSS.meanZ,3)');
colormap(gca,'jet')
set(h,'AlphaData',~isnan(nanmean(WTPSSBehavior.pupcyclePSS.meanZ,3)'));
hold on
plot(WTPSSBehavior.pupcyclePSS.Xbins([1 end],:,1),[0 0],'k--')
%ColorbarWithAxis([min(min(pupcyclePSS.meanZ)) max(max(pupcyclePSS.meanZ))],['PSS (au)'])
caxis([cmin cmax]);
LogScale('x',10)
axis xy
xlabel('Pupil Area (med^-^1)')
ylabel('dp/dt')

subplot(2,4,6)    
h = imagesc(KOPSSBehavior.pupcyclePSS.Xbins(:,:,1),KOPSSBehavior.pupcyclePSS.Ybins(:,:,1),...
    nanmean(KOPSSBehavior.pupcyclePSS.meanZ,3)');
colormap(gca,'jet')
set(h,'AlphaData',~isnan(nanmean(KOPSSBehavior.pupcyclePSS.meanZ,3)'));
hold on
plot(KOPSSBehavior.pupcyclePSS.Xbins([1 end],:,1),[0 0],'k--')
%ColorbarWithAxis([min(min(pupcyclePSS.meanZ)) max(max(pupcyclePSS.meanZ))],['PSS (au)'])
caxis([cmin cmax]);
LogScale('x',10)
axis xy
xlabel('Pupil Area (med^-^1)')
ylabel('dp/dt')

subplot(2,4,3)
imagesc(WTPSSBehavior.PSSPupdynhist.bins{:,1,1},WTPSSBehavior.PSSPupdynhist.bins{:,2,1},...
    log10(nanmean(WTPSSBehavior.PSSPupdynhist.Pupcounts,3))'); hold on;
colormap(gca,'jet')
%ColorbarWithAxis([min(min(log10(PSSPupdynhist.Pupcounts))) max(max(log10(PSSPupdynhist.Pupcounts)))/2],['counts (au)'])
LogScale('c',10)
xlabel('PSS (au)'); ylabel('Pupil Area (med^-^1)')
%axis square
axis xy

subplot(2,4,7)
imagesc(KOPSSBehavior.PSSPupdynhist.bins{:,1,1},KOPSSBehavior.PSSPupdynhist.bins{:,2,1},...
    log10(nanmean(KOPSSBehavior.PSSPupdynhist.Pupcounts,3))'); hold on;
colormap(gca,'jet')
%ColorbarWithAxis([min(min(log10(PSSPupdynhist.Pupcounts))) max(max(log10(PSSPupdynhist.Pupcounts)))/2],['counts (au)'])
LogScale('c',10)
xlabel('PSS (au)'); ylabel('Pupil Area (med^-^1)')
%axis square
axis xy

subplot(2,4,4)
imagesc(WTPSSBehavior.PSSPupdynhist.bins{:,1,1},WTPSSBehavior.PSSPupdynhist.bins{:,2,1},...
    log10(nanmean(WTPSSBehavior.PSSPupdynhist.dPcounts,3))'); hold on;
colormap(gca,'jet')
%ColorbarWithAxis([min(min(log10(PSSPupdynhist.dPcounts))) max(max(log10(PSSPupdynhist.dPcounts)))/2],['counts (au)'])
LogScale('c',10)
xlabel('PSS (au)'); ylabel('dPdt (med^-^1s^-^1)')
ylim([-0.4 0.4])
axis xy

subplot(2,4,8)
imagesc(KOPSSBehavior.PSSPupdynhist.bins{:,1,1},KOPSSBehavior.PSSPupdynhist.bins{:,2,1},...
    log10(nanmean(KOPSSBehavior.PSSPupdynhist.dPcounts,3))'); hold on;
colormap(gca,'jet')
%ColorbarWithAxis([min(min(log10(PSSPupdynhist.dPcounts))) max(max(log10(PSSPupdynhist.dPcounts)))/2],['counts (au)'])
LogScale('c',10)
xlabel('PSS (au)'); ylabel('dPdt (med^-^1s^-^1)')
ylim([-0.4 0.4])
axis xy

NiceSave('PSS_Phase_dPdt',figfolder,baseName)

%% FIGURE:
cosx = linspace(-pi,3*pi,100);

cmin = min([min(min(nanmean(WTPSSBehavior.pupildist.joint,3)))...
    min(min(nanmean(KOPSSBehavior.pupildist.joint,3)))]);
cmax = max([max(max(nanmean(WTPSSBehavior.pupildist.joint,3)))...
    max(max(nanmean(KOPSSBehavior.pupildist.joint,3)))]);

figure;
subplot(2,2,1);
imagesc(WTPSSBehavior.pupildist.bins{:,1,1},WTPSSBehavior.pupildist.bins{:,2,1},...
    nanmean(WTPSSBehavior.pupildist.joint,3)')
hold on
imagesc(WTPSSBehavior.pupildist.bins{:,1,1}+2*pi,WTPSSBehavior.pupildist.bins{:,2,1},...
    nanmean(WTPSSBehavior.pupildist.joint,3)')
plot(cosx,cos(cosx)./3-1.1,'w','linewidth',2)
xlim([-pi 3*pi])
caxis([cmin cmax]);
axis xy
xlabel('Pupil Phase');ylabel('Pupil Power')
title('P(Pupil power,phase)')

subplot(2,2,3);
imagesc(KOPSSBehavior.pupildist.bins{:,1,1},KOPSSBehavior.pupildist.bins{:,2,1},...
    nanmean(KOPSSBehavior.pupildist.joint,3)')
hold on
imagesc(KOPSSBehavior.pupildist.bins{:,1,1}+2*pi,KOPSSBehavior.pupildist.bins{:,2,1},...
    nanmean(KOPSSBehavior.pupildist.joint,3)')
plot(cosx,cos(cosx)./3-1.1,'w','linewidth',2)
xlim([-pi 3*pi])
caxis([cmin cmax]);
axis xy
xlabel('Pupil Phase');ylabel('Pupil Power')
title('P(Pupil power,phase)')

cmin = min([min(min(nanmean(WTPSSBehavior.pupildist.conditional,3)))...
    min(min(nanmean(KOPSSBehavior.pupildist.conditional,3)))]);
cmax = max([max(max(nanmean(WTPSSBehavior.pupildist.conditional,3)))...
    max(max(nanmean(KOPSSBehavior.pupildist.conditional,3)))]);

subplot(2,2,2);
imagesc(WTPSSBehavior.pupildist.bins{:,1,1},WTPSSBehavior.pupildist.bins{:,2,1},...
    nanmean(WTPSSBehavior.pupildist.conditional,3)')
hold on
imagesc(WTPSSBehavior.pupildist.bins{:,1,1}+2*pi,WTPSSBehavior.pupildist.bins{:,2,1},...
    nanmean(WTPSSBehavior.pupildist.conditional,3)')
plot(cosx,cos(cosx)./3-1.1,'w','linewidth',2)
xlim([-pi 3*pi])
caxis([cmin cmax]);
axis xy
xlabel('Pupil Phase');ylabel('Pupil Power')
title('P(Pupil power|phase)')

subplot(2,2,4);
imagesc(KOPSSBehavior.pupildist.bins{:,1,1},KOPSSBehavior.pupildist.bins{:,2,1},...
    nanmean(KOPSSBehavior.pupildist.conditional,3)')
hold on
imagesc(KOPSSBehavior.pupildist.bins{:,1,1}+2*pi,KOPSSBehavior.pupildist.bins{:,2,1},...
    nanmean(KOPSSBehavior.pupildist.conditional,3)')
plot(cosx,cos(cosx)./3-1.1,'w','linewidth',2)
xlim([-pi 3*pi])
caxis([cmin cmax]);
axis xy
xlabel('Pupil Phase');ylabel('Pupil Power')
title('P(Pupil power|phase)')

NiceSave('PSS_PhasePow',figfolder,baseName)

%% FIGURE:
cmin = min([min(min(nanmean(WTPSSBehavior.pupilPSSdist.joint,3)))...
    min(min(nanmean(KOPSSBehavior.pupilPSSdist.joint,3)))]);
cmax = max([max(max(nanmean(WTPSSBehavior.pupilPSSdist.joint,3)))...
    max(max(nanmean(KOPSSBehavior.pupilPSSdist.joint,3)))]);

figure;
subplot(2,4,1);
imagesc(WTPSSBehavior.pupilPSSdist.bins{:,1,1},WTPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(WTPSSBehavior.pupilPSSdist.joint,3)')
hold on
imagesc(WTPSSBehavior.pupilPSSdist.bins{:,1,1}+2*pi,WTPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(WTPSSBehavior.pupilPSSdist.joint,3)')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil Phase');ylabel('PSS')
title('P(PSS,phase)')

subplot(2,4,5);
imagesc(KOPSSBehavior.pupilPSSdist.bins{:,1,1},KOPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(KOPSSBehavior.pupilPSSdist.joint,3)')
hold on
imagesc(KOPSSBehavior.pupilPSSdist.bins{:,1,1}+2*pi,KOPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(KOPSSBehavior.pupilPSSdist.joint,3)')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil Phase');ylabel('PSS')
title('P(PSS,phase)')

cmin = min([min(min(nanmean(WTPSSBehavior.pupilPSSdist.conditional,3)))...
    min(min(nanmean(KOPSSBehavior.pupilPSSdist.conditional,3)))]);
cmax = max([max(max(nanmean(WTPSSBehavior.pupilPSSdist.conditional,3)))...
    max(max(nanmean(KOPSSBehavior.pupilPSSdist.conditional,3)))]);

subplot(2,4,2);
imagesc(WTPSSBehavior.pupilPSSdist.bins{:,1,1},WTPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(WTPSSBehavior.pupilPSSdist.conditional,3)')
hold on
imagesc(WTPSSBehavior.pupilPSSdist.bins{:,1,1}+2*pi,WTPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(WTPSSBehavior.pupilPSSdist.conditional,3)')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('PSS')
title('P(PSS|phase)')

subplot(2,4,6);
imagesc(KOPSSBehavior.pupilPSSdist.bins{:,1,1},KOPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(KOPSSBehavior.pupilPSSdist.conditional,3)')
hold on
imagesc(KOPSSBehavior.pupilPSSdist.bins{:,1,1}+2*pi,KOPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(KOPSSBehavior.pupilPSSdist.conditional,3)')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('PSS')
title('P(PSS|phase)')

cmin = min([min(min(nanmean(WTPSSBehavior.pupilPSSdist.conditional_high,3)))...
    min(min(nanmean(KOPSSBehavior.pupilPSSdist.conditional_high,3)))]);
cmax = max([max(max(nanmean(WTPSSBehavior.pupilPSSdist.conditional_high,3)))...
    max(max(nanmean(KOPSSBehavior.pupilPSSdist.conditional_high,3)))]);

subplot(2,4,3);
imagesc(WTPSSBehavior.pupilPSSdist.bins{:,1,1},WTPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(WTPSSBehavior.pupilPSSdist.conditional_high,3)')
hold on
imagesc(WTPSSBehavior.pupilPSSdist.bins{:,1,1}+2*pi,WTPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(WTPSSBehavior.pupilPSSdist.conditional_high,3)')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('PSS')
title('>median pupil');

subplot(2,4,7);
imagesc(KOPSSBehavior.pupilPSSdist.bins{:,1,1},KOPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(KOPSSBehavior.pupilPSSdist.conditional_high,3)')
hold on
imagesc(KOPSSBehavior.pupilPSSdist.bins{:,1,1}+2*pi,KOPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(KOPSSBehavior.pupilPSSdist.conditional_high,3)')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('PSS')
title('>median pupil');

cmin = min([min(min(nanmean(WTPSSBehavior.pupilPSSdist.conditional_low,3)))...
    min(min(nanmean(KOPSSBehavior.pupilPSSdist.conditional_low,3)))]);
cmax = max([max(max(nanmean(WTPSSBehavior.pupilPSSdist.conditional_low,3)))...
    max(max(nanmean(KOPSSBehavior.pupilPSSdist.conditional_low,3)))]);

subplot(2,4,4);
imagesc(WTPSSBehavior.pupilPSSdist.bins{:,1,1},WTPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(WTPSSBehavior.pupilPSSdist.conditional_low,3)')
hold on
imagesc(WTPSSBehavior.pupilPSSdist.bins{:,1,1}+2*pi,WTPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(WTPSSBehavior.pupilPSSdist.conditional_low,3)')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil Phase');ylabel('PSS')
title('<median pupil');

subplot(2,4,8);
imagesc(KOPSSBehavior.pupilPSSdist.bins{:,1,1},KOPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(KOPSSBehavior.pupilPSSdist.conditional_low,3)')
hold on
imagesc(KOPSSBehavior.pupilPSSdist.bins{:,1,1}+2*pi,KOPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(KOPSSBehavior.pupilPSSdist.conditional_low,3)')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil Phase');ylabel('PSS')
title('<median pupil');

NiceSave('PSSbyPupilPhase',figfolder,baseName)

%% FIGURE:
cmin = min([min(min(nanmean(WTPSSBehavior.pupilPSSdist.conditional_highw,3)))...
    min(min(nanmean(KOPSSBehavior.pupilPSSdist.conditional_highw,3)))]);
cmax = max([max(max(nanmean(WTPSSBehavior.pupilPSSdist.conditional_highw,3)))...
    max(max(nanmean(KOPSSBehavior.pupilPSSdist.conditional_highw,3)))]);

figure;
subplot(2,4,1);
imagesc(WTPSSBehavior.pupilPSSdist.bins{:,1,1},WTPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(WTPSSBehavior.pupilPSSdist.conditional_highw,3)')
hold on
imagesc(WTPSSBehavior.pupilPSSdist.bins{:,1,1}+2*pi,WTPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(WTPSSBehavior.pupilPSSdist.conditional_highw,3)')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('PSS')
title('>median pupil, wh');

subplot(2,4,5);
imagesc(KOPSSBehavior.pupilPSSdist.bins{:,1,1},KOPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(KOPSSBehavior.pupilPSSdist.conditional_highw,3)')
hold on
imagesc(KOPSSBehavior.pupilPSSdist.bins{:,1,1}+2*pi,KOPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(KOPSSBehavior.pupilPSSdist.conditional_highw,3)')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('PSS')
title('>median pupil, wh');

cmin = min([min(min(nanmean(WTPSSBehavior.pupilPSSdist.conditional_loww,3)))...
    min(min(nanmean(KOPSSBehavior.pupilPSSdist.conditional_loww,3)))]);
cmax = max([max(max(nanmean(WTPSSBehavior.pupilPSSdist.conditional_loww,3)))...
    max(max(nanmean(KOPSSBehavior.pupilPSSdist.conditional_loww,3)))]);

subplot(2,4,2);
imagesc(WTPSSBehavior.pupilPSSdist.bins{:,1,1},WTPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(WTPSSBehavior.pupilPSSdist.conditional_loww,3)')
hold on
imagesc(WTPSSBehavior.pupilPSSdist.bins{:,1,1}+2*pi,WTPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(WTPSSBehavior.pupilPSSdist.conditional_loww,3)')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('PSS')
title('<median pupil, wh');

subplot(2,4,6);
imagesc(KOPSSBehavior.pupilPSSdist.bins{:,1,1},KOPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(KOPSSBehavior.pupilPSSdist.conditional_loww,3)')
hold on
imagesc(KOPSSBehavior.pupilPSSdist.bins{:,1,1}+2*pi,KOPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(KOPSSBehavior.pupilPSSdist.conditional_loww,3)')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('PSS')
title('<median pupil, wh');

cmin = min([min(min(nanmean(WTPSSBehavior.pupilPSSdist.conditional_highnow,3)))...
    min(min(nanmean(KOPSSBehavior.pupilPSSdist.conditional_highnow,3)))]);
cmax = max([max(max(nanmean(WTPSSBehavior.pupilPSSdist.conditional_highnow,3)))...
    max(max(nanmean(KOPSSBehavior.pupilPSSdist.conditional_highnow,3)))]);

subplot(2,4,3);
imagesc(WTPSSBehavior.pupilPSSdist.bins{:,1,1},WTPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(WTPSSBehavior.pupilPSSdist.conditional_highnow,3)')
hold on
imagesc(WTPSSBehavior.pupilPSSdist.bins{:,1,1}+2*pi,WTPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(WTPSSBehavior.pupilPSSdist.conditional_highnow,3)')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('PSS')
title('>median pupil, no whisky');

subplot(2,4,7);
imagesc(KOPSSBehavior.pupilPSSdist.bins{:,1,1},KOPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(KOPSSBehavior.pupilPSSdist.conditional_highnow,3)')
hold on
imagesc(KOPSSBehavior.pupilPSSdist.bins{:,1,1}+2*pi,KOPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(KOPSSBehavior.pupilPSSdist.conditional_highnow,3)')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('PSS')
title('>median pupil, no whisky');

cmin = min([min(min(nanmean(WTPSSBehavior.pupilPSSdist.conditional_lownow,3)))...
    min(min(nanmean(KOPSSBehavior.pupilPSSdist.conditional_lownow,3)))]);
cmax = max([max(max(nanmean(WTPSSBehavior.pupilPSSdist.conditional_lownow,3)))...
    max(max(nanmean(KOPSSBehavior.pupilPSSdist.conditional_lownow,3)))]);

subplot(2,4,4);
imagesc(WTPSSBehavior.pupilPSSdist.bins{:,1,1},WTPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(WTPSSBehavior.pupilPSSdist.conditional_lownow,3)')
hold on
imagesc(WTPSSBehavior.pupilPSSdist.bins{:,1,1}+2*pi,WTPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(WTPSSBehavior.pupilPSSdist.conditional_lownow,3)')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('PSS')
title('<median pupil, no whisky');

subplot(2,4,8);
imagesc(KOPSSBehavior.pupilPSSdist.bins{:,1,1},KOPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(KOPSSBehavior.pupilPSSdist.conditional_lownow,3)')
hold on
imagesc(KOPSSBehavior.pupilPSSdist.bins{:,1,1}+2*pi,KOPSSBehavior.pupilPSSdist.bins{:,2,1},...
    nanmean(KOPSSBehavior.pupilPSSdist.conditional_lownow,3)')
plot(cosx,cos(cosx)./3-3.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('PSS')
title('<median pupil, no whisky');

NiceSave('PSSbyPupilPhase_Whisksep',figfolder,baseName)

%%
load(fullfile(basePath,['WT_EM1M3.BehaviorAnalysis2.mat']));
WTthresh = nanmean(groupBehaviorInts.whthresh);

load(fullfile(basePath,['KO_EM1M3.BehaviorAnalysis2.mat']));
KOthresh = nanmean(groupBehaviorInts.whthresh);

%% FIGURE:
cosx = linspace(-pi,3*pi,100);

cmin = min([min(min(nanmean(WTPSSBehavior.pupilEMGdist.joint,3)))...
    min(min(nanmean(KOPSSBehavior.pupilEMGdist.joint,3)))]);
cmax = max([max(max(nanmean(WTPSSBehavior.pupilEMGdist.joint,3)))...
    max(max(nanmean(KOPSSBehavior.pupilEMGdist.joint,3)))]);

figure;
subplot(2,2,1);
imagesc(WTPSSBehavior.pupilEMGdist.bins{:,1,1},WTPSSBehavior.pupilEMGdist.bins{:,2,1},...
    nanmean(WTPSSBehavior.pupilEMGdist.joint,3)')
hold on
imagesc(WTPSSBehavior.pupilEMGdist.bins{:,1,1}+2*pi,WTPSSBehavior.pupilEMGdist.bins{:,2,1},...
    nanmean(WTPSSBehavior.pupilEMGdist.joint,3)')
plot(cosx,cos(cosx)./3-1.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('EMG')
title('P(EMG,phase)')

subplot(2,2,3);
imagesc(KOPSSBehavior.pupilEMGdist.bins{:,1,1},KOPSSBehavior.pupilEMGdist.bins{:,2,1},...
    nanmean(KOPSSBehavior.pupilEMGdist.joint,3)')
hold on
imagesc(KOPSSBehavior.pupilEMGdist.bins{:,1,1}+2*pi,KOPSSBehavior.pupilEMGdist.bins{:,2,1},...
    nanmean(KOPSSBehavior.pupilEMGdist.joint,3)')
plot(cosx,cos(cosx)./3-1.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('EMG')
title('P(EMG,phase)')

cmin = min([min(min(nanmean(WTPSSBehavior.pupilEMGdist.conditional,3)))...
    min(min(nanmean(KOPSSBehavior.pupilEMGdist.conditional,3)))]);
cmax = max([max(max(nanmean(WTPSSBehavior.pupilEMGdist.conditional,3)))...
    max(max(nanmean(KOPSSBehavior.pupilEMGdist.conditional,3)))]);

subplot(2,2,2);
imagesc(WTPSSBehavior.pupilEMGdist.bins{:,1,1},WTPSSBehavior.pupilEMGdist.bins{:,2,1},...
    nanmean(WTPSSBehavior.pupilEMGdist.conditional,3)')
hold on
imagesc(WTPSSBehavior.pupilEMGdist.bins{:,1,1}+2*pi,WTPSSBehavior.pupilEMGdist.bins{:,2,1},...
    nanmean(WTPSSBehavior.pupilEMGdist.conditional,3)')
plot(cosx,cos(cosx)./3-1.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('EMG')
title('P(EMG|phase)')

subplot(2,2,4);
imagesc(KOPSSBehavior.pupilEMGdist.bins{:,1,1},KOPSSBehavior.pupilEMGdist.bins{:,2,1},...
    nanmean(KOPSSBehavior.pupilEMGdist.conditional,3)')
hold on
imagesc(KOPSSBehavior.pupilEMGdist.bins{:,1,1}+2*pi,KOPSSBehavior.pupilEMGdist.bins{:,2,1},...
    nanmean(KOPSSBehavior.pupilEMGdist.conditional,3)')
plot(cosx,cos(cosx)./3-1.6,'w','linewidth',2)
xlim([-pi 3*pi])
axis xy
caxis([cmin cmax]);
xlabel('Pupil phase');ylabel('EMG')
title('P(EMG|phase)')

NiceSave('EMGbyPupilPhase',figfolder,baseName)

%% FIGURE:
cmin = min([min(min(nanmean(WTPSSBehavior.PSSEMGdist.joint,3)))...
    min(min(nanmean(KOPSSBehavior.PSSEMGdist.joint,3)))]);
cmax = max([max(max(nanmean(WTPSSBehavior.PSSEMGdist.joint,3)))...
    max(max(nanmean(KOPSSBehavior.PSSEMGdist.joint,3)))]);

figure;
subplot(2,2,1);
imagesc(WTPSSBehavior.PSSEMGdist.edges{:,1,1},WTPSSBehavior.PSSEMGdist.edges{:,2,1},...
    nanmean(WTPSSBehavior.PSSEMGdist.joint,3)')
hold on
axis square
caxis([cmin cmax]);
xlabel('EMG');ylabel('PSS')
title('P(EMG,PSS)')

subplot(2,2,3);
imagesc(KOPSSBehavior.PSSEMGdist.edges{:,1,1},KOPSSBehavior.PSSEMGdist.edges{:,2,1},...
    nanmean(KOPSSBehavior.PSSEMGdist.joint,3)')
hold on
axis square
caxis([cmin cmax]);
xlabel('EMG');ylabel('PSS')
title('P(EMG,PSS)')

cmin = min([min(min(nanmean(WTPSSBehavior.PSSEMGdist.conditional,3)))...
    min(min(nanmean(KOPSSBehavior.PSSEMGdist.conditional,3)))]);
cmax = max([max(max(nanmean(WTPSSBehavior.PSSEMGdist.conditional,3)))...
    max(max(nanmean(KOPSSBehavior.PSSEMGdist.conditional,3)))]);

subplot(2,2,2);
imagesc(WTPSSBehavior.PSSEMGdist.edges{:,1,1},WTPSSBehavior.PSSEMGdist.edges{:,2,1},...
    nanmean(WTPSSBehavior.PSSEMGdist.conditional,3)')
hold on
axis square
caxis([cmin cmax]);
xlabel('EMG');ylabel('PSS')
title('P(EMG|PSS)')

subplot(2,2,4);
imagesc(KOPSSBehavior.PSSEMGdist.edges{:,1,1},KOPSSBehavior.PSSEMGdist.edges{:,2,1},...
    nanmean(KOPSSBehavior.PSSEMGdist.conditional,3)')
hold on
axis square
caxis([cmin cmax]);
xlabel('EMG');ylabel('PSS')
title('P(EMG|PSS)')

NiceSave('EMGPSS',figfolder,baseName)

%% FIGURE:
cmin = min([min(min(nanmean(WTPSSBehavior.EMGPupPSS.meanZ,3)))...
    min(min(nanmean(KOPSSBehavior.EMGPupPSS.meanZ,3)))]);
cmax = max([max(max(nanmean(WTPSSBehavior.EMGPupPSS.meanZ,3)))...
    max(max(nanmean(KOPSSBehavior.EMGPupPSS.meanZ,3)))]);

figure;
subplot(2,3,1);
a = imagesc([WTPSSBehavior.EMGPupPSS.Xbins(:,:,1) WTPSSBehavior.EMGPupPSS.Xbins(:,:,1)+2*pi],...
    WTPSSBehavior.EMGPupPSS.Ybins(:,:,1),...
    [nanmean(WTPSSBehavior.EMGPupPSS.meanZ,3); nanmean(WTPSSBehavior.EMGPupPSS.meanZ,3)]');
colormap(gca,'jet')
alpha(a,double(~isnan([nanmean(WTPSSBehavior.EMGPupPSS.meanZ,3); nanmean(WTPSSBehavior.EMGPupPSS.meanZ,3)]')))
caxis([cmin cmax]);
ylim([-3 0]);
xlabel('Pupil phase');ylabel('PSS')
axis xy
title('EMG by PSS/Pupil');

subplot(2,3,4);
a = imagesc([KOPSSBehavior.EMGPupPSS.Xbins(:,:,1) KOPSSBehavior.EMGPupPSS.Xbins(:,:,1)+2*pi],...
    KOPSSBehavior.EMGPupPSS.Ybins(:,:,1),...
    [nanmean(KOPSSBehavior.EMGPupPSS.meanZ,3); nanmean(KOPSSBehavior.EMGPupPSS.meanZ,3)]');
colormap(gca,'jet')
alpha(a,double(~isnan([nanmean(KOPSSBehavior.EMGPupPSS.meanZ,3); nanmean(KOPSSBehavior.EMGPupPSS.meanZ,3)]')))
caxis([cmin cmax]);
ylim([-3 0]);
xlabel('Pupil phase');ylabel('PSS')
axis xy
title('EMG by PSS/Pupil');

cmin = min([min(min(nanmean(WTPSSBehavior.PSSbyhiPup.meanZ,3)))...
    min(min(nanmean(KOPSSBehavior.PSSbyhiPup.meanZ,3)))...
    min(min(nanmean(WTPSSBehavior.PSSbyloPup.meanZ,3)))...
    min(min(nanmean(KOPSSBehavior.PSSbyloPup.meanZ,3)))]);
cmax = max([max(max(nanmean(WTPSSBehavior.PSSbyhiPup.meanZ,3)))...
    max(max(nanmean(KOPSSBehavior.PSSbyhiPup.meanZ,3)))...
    max(max(nanmean(WTPSSBehavior.PSSbyloPup.meanZ,3)))...
    max(max(nanmean(KOPSSBehavior.PSSbyloPup.meanZ,3)))]);

subplot(2,3,2);
a = imagesc([WTPSSBehavior.PSSbyhiPup.Xbins(:,:,1) WTPSSBehavior.PSSbyhiPup.Xbins(:,:,1)+2*pi],...
    WTPSSBehavior.PSSbyhiPup.Ybins(:,:,1),...
    [nanmean(WTPSSBehavior.PSSbyhiPup.meanZ,3); nanmean(WTPSSBehavior.PSSbyhiPup.meanZ,3)]');
colormap(gca,'jet')
alpha(a,double(~isnan([nanmean(WTPSSBehavior.PSSbyhiPup.meanZ,3); nanmean(WTPSSBehavior.PSSbyhiPup.meanZ,3)]')))
% ColorbarWithAxis([min(min(PSSbyhiPup.meanZ)) max(max(PSSbyhiPup.meanZ))],'PSS (au)')
caxis([cmin cmax]);
ylim([-2 1.5]);
xlabel('Pupil phase');ylabel('EMG')
axis xy
title('>median pupil');

subplot(2,3,5);
a = imagesc([KOPSSBehavior.PSSbyhiPup.Xbins(:,:,1) KOPSSBehavior.PSSbyhiPup.Xbins(:,:,1)+2*pi],...
    KOPSSBehavior.PSSbyhiPup.Ybins(:,:,1),...
    [nanmean(KOPSSBehavior.PSSbyhiPup.meanZ,3); nanmean(KOPSSBehavior.PSSbyhiPup.meanZ,3)]');
colormap(gca,'jet')
alpha(a,double(~isnan([nanmean(KOPSSBehavior.PSSbyhiPup.meanZ,3); nanmean(KOPSSBehavior.PSSbyhiPup.meanZ,3)]')))
% ColorbarWithAxis([min(min(PSSbyhiPup.meanZ)) max(max(PSSbyhiPup.meanZ))],'PSS (au)')
caxis([cmin cmax]);
ylim([-2 1.5]);
xlabel('Pupil phase');ylabel('EMG')
axis xy
title('>median pupil');

subplot(2,3,3);
a = imagesc([WTPSSBehavior.PSSbyloPup.Xbins(:,:,1) WTPSSBehavior.PSSbyloPup.Xbins(:,:,1)+2*pi],...
    WTPSSBehavior.PSSbyloPup.Ybins(:,:,1),...
    [nanmean(WTPSSBehavior.PSSbyloPup.meanZ,3); nanmean(WTPSSBehavior.PSSbyloPup.meanZ,3)]');
colormap(gca,'jet')
alpha(a,double(~isnan([nanmean(WTPSSBehavior.PSSbyloPup.meanZ,3); nanmean(WTPSSBehavior.PSSbyloPup.meanZ,3)]')))
% ColorbarWithAxis([min(min(PSSbyloPup.meanZ)) max(max(PSSbyloPup.meanZ))],'PSS (au)')
caxis([cmin cmax]);
ylim([-2 1.5]);
xlabel('Pupil phase');ylabel('EMG')
axis xy
title('<median pupil');

subplot(2,3,6);
a = imagesc([KOPSSBehavior.PSSbyloPup.Xbins(:,:,1) KOPSSBehavior.PSSbyloPup.Xbins(:,:,1)+2*pi],...
    KOPSSBehavior.PSSbyloPup.Ybins(:,:,1),...
    [nanmean(KOPSSBehavior.PSSbyloPup.meanZ,3); nanmean(KOPSSBehavior.PSSbyloPup.meanZ,3)]');
colormap(gca,'jet')
alpha(a,double(~isnan([nanmean(KOPSSBehavior.PSSbyloPup.meanZ,3); nanmean(KOPSSBehavior.PSSbyloPup.meanZ,3)]')))
% ColorbarWithAxis([min(min(PSSbyloPup.meanZ)) max(max(PSSbyloPup.meanZ))],'PSS (au)')
caxis([cmin cmax]);
ylim([-2 1.5]);
xlabel('Pupil phase');ylabel('EMG')
axis xy
title('<median pupil');

NiceSave('PupilEMGPSS',figfolder,baseName)

%% FIGURE:
cmin = min([min(min(nanmean(WTPSSBehavior.phasePETH.high.mean,3)))...
    min(min(nanmean(WTPSSBehavior.phasePETH.low.mean,3)))...
    min(min(nanmean(KOPSSBehavior.phasePETH.high.mean,3)))...
    min(min(nanmean(KOPSSBehavior.phasePETH.low.mean,3)))]);
cmax = max([max(max(nanmean(WTPSSBehavior.phasePETH.high.mean,3)))...
    max(max(nanmean(WTPSSBehavior.phasePETH.low.mean,3)))...
    max(max(nanmean(KOPSSBehavior.phasePETH.high.mean,3)))...
    max(max(nanmean(KOPSSBehavior.phasePETH.low.mean,3)))]);

figure;
subplot(2,2,1);
imagesc(WTPSSBehavior.phasePETH.high.bincenters(:,:,1),WTPSSBehavior.phasePETH.high.bincenters(:,:,1),...
    nanmean(WTPSSBehavior.phasePETH.high.mean,3)'); hold on
imagesc(WTPSSBehavior.phasePETH.high.bincenters(:,:,1),WTPSSBehavior.phasePETH.high.bincenters(:,:,1)+2*pi,...
    nanmean(WTPSSBehavior.phasePETH.high.mean,3)')
plot(cos(cosx),cosx,'w','linewidth',2)
plot([0 0],[-pi 3*pi],'r')
colorbar
axis xy
caxis([cmin cmax]);
xlim([-1 4]);ylim([-pi 3*pi])
xlabel('t (s, aligned to Wh Onset)');ylabel('Pupil Phase')
title('>median pupil')

subplot(2,2,3);
imagesc(KOPSSBehavior.phasePETH.high.bincenters(:,:,1),KOPSSBehavior.phasePETH.high.bincenters(:,:,1),...
    nanmean(KOPSSBehavior.phasePETH.high.mean,3)'); hold on
imagesc(KOPSSBehavior.phasePETH.high.bincenters(:,:,1),KOPSSBehavior.phasePETH.high.bincenters(:,:,1)+2*pi,...
    nanmean(KOPSSBehavior.phasePETH.high.mean,3)')
plot(cos(cosx),cosx,'w','linewidth',2)
plot([0 0],[-pi 3*pi],'r')
colorbar
axis xy
caxis([cmin cmax]);
xlim([-1 4]);ylim([-pi 3*pi])
xlabel('t (s, aligned to Wh Onset)');ylabel('Pupil Phase')
title('>median pupil')

subplot(2,2,2);
imagesc(WTPSSBehavior.phasePETH.low.bincenters(:,:,1),WTPSSBehavior.phasePETH.low.bincenters(:,:,1),...
    nanmean(WTPSSBehavior.phasePETH.low.mean,3)'); hold on
imagesc(WTPSSBehavior.phasePETH.low.bincenters(:,:,1),WTPSSBehavior.phasePETH.low.bincenters(:,:,1)+2*pi,...
    nanmean(WTPSSBehavior.phasePETH.low.mean,3)')
plot(cos(cosx),cosx,'w','linewidth',2)
plot([0 0],[-pi 3*pi],'r')
colorbar
axis xy
caxis([cmin cmax]);
xlim([-1 4]);ylim([-pi 3*pi])
xlabel('t (s, aligned to Wh Onset)');ylabel('Pupil Phase')
title('<median pupil')

subplot(2,2,4);
imagesc(KOPSSBehavior.phasePETH.low.bincenters(:,:,1),KOPSSBehavior.phasePETH.low.bincenters(:,:,1),...
    nanmean(KOPSSBehavior.phasePETH.low.mean,3)'); hold on
imagesc(KOPSSBehavior.phasePETH.low.bincenters(:,:,1),KOPSSBehavior.phasePETH.low.bincenters(:,:,1)+2*pi,...
    nanmean(KOPSSBehavior.phasePETH.low.mean,3)')
plot(cos(cosx),cosx,'w','linewidth',2)
plot([0 0],[-pi 3*pi],'r')
colorbar
axis xy
caxis([cmin cmax]);
xlim([-1 4]);ylim([-pi 3*pi])
xlabel('t (s, aligned to Wh Onset)');ylabel('Pupil Phase')
title('<median pupil')

NiceSave('zPSS_PETHbyPhase',figfolder,baseName)

%%
[~,WTwhisksorts.phaseidx] = sort(WTwhisksorts.phase);
[~,WTwhisksorts.duridx] = sort(WTwhisksorts.dur);
[~,WTwhisksorts.ampidx] = sort(WTwhisksorts.amp);
[~,WTwhisksorts.pssidx] = sort(WTwhisksorts.maxpss);

[~,KOwhisksorts.phaseidx] = sort(KOwhisksorts.phase);
[~,KOwhisksorts.duridx] = sort(KOwhisksorts.dur);
[~,KOwhisksorts.ampidx] = sort(KOwhisksorts.amp);
[~,KOwhisksorts.pssidx] = sort(KOwhisksorts.maxpss);

%% FIGURE:
figure;
subplot(2,3,1);
imagesc(WTtimelockedPSS.timestamps(:,1),[1 size(WTtimelockedPSS.data,2)],...
    WTtimelockedPSS.data(:,WTwhisksorts.phaseidx)')
hold on
plot(WTwhisksorts.dur(WTwhisksorts.phaseidx),[1:size(WTtimelockedPSS.data,2)],'.r','markersize',2)
plot([0 0],[0 size(WTtimelockedPSS.data,2)],'b')
axis tight
xlim([-2 5])
colorbar; caxis([-3 -1])
xlabel('t (s, aligned to Wh Onset)');ylabel('trial no.')
title('Epochs sorted by pupil phase');

subplot(2,3,4);
imagesc(KOtimelockedPSS.timestamps(:,1),[1 size(KOtimelockedPSS.data,2)],...
    KOtimelockedPSS.data(:,KOwhisksorts.phaseidx)')
hold on
plot(KOwhisksorts.dur(KOwhisksorts.phaseidx),[1:size(KOtimelockedPSS.data,2)],'.r','markersize',2)
plot([0 0],[0 size(KOtimelockedPSS.data,2)],'b')
axis tight
xlim([-2 5])
colorbar; caxis([-3 -1])
xlabel('t (s, aligned to Wh Onset)');ylabel('trial no.')
title('Epochs sorted by pupil phase');

subplot(2,3,2);
imagesc(WTtimelockedPSS.timestamps(:,1),[1 size(WTtimelockedPSS.data,2)],...
    WTtimelockedPSS.data(:,WTwhisksorts.duridx)'); 
hold on
plot(WTwhisksorts.dur(WTwhisksorts.duridx),[1:size(WTtimelockedPSS.data,2)],'.r','markersize',1)
plot([0 0],[0 size(WTtimelockedPSS.data,2)],'b')
axis tight
xlim([-2 5])
colorbar; caxis([-3 -1])
xlabel('t (s, aligned to Wh Onset)');ylabel('trial no.')
title('Epochs sorted by Wh duration');

subplot(2,3,5);
imagesc(KOtimelockedPSS.timestamps(:,1),[1 size(KOtimelockedPSS.data,2)],...
    KOtimelockedPSS.data(:,KOwhisksorts.duridx)')
hold on
plot(KOwhisksorts.dur(KOwhisksorts.duridx),[1:size(KOtimelockedPSS.data,2)],'.r','markersize',1)
plot([0 0],[0 size(KOtimelockedPSS.data,2)],'b')
axis tight
xlim([-2 5])
colorbar; caxis([-3 -1])
xlabel('t (s, aligned to Wh Onset)');ylabel('trial no.')
title('Epochs sorted by Wh duration');

subplot(2,3,3);
imagesc(WTtimelockedPSS.timestamps(:,1),[1 size(WTtimelockedPSS.data,2)],...
    WTtimelockedPSS.data(:,WTwhisksorts.ampidx)'); 
hold on
plot(WTwhisksorts.dur(WTwhisksorts.ampidx),[1:size(WTtimelockedPSS.data,2)],'.r','markersize',1)
plot([0 0],[0 size(WTtimelockedPSS.data,2)],'b')
axis tight
xlim([-2 5])
colorbar; caxis([-3 -1])
xlabel('t (s, aligned to Wh Onset)');ylabel('trial no.')
title('Epochs sorted by Wh amp');

subplot(2,3,6);
imagesc(KOtimelockedPSS.timestamps(:,1),[1 size(KOtimelockedPSS.data,2)],...
    KOtimelockedPSS.data(:,KOwhisksorts.ampidx)')
hold on
plot(KOwhisksorts.dur(KOwhisksorts.ampidx),[1:size(KOtimelockedPSS.data,2)],'.r','markersize',1)
plot([0 0],[0 size(KOtimelockedPSS.data,2)],'b')
axis tight
xlim([-2 5])
colorbar; caxis([-3 -1])
xlabel('t (s, aligned to Wh Onset)');ylabel('trial no.')
title('Epochs sorted by Wh amp');

NiceSave('PSS_sortedWhisks',figfolder,baseName)

%%
[WTEMGPhase.meanZ,WTEMGPhase.N,WTEMGPhase.Xbins,WTEMGPhase.Ybins ] = ConditionalHist3( log10(WTwhisksorts.dur),...
    log10(WTwhisksorts.amp),WTwhisksorts.phase,...
    'minXY',0,'Xbounds',[-1 1.5],'Ybounds',[-1 1.5],...
    'numXbins',100,'numYbins',100);

[WTEMGPSS.meanZ,WTEMGPSS.N,WTEMGPSS.Xbins,WTEMGPSS.Ybins ] = ConditionalHist3( log10(WTwhisksorts.dur),...
    log10(WTwhisksorts.amp),WTwhisksorts.maxpss,...
    'minXY',0,'Xbounds',[-1 1.5],'Ybounds',[-1 1.5],...
    'numXbins',100,'numYbins',100);

[KOEMGPhase.meanZ,KOEMGPhase.N,KOEMGPhase.Xbins,KOEMGPhase.Ybins ] = ConditionalHist3( log10(KOwhisksorts.dur),...
    log10(KOwhisksorts.amp),KOwhisksorts.phase,...
    'minXY',0,'Xbounds',[-1 1.5],'Ybounds',[-1 1.5],...
    'numXbins',100,'numYbins',100);

[KOEMGPSS.meanZ,KOEMGPSS.N,KOEMGPSS.Xbins,KOEMGPSS.Ybins ] = ConditionalHist3( log10(KOwhisksorts.dur),...
    log10(KOwhisksorts.amp),KOwhisksorts.maxpss,...
    'minXY',0,'Xbounds',[-1 1.5],'Ybounds',[-1 1.5],...
    'numXbins',100,'numYbins',100);

%% FIGURE:
cmin = min([min(min(WTEMGPhase.meanZ))...
    min(min(KOEMGPhase.meanZ))]);
cmax = max([max(max(WTEMGPhase.meanZ))...
    max(max(KOEMGPhase.meanZ))]);

figure;    
subplot(2,2,1);
h = imagesc(WTEMGPhase.Xbins,WTEMGPhase.Ybins,WTEMGPhase.meanZ');
colormap(gca,'jet')
set(h,'AlphaData',~isnan(WTEMGPhase.meanZ'));
LogScale('x',10); LogScale('y',10);
axis xy
caxis([cmin cmax]);
colorbar
xlabel('Wh duration (s)')
ylabel('EMG amplitude (modZ)')
title('Pupil phase')

subplot(2,2,3);
h = imagesc(KOEMGPhase.Xbins,KOEMGPhase.Ybins,KOEMGPhase.meanZ');
colormap(gca,'jet')
set(h,'AlphaData',~isnan(KOEMGPhase.meanZ'));
LogScale('x',10); LogScale('y',10);
axis xy
caxis([cmin cmax]);
colorbar
xlabel('Wh duration (s)')
ylabel('EMG amplitude (modZ)')
title('Pupil phase')

cmin = min([min(min(WTEMGPSS.meanZ))...
    min(min(KOEMGPSS.meanZ))]);
cmax = max([max(max(WTEMGPSS.meanZ))...
    max(max(KOEMGPSS.meanZ))]);

subplot(2,2,2);
h = imagesc(WTEMGPSS.Xbins,WTEMGPSS.Ybins,WTEMGPSS.meanZ');
colormap(gca,'jet')
set(h,'AlphaData',~isnan(WTEMGPSS.meanZ'));
LogScale('x',10); LogScale('y',10);
axis xy
caxis([cmin cmax]);
colorbar
xlabel('Wh duration (s)')
ylabel('EMG amplitude (modZ)')
title('PSS')

subplot(2,2,4);
h = imagesc(KOEMGPSS.Xbins,KOEMGPSS.Ybins,KOEMGPSS.meanZ');
colormap(gca,'jet')
set(h,'AlphaData',~isnan(KOEMGPSS.meanZ'));
LogScale('x',10); LogScale('y',10);
axis xy
caxis([cmin cmax]);
colorbar
xlabel('Wh duration (s)')
ylabel('EMG amplitude (modZ)')
title('PSS')

NiceSave('PSS_EMGSpace',figfolder,baseName)

%%
[~,WTpupsorts.Pupwhdiffidx] = sort(WTpupsorts.Pupwhdiff);
[~,WTpupsorts.dPidx] = sort(WTpupsorts.dP);
[~,WTpupsorts.duridx] = sort(WTpupsorts.dur);
[~,WTpupsorts.peakidx] = sort(WTpupsorts.peak);
[~,WTpupsorts.phaseidx] = sort(WTpupsorts.phase);
[~,WTpupsorts.powidx] = sort(WTpupsorts.pow);

[~,KOpupsorts.Pupwhdiffidx] = sort(KOpupsorts.Pupwhdiff);
[~,KOpupsorts.dPidx] = sort(KOpupsorts.dP);
[~,KOpupsorts.duridx] = sort(KOpupsorts.dur);
[~,KOpupsorts.peakidx] = sort(KOpupsorts.peak);
[~,KOpupsorts.phaseidx] = sort(KOpupsorts.phase);
[~,KOpupsorts.powidx] = sort(KOpupsorts.pow);

%% FIGURE:
figure;
subplot(2,3,1);
imagesc(WTpuplockedPSS.timestamps(:,1),[1 size(WTpuplockedPSS.data,2)],...
    WTpuplockedPSS.data(:,WTpupsorts.phaseidx)')
hold on
%plot(WTwhisksorts.dur(WTwhisksorts.phaseidx),[1:size(WTpuplockedPSS.data,2)],'.r','markersize',2)
plot([0 0],[0 size(WTpuplockedPSS.data,2)],'b')
axis tight
xlim([-10 10])
colorbar; caxis([-3 -1])
xlabel('t (s, aligned to Pup Onset)'); ylabel('trial no.')
title('Sorted by pupil phase');

subplot(2,3,4);
imagesc(KOpuplockedPSS.timestamps(:,1),[1 size(KOpuplockedPSS.data,2)],...
    KOpuplockedPSS.data(:,KOpupsorts.phaseidx)')
hold on
%plot(KOwhisksorts.dur(KOwhisksorts.phaseidx),[1:size(KOpuplockedPSS.data,2)],'.r','markersize',2)
plot([0 0],[0 size(KOpuplockedPSS.data,2)],'b')
axis tight
xlim([-10 10])
colorbar; caxis([-3 -1])
xlabel('t (s, aligned to Pup Onset)'); ylabel('trial no.')
title('Sorted by pupil phase');

subplot(2,3,2);
imagesc(WTpuplockedPSS.timestamps(:,1),[1 size(WTpuplockedPSS.data,2)],...
    WTpuplockedPSS.data(:,WTpupsorts.powidx)')
hold on
%plot(WTwhisksorts.dur(WTwhisksorts.powidx),[1:size(WTpuplockedPSS.data,2)],'.r','markersize',2)
plot([0 0],[0 size(WTpuplockedPSS.data,2)],'b')
axis tight
xlim([-10 10])
colorbar; caxis([-3 -1])
xlabel('t (s, aligned to Pup Onset)'); ylabel('trial no.')
title('Sorted by pupil power');

subplot(2,3,5);
imagesc(KOpuplockedPSS.timestamps(:,1),[1 size(KOpuplockedPSS.data,2)],...
    KOpuplockedPSS.data(:,KOpupsorts.powidx)')
hold on
%plot(KOwhisksorts.dur(KOwhisksorts.powidx),[1:size(KOpuplockedPSS.data,2)],'.r','markersize',2)
plot([0 0],[0 size(KOpuplockedPSS.data,2)],'b')
axis tight
xlim([-10 10])
colorbar; caxis([-3 -1])
xlabel('t (s, aligned to Pup Onset)'); ylabel('trial no.')
title('Sorted by pupil power');

subplot(2,3,3);
imagesc(WTpuplockedPSS.timestamps(:,1),[1 size(WTpuplockedPSS.data,2)],...
    WTpuplockedPSS.data(:,WTpupsorts.Pupwhdiffidx)')
hold on
plot(WTpupsorts.Pupwhdiff(WTpupsorts.Pupwhdiffidx),[1:size(WTpuplockedPSS.data,2)],'.r','markersize',2)
plot([0 0],[0 size(WTpuplockedPSS.data,2)],'b')
axis tight
xlim([-10 10])
colorbar; caxis([-3 -1])
xlabel('t (s, aligned to Pup Onset)'); ylabel('trial no.')
title('Sorted by Wh-Pup interval');

subplot(2,3,6);
imagesc(KOpuplockedPSS.timestamps(:,1),[1 size(KOpuplockedPSS.data,2)],...
    KOpuplockedPSS.data(:,KOpupsorts.Pupwhdiffidx)')
hold on
plot(KOpupsorts.Pupwhdiff(KOpupsorts.Pupwhdiffidx),[1:size(KOpuplockedPSS.data,2)],'.r','markersize',2)
plot([0 0],[0 size(KOpuplockedPSS.data,2)],'b')
axis tight
xlim([-10 10])
colorbar; caxis([-3 -1])
xlabel('t (s, aligned to Pup Onset)'); ylabel('trial no.')
title('Sorted by Wh-Pup interval');

NiceSave('PSS_sortedPupil',figfolder,baseName)
