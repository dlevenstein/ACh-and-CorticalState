%%
basePath = pwd;
baseName = 'EM1M3';

figfolder = fullfile(basePath,'SummaryFigures');
if (~exist(figfolder,'dir'))
    mkdir(figfolder)
end

load(fullfile(basePath,['WT_EM1M3.BehaviorAnalysis2.mat']));
WTBehavior = groupBehavior;
WTBehaviorInts = groupBehaviorInts;

load(fullfile(basePath,['KO_EM1M3.BehaviorAnalysis2.mat']));
KOBehavior = groupBehavior;
KOBehaviorInts = groupBehaviorInts;

%% FIGURE:
figure;
subplot(2,2,1);
bar(WTBehavior.puphist.bins(:,:,1),nanmean(WTBehavior.puphist.counts,3),...
'facecolor','k','facealpha',0.85); hold on;    
errorbar(WTBehavior.puphist.bins(:,:,1),nanmean(WTBehavior.puphist.counts,3),...
    nanstd(WTBehavior.puphist.counts,0,3),'k.');
bar(KOBehavior.puphist.bins(:,:,1),nanmean(KOBehavior.puphist.counts,3),...
    'facecolor','r','facealpha',0.65);
errorbar(KOBehavior.puphist.bins(:,:,1),nanmean(KOBehavior.puphist.counts,3),...
    nanstd(KOBehavior.puphist.counts,0,3),'r.');
axis tight;
xlim([0 3]); ylim([0 0.5])
xlabel('pupil area'); ylabel('counts (norm)');

subplot(2,2,2);
shadederror(log10(WTBehavior.pupPSD.freqs(:,:,1)),...
    nanmean(WTBehavior.pupPSD.psd,3)',...
    nanstd(WTBehavior.pupPSD.psd,0,3)','k','none',1,'k');
shadederror(log10(KOBehavior.pupPSD.freqs(:,:,1)),...
    nanmean(KOBehavior.pupPSD.psd,3)',...
    nanstd(KOBehavior.pupPSD.psd,0,3)','r','none',1,'r');
LogScale('x',10)
xlabel('f (Hz)'); ylabel('power (dB)')
axis tight
ylim([-1.5 1.5]);
title('Pupil diameter power spectrum');

subplot(2,2,3);
imagesc(WTBehavior.pupdthist.bins{:,1,1},WTBehavior.pupdthist.bins{:,2,1},...
    nanmean(WTBehavior.pupdthist.counts,3)'); hold on;
colormap(gca,[1 1 1; colormap('jet')])
plot(get(gca,'xlim'),[0 0],'r-')
% ColorbarWithAxis([-0.1 1],['counts (au)'])
caxis([-0.0005 0.015]);
xlabel('diameter (norm)'); ylabel('dP/dt')
LogScale('x',10);
axis square
ylim([-0.25 0.25]); xlim([-0.25 0.25]);
title('floxed M1M3');

subplot(2,2,4);
imagesc(KOBehavior.pupdthist.bins{:,1,1},KOBehavior.pupdthist.bins{:,2,1},...
    nanmean(KOBehavior.pupdthist.counts,3)'); hold on;
colormap(gca,[1 1 1; colormap('jet')])
plot(get(gca,'xlim'),[0 0],'r-')
% ColorbarWithAxis([-0.1 0.025],['counts (au)'])
caxis([-0.0005 0.015]);
xlabel('diameter (norm)'); ylabel('dP/dt')
LogScale('x',10);
axis square
ylim([-0.25 0.25]); xlim([-0.25 0.25]);
title('Emx1;M1M3');

NiceSave('PupilStats',figfolder,baseName)

%% FIGURE:
figure;
subplot(1,2,1)
a = imagesc(WTBehavior.phasedynamics.Xbins(:,:,1),WTBehavior.phasedynamics.Ybins(:,:,1),...
    circ_mean(WTBehavior.phasedynamics.meanZ,[],3)');
colormap(gca,'hsv')
alpha(a,double(~isnan(circ_mean(WTBehavior.phasedynamics.meanZ,[],3)')))
ColorbarWithAxis([-pi pi],'phase')
axis xy
xlabel('diameter');ylabel('dP/dt')
axis square
title('floxed M1M3');

subplot(1,2,2)
a = imagesc(KOBehavior.phasedynamics.Xbins(:,:,1),KOBehavior.phasedynamics.Ybins(:,:,1),...
    circ_mean(KOBehavior.phasedynamics.meanZ,[],3)');
colormap(gca,'hsv')
alpha(a,double(~isnan(circ_mean(KOBehavior.phasedynamics.meanZ,[],3)')))
ColorbarWithAxis([-pi pi],'phase')
axis xy
xlabel('diameter');ylabel('dP/dt')
axis square
title('Emx1;M1M3');

% subplot(2,2,3);
% a = imagesc(WTBehavior.ampdynamics.Xbins(:,:,1),WTBehavior.ampdynamics.Ybins(:,:,1),...
%     nanmean(WTBehavior.ampdynamics.meanZ,3)');
% colormap(gca,'jet')
% alpha(a,double(~isnan(nanmean(WTBehavior.ampdynamics.meanZ,3)')))
% %ColorbarWithAxis([0 1.2],'power')
% caxis([0 1.2])
% axis xy
% xlabel('diameter');ylabel('dP/dt')
% axis square
% 
% subplot(2,2,4);
% a = imagesc(KOBehavior.ampdynamics.Xbins(:,:,1),KOBehavior.ampdynamics.Ybins(:,:,1),...
%     nanmean(KOBehavior.ampdynamics.meanZ,3)');
% colormap(gca,'jet')
% alpha(a,double(~isnan(nanmean(KOBehavior.ampdynamics.meanZ,3)')))
% %ColorbarWithAxis([0 1.2],'power')
% caxis([0 1.2])
% axis xy
% xlabel('diameter');ylabel('dP/dt')
% axis square

NiceSave('Pupil_PupSpace1',figfolder,baseName)

%% FIGURE:
cmin = min([min(min(nanmean(WTBehavior.areadynamics.meanZ,3)))...
    min(min(nanmean(KOBehavior.areadynamics.meanZ,3)))]);
cmax = max([max(max(nanmean(WTBehavior.areadynamics.meanZ,3)))...
    max(max(nanmean(KOBehavior.areadynamics.meanZ,3)))]);

figure;
subplot(2,2,1);
a = imagesc([WTBehavior.areadynamics.Xbins(:,:,1) WTBehavior.areadynamics.Xbins(:,:,1)+2*pi],...
    WTBehavior.areadynamics.Ybins(:,:,1),[nanmean(WTBehavior.areadynamics.meanZ,3);...
    nanmean(WTBehavior.areadynamics.meanZ,3)]');
colormap(gca,'jet')
alpha(a,double(~isnan([nanmean(WTBehavior.areadynamics.meanZ,3);...
    nanmean(WTBehavior.areadynamics.meanZ,3)]')))
%ColorbarWithAxis([min(min(areadynamics.meanZ)) max(max(areadynamics.meanZ))],'pupil area (au)')
caxis([cmin cmax/2])
ylim([-1.5 0]);
xlabel('Phase');ylabel('Amplitude')
axis xy

subplot(2,2,2);
a = imagesc([KOBehavior.areadynamics.Xbins(:,:,1) KOBehavior.areadynamics.Xbins(:,:,1)+2*pi],...
    KOBehavior.areadynamics.Ybins(:,:,1),[nanmean(KOBehavior.areadynamics.meanZ,3);...
    nanmean(KOBehavior.areadynamics.meanZ,3)]');
colormap(gca,'jet')
alpha(a,double(~isnan([nanmean(KOBehavior.areadynamics.meanZ,3);...
    nanmean(KOBehavior.areadynamics.meanZ,3)]')))
%ColorbarWithAxis([min(min(areadynamics.meanZ)) max(max(areadynamics.meanZ))],'pupil area (au)')
caxis([cmin cmax/2])
ylim([-1.5 0]);
xlabel('Phase');ylabel('Amplitude')
axis xy

cmin = min([min(min(nanmean(WTBehavior.dpdtdynamics.meanZ,3)))...
    min(min(nanmean(KOBehavior.dpdtdynamics.meanZ,3)))]);
cmax = max([max(max(nanmean(WTBehavior.dpdtdynamics.meanZ,3)))...
    max(max(nanmean(KOBehavior.dpdtdynamics.meanZ,3)))]);

subplot(2,2,3);
a = imagesc([WTBehavior.dpdtdynamics.Xbins(:,:,1) WTBehavior.dpdtdynamics.Xbins(:,:,1)+2*pi],...
    WTBehavior.dpdtdynamics.Ybins(:,:,1),[nanmean(WTBehavior.dpdtdynamics.meanZ,3);...
    nanmean(WTBehavior.dpdtdynamics.meanZ,3)]');
colormap(gca,'jet')
alpha(a,double(~isnan([nanmean(WTBehavior.dpdtdynamics.meanZ,3);...
    nanmean(WTBehavior.dpdtdynamics.meanZ,3)]')))
%ColorbarWithAxis([min(min(dpdtdynamics.meanZ)) max(max(dpdtdynamics.meanZ))],'pupil area (au)')
caxis([cmin/2 cmax])
ylim([-1.5 0]);
xlabel('Phase');ylabel('Amplitude')
axis xy

subplot(2,2,4);
a = imagesc([KOBehavior.dpdtdynamics.Xbins(:,:,1) KOBehavior.dpdtdynamics.Xbins(:,:,1)+2*pi],...
    KOBehavior.dpdtdynamics.Ybins(:,:,1),[nanmean(KOBehavior.dpdtdynamics.meanZ,3);...
    nanmean(KOBehavior.dpdtdynamics.meanZ,3)]');
colormap(gca,'jet')
alpha(a,double(~isnan([nanmean(KOBehavior.dpdtdynamics.meanZ,3);...
    nanmean(KOBehavior.dpdtdynamics.meanZ,3)]')))
%ColorbarWithAxis([min(min(dpdtdynamics.meanZ)) max(max(dpdtdynamics.meanZ))],'pupil area (au)')
caxis([cmin/2 cmax])
ylim([-1.5 0]);
xlabel('Phase');ylabel('Amplitude')
axis xy

NiceSave('Pupil_PupSpace2',figfolder,baseName)

%% FIGURE:
figure;
subplot(2,2,1);
bar(WTBehavior.EMGhist.logbins(:,:,1),...
    nanmean(WTBehavior.EMGhist.logcounts,3),...
    'facecolor','k','facealpha',0.85); hold on;
h1 = errorbar(WTBehavior.EMGhist.logbins(:,:,1),...
    nanmean(WTBehavior.EMGhist.logcounts,3),...
    nanstd(WTBehavior.EMGhist.logcounts,0,3),'k.');
bar(KOBehavior.EMGhist.logbins(:,:,1),...
    nanmean(KOBehavior.EMGhist.logcounts,3),...
    'facecolor','r','facealpha',0.65); hold on;
h2 = errorbar(KOBehavior.EMGhist.logbins(:,:,1),...
    nanmean(KOBehavior.EMGhist.logcounts,3),...
    nanstd(KOBehavior.EMGhist.logcounts,0,3),'r.');
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
LogScale('x',10);
axis tight
ylim([0 0.25]);
xlabel('EMG (norm.)'); ylabel('counts (norm.)');
legend({'floxed M1/M3','Emx1;M1/M3'},'location','northeast');

subplot(2,2,2);
shadederror(WTBehavior.Whdurhist.bins(:,:,1),...
    nanmean(WTBehavior.Whdurhist.Whdurs,3)',...
    nanstd(WTBehavior.Whdurhist.Whdurs,0,3)','k','none',1,'k');
shadederror(KOBehavior.Whdurhist.bins(:,:,1),...
    nanmean(KOBehavior.Whdurhist.Whdurs,3)',...
    nanstd(KOBehavior.Whdurhist.Whdurs,0,3)','r','none',1,'r');
LogScale('x',10)
%legend({'floxed M1/M3','Emx1;M1/M3'},'location','northeast');
xlabel('duration (s)'); ylabel('counts (norm.)')
title('Wh durations')

subplot(2,2,4);
shadederror(WTBehavior.Whdurhist.bins(:,:,1),...
    nanmean(WTBehavior.Whdurhist.InterWhdurs,3)',...
    nanstd(WTBehavior.Whdurhist.InterWhdurs,0,3)','k','none',1,'k');
shadederror(KOBehavior.Whdurhist.bins(:,:,1),...
    nanmean(KOBehavior.Whdurhist.InterWhdurs,3)',...
    nanstd(KOBehavior.Whdurhist.InterWhdurs,0,3)','r','none',1,'r');
LogScale('x',10)
%legend({'floxed M1/M3','Emx1;M1/M3'},'location','northeast');
xlabel('duration (s)'); ylabel('counts (norm.)')
title('NWh durations')

subplot(2,2,3);
shadederror(log10(WTBehavior.EMGPSD.freqs(:,:,1)),nanmean(WTBehavior.EMGPSD.psd,3)',...
    nanstd(WTBehavior.EMGPSD.psd,0,3)','k','none',1,'k');
shadederror(log10(KOBehavior.EMGPSD.freqs(:,:,1)),nanmean(KOBehavior.EMGPSD.psd,3)',...
    nanstd(KOBehavior.EMGPSD.psd,0,3)','r','none',1,'r');
LogScale('x',10)
xlabel('f (Hz)'); ylabel('power (dB)')
axis tight
ylim([-1 2.5]);
title('EMG power spectrum');

NiceSave('EMGStats',figfolder,baseName)

%% FIGURE:
cmin = min([min(min(nanmean(WTBehavior.pupilEMGdist.counts,3)))...
    min(min(nanmean(KOBehavior.pupilEMGdist.counts,3)))]);
cmax = max([max(max(nanmean(WTBehavior.pupilEMGdist.counts,3)))...
    max(max(nanmean(KOBehavior.pupilEMGdist.counts,3)))]);

figure;
subplot(2,3,1);
colormap(gca,[1 1 1;colormap('jet')]);
imagesc(WTBehavior.pupilEMGdist.bins{:,1,1},WTBehavior.pupilEMGdist.bins{:,2,1},...
    nanmean(WTBehavior.pupilEMGdist.counts,3)'); hold on
%plot(get(gca,'xlim'),log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],'w--')
LogScale('y',10); axis xy
%ColorbarWithAxis([0 0.6],['counts (au)'])
caxis([cmin cmax/1.5]);
xlabel('Pupil diameter (med norm.)'); ylabel('EMG');

subplot(2,3,4);
colormap(gca,[1 1 1;colormap('jet')]);
imagesc(KOBehavior.pupilEMGdist.bins{:,1,1},KOBehavior.pupilEMGdist.bins{:,2,1},...
    nanmean(KOBehavior.pupilEMGdist.counts,3)'); hold on
%plot(get(gca,'xlim'),log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],'w--')
LogScale('y',10); axis xy
%ColorbarWithAxis([0 0.6],['counts (au)'])
caxis([cmin cmax/1.5]);
xlabel('Pupil diameter (med norm.)'); ylabel('EMG');

temp = nanmean(WTBehavior.pupildynamicsEMG.meanZ,3);
temp(nanmean(WTBehavior.pupildynamicsEMG.N,3)<2) = NaN;
subplot(2,3,2);
imagesc(WTBehavior.pupildynamicsEMG.Xbins(:,:,1),WTBehavior.pupildynamicsEMG.Ybins(:,:,1),...
    temp');
colormap(gca,[1 1 1; colormap('jet')])
axis xy
xlim([-0.5 0.5]);ylim([-0.5 0.5])
%colorbar
ylabel('dP/dt')
title('EMG')
hold on
%caxis([min(min(pupildynamicsEMG.meanZ)) max(max(pupildynamicsEMG.meanZ))])
%LogScale('c',10)
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color','k')

temp = nanmean(KOBehavior.pupildynamicsEMG.meanZ,3);
temp(nanmean(KOBehavior.pupildynamicsEMG.N,3)<2) = NaN;
subplot(2,3,5);
imagesc(KOBehavior.pupildynamicsEMG.Xbins(:,:,1),KOBehavior.pupildynamicsEMG.Ybins(:,:,1),...
    temp');
colormap(gca,[1 1 1; colormap('jet')])
axis xy
xlim([-0.5 0.5]);ylim([-0.5 0.5])
%colorbar
ylabel('dP/dt')
title('EMG')
hold on
%caxis([min(min(pupildynamicsEMG.meanZ)) max(max(pupildynamicsEMG.meanZ))])
%LogScale('c',10)
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color','k')

temp = nanmean(WTBehavior.pupildynamicsEMG.pWhisk,3);
temp(nanmean(WTBehavior.pupildynamicsEMG.N,3)<2) = NaN;
subplot(2,3,3);
imagesc(WTBehavior.pupildynamicsEMG.Xbins(:,:,1),WTBehavior.pupildynamicsEMG.Ybins(:,:,1),...
    temp');
colormap(gca,[1 1 1; colormap('jet')])
axis xy
xlim([-0.5 0.5]);ylim([-0.5 0.5])
caxis([-0.01 1])
%colorbar
ylabel('dP/dt')
title('p(Whisking)')
hold on
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color','k')

temp = nanmean(KOBehavior.pupildynamicsEMG.pWhisk,3);
temp(nanmean(KOBehavior.pupildynamicsEMG.N,3)<2) = NaN;
subplot(2,3,6);
imagesc(KOBehavior.pupildynamicsEMG.Xbins(:,:,1),KOBehavior.pupildynamicsEMG.Ybins(:,:,1),...
    temp');
colormap(gca,[1 1 1; colormap('jet')])
axis xy
xlim([-0.5 0.5]);ylim([-0.5 0.5])
caxis([-0.01 1])
%colorbar
ylabel('dP/dt')
title('p(Whisking)')
hold on
plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color','k')

NiceSave('EMG_PupSpace1',figfolder,baseName);

%%
X = WTBehaviorInts.whints_pupilphase(:,1);
Y = log10(WTBehaviorInts.whints_pupilamp(:,1));
Y = Y(~isnan(X));
Z = WTBehaviorInts.whints_durs(:,1);
Z = Z(~isnan(X));
X = X(~isnan(X));

[WTpupilphaseEMG.meanWhdur,WTpupilphaseEMG.numWhstarts,WTpupilphaseEMG.Xbins,WTpupilphaseEMG.Ybins] = ConditionalHist3( X,Y,Z,...
    'minXY',0.5,'Xbounds',[-pi pi],'Ybounds',[-3 1],...
    'sigma',0.05,'numXbins',100,'numYbins',100);

WTpupilphaseEMG.meanWhdur(WTpupilphaseEMG.numWhstarts < 0.02) = NaN;
WTpupilphaseEMG.numWhstarts(WTpupilphaseEMG.numWhstarts < 0.02) = NaN;
WTpupilphaseEMG.Whstartrate = WTpupilphaseEMG.numWhstarts./nanmean(WTBehavior.pupilphaseEMG.occupancy,3);

X = KOBehaviorInts.whints_pupilphase(:,1);
Y = log10(KOBehaviorInts.whints_pupilamp(:,1));
Y = Y(~isnan(X));
Z = KOBehaviorInts.whints_durs(:,1);
Z = Z(~isnan(X));
X = X(~isnan(X));

[KOpupilphaseEMG.meanWhdur,KOpupilphaseEMG.numWhstarts,KOpupilphaseEMG.Xbins,KOpupilphaseEMG.Ybins] = ConditionalHist3( X,Y,Z,...
    'minXY',0.5,'Xbounds',[-pi pi],'Ybounds',[-3 1],...
    'sigma',0.05,'numXbins',100,'numYbins',100);

KOpupilphaseEMG.meanWhdur(KOpupilphaseEMG.numWhstarts < 0.02) = NaN;
KOpupilphaseEMG.numWhstarts(KOpupilphaseEMG.numWhstarts < 0.02) = NaN;
KOpupilphaseEMG.Whstartrate = KOpupilphaseEMG.numWhstarts./nanmean(KOBehavior.pupilphaseEMG.occupancy,3);

%% FIGURE:
cmin = min([min(min(log10(WTpupilphaseEMG.meanWhdur)))...
    min(min(log10(KOpupilphaseEMG.meanWhdur)))]);
cmax = max([max(max(log10(WTpupilphaseEMG.meanWhdur)))...
    max(max(log10(KOpupilphaseEMG.meanWhdur)))]);

figure;
subplot(2,2,1);
a = imagesc([WTpupilphaseEMG.Xbins WTpupilphaseEMG.Xbins+2*pi],...
    WTpupilphaseEMG.Ybins,log10([WTpupilphaseEMG.meanWhdur; WTpupilphaseEMG.meanWhdur])');
colormap(gca,'jet')
alpha(a,double(~isnan([WTpupilphaseEMG.meanWhdur; WTpupilphaseEMG.meanWhdur]')))
%ColorbarWithAxis(log10([min(min(pupilphaseEMG.meanWhdur)) max(max(pupilphaseEMG.meanWhdur))]),'EMG Envelope')
ylim([-1.5 0]);  xlim([-pi 3*pi])
caxis([cmin cmax])
xlabel('Phase');ylabel('Amplitude')
colorbar
axis xy
LogScale('c',10)
title('Whisk duration (s)')

subplot(2,2,3);
a = imagesc([KOpupilphaseEMG.Xbins KOpupilphaseEMG.Xbins+2*pi],...
    KOpupilphaseEMG.Ybins,log10([KOpupilphaseEMG.meanWhdur; KOpupilphaseEMG.meanWhdur])');
colormap(gca,'jet')
alpha(a,double(~isnan([KOpupilphaseEMG.meanWhdur; KOpupilphaseEMG.meanWhdur]')))
%ColorbarWithAxis(log10([min(min(pupilphaseEMG.meanWhdur)) max(max(pupilphaseEMG.meanWhdur))]),'EMG Envelope')
ylim([-1.5 0]);  xlim([-pi 3*pi])
caxis([cmin cmax])
LogScale('c',10)
xlabel('Phase');ylabel('Amplitude')
colorbar
axis xy
LogScale('c',10)
title('Whisk duration (s)')

cmin = min([min(min(log10(WTpupilphaseEMG.Whstartrate)))...
    min(min(log10(KOpupilphaseEMG.Whstartrate)))]);
cmax = max([max(max(log10(WTpupilphaseEMG.Whstartrate)))...
    max(max(log10(KOpupilphaseEMG.Whstartrate)))]);

subplot(2,2,2);
a = imagesc([WTpupilphaseEMG.Xbins WTpupilphaseEMG.Xbins+2*pi],...
    WTpupilphaseEMG.Ybins,log10([WTpupilphaseEMG.Whstartrate; WTpupilphaseEMG.Whstartrate])');
colormap(gca,'jet')
alpha(a,double(~isnan([WTpupilphaseEMG.Whstartrate; WTpupilphaseEMG.Whstartrate]')))
%ColorbarWithAxis(log10([min(min(pupilphaseEMG.Whstartrate)) max(max(pupilphaseEMG.Whstartrate))]),'EMG Envelope')
ylim([-1.5 0]); xlim([-pi 3*pi])
xlabel('Phase');ylabel('Amplitude')
colorbar
axis xy
caxis([cmin 1.5])
LogScale('c',10)
title('Whisk start rate')

subplot(2,2,4);
a = imagesc([KOpupilphaseEMG.Xbins KOpupilphaseEMG.Xbins+2*pi],...
    KOpupilphaseEMG.Ybins,log10([KOpupilphaseEMG.Whstartrate; KOpupilphaseEMG.Whstartrate])');
colormap(gca,'jet')
alpha(a,double(~isnan([KOpupilphaseEMG.Whstartrate; KOpupilphaseEMG.Whstartrate]')))
%ColorbarWithAxis(log10([min(min(pupilphaseEMG.Whstartrate)) max(max(pupilphaseEMG.Whstartrate))]),'EMG Envelope')
ylim([-1.5 0]); xlim([-pi 3*pi])
xlabel('Phase');ylabel('Amplitude')
colorbar
axis xy
caxis([cmin 1.5])
LogScale('c',10)
title('Whisk start rate')

NiceSave('EMG_PupSpace2',figfolder,baseName);

%% FIGURE:
cmin = min([min(min(nanmean(WTBehavior.pupilphaseEMG.meanZ,3)))...
    min(min(nanmean(KOBehavior.pupilphaseEMG.meanZ,3)))]);
cmax = max([max(max(nanmean(WTBehavior.pupilphaseEMG.meanZ,3)))...
    max(max(nanmean(KOBehavior.pupilphaseEMG.meanZ,3)))]);

figure;
subplot(2,2,1);
a = imagesc([WTBehavior.pupilphaseEMG.Xbins(:,:,1) WTBehavior.pupilphaseEMG.Xbins(:,:,1)+2*pi],...
    WTBehavior.pupilphaseEMG.Ybins(:,:,1),...
    [nanmean(WTBehavior.pupilphaseEMG.meanZ,3); nanmean(WTBehavior.pupilphaseEMG.meanZ,3)]');
colormap(gca,'jet')
alpha(a,double(~isnan([nanmean(WTBehavior.pupilphaseEMG.meanZ,3); nanmean(WTBehavior.pupilphaseEMG.meanZ,3)]')))
colorbar
%ColorbarWithAxis([min(min(pupilphaseEMG.meanZ)) max(max(pupilphaseEMG.meanZ))],'EMG')
caxis([cmin cmax])
ylim([-1.5 0]);
xlabel('Phase');ylabel('Amplitude')
title('EMG')
axis xy

subplot(2,2,3);
a = imagesc([KOBehavior.pupilphaseEMG.Xbins(:,:,1) KOBehavior.pupilphaseEMG.Xbins(:,:,1)+2*pi],...
    KOBehavior.pupilphaseEMG.Ybins(:,:,1),...
    [nanmean(KOBehavior.pupilphaseEMG.meanZ,3); nanmean(KOBehavior.pupilphaseEMG.meanZ,3)]');
colormap(gca,'jet')
alpha(a,double(~isnan([nanmean(KOBehavior.pupilphaseEMG.meanZ,3); nanmean(KOBehavior.pupilphaseEMG.meanZ,3)]')))
colorbar
%ColorbarWithAxis([min(min(pupilphaseEMG.meanZ)) max(max(pupilphaseEMG.meanZ))],'EMG')
caxis([cmin cmax])
ylim([-1.5 0]);
xlabel('Phase');ylabel('Amplitude')
title('EMG')
axis xy

cmin = min([min(min(nanmean(WTBehavior.pupilphaseEMG.pWhisk,3)))...
    min(min(nanmean(KOBehavior.pupilphaseEMG.pWhisk,3)))]);
cmax = max([max(max(nanmean(WTBehavior.pupilphaseEMG.pWhisk,3)))...
    max(max(nanmean(KOBehavior.pupilphaseEMG.pWhisk,3)))]);

subplot(2,2,2);
a = imagesc([WTBehavior.pupilphaseEMG.Xbins(:,:,1) WTBehavior.pupilphaseEMG.Xbins(:,:,1)+2*pi],...
    WTBehavior.pupilphaseEMG.Ybins(:,:,1),...
    [nanmean(WTBehavior.pupilphaseEMG.pWhisk,3); nanmean(WTBehavior.pupilphaseEMG.pWhisk,3)]');
colormap(gca,'jet')
alpha(a,double(~isnan([nanmean(WTBehavior.pupilphaseEMG.pWhisk,3); nanmean(WTBehavior.pupilphaseEMG.pWhisk,3)]')))
colorbar
%ColorbarWithAxis([min(min(pupilphaseEMG.pWhisk)) max(max(pupilphaseEMG.pWhisk))],'EMG')
%caxis([min(min(pupilphaseEMG.pWhisk)) max(max(pupilphaseEMG.pWhisk))])
ylim([-1.5 0]);
xlabel('Phase');ylabel('Amplitude')
title('EMG')
axis xy

subplot(2,2,4);
a = imagesc([KOBehavior.pupilphaseEMG.Xbins(:,:,1) KOBehavior.pupilphaseEMG.Xbins(:,:,1)+2*pi],...
    KOBehavior.pupilphaseEMG.Ybins(:,:,1),...
    [nanmean(KOBehavior.pupilphaseEMG.pWhisk,3); nanmean(KOBehavior.pupilphaseEMG.pWhisk,3)]');
colormap(gca,'jet')
alpha(a,double(~isnan([nanmean(KOBehavior.pupilphaseEMG.pWhisk,3); nanmean(KOBehavior.pupilphaseEMG.pWhisk,3)]')))
colorbar
%ColorbarWithAxis([min(min(pupilphaseEMG.pWhisk)) max(max(pupilphaseEMG.pWhisk))],'EMG')
%caxis([min(min(pupilphaseEMG.pWhisk)) max(max(pupilphaseEMG.pWhisk))])
ylim([-1.5 0]);
xlabel('Phase');ylabel('Amplitude')
title('EMG')
axis xy

NiceSave('EMG_PupSpace3',figfolder,baseName);

%% FIGURE:
lagwin = [-5 5];

figure;
subplot(3,2,1);
shadederror(WTBehavior.pupACG.tlag(:,:,1),...
    nanmean(WTBehavior.pupACG.ACG,3)',...
    nanstd(WTBehavior.pupACG.ACG,0,3)','k','none',1,'k');
shadederror(KOBehavior.pupACG.tlag(:,:,1),...
    nanmean(KOBehavior.pupACG.ACG,3)',...
    nanstd(KOBehavior.pupACG.ACG,0,3)','r','none',1,'r');
axis tight
%plot(get(gca,'xlim'),[0 0],'r')
xlabel('lag (s)')
ylabel('autocovariance')
title('Pupil diameter')

subplot(3,2,2);
shadederror(WTBehavior.pupilEMGcorr.corrlags(:,:,1),...
    nanmean(WTBehavior.pupilEMGcorr.xcorr,3)',...
    nanstd(WTBehavior.pupilEMGcorr.xcorr,0,3)','k','none',1,'k');
shadederror(KOBehavior.pupilEMGcorr.corrlags(:,:,1),...
    nanmean(KOBehavior.pupilEMGcorr.xcorr,3)',...
    nanstd(KOBehavior.pupilEMGcorr.xcorr,0,3)','r','none',1,'r');
axis tight
plot([0 0],get(gca,'ylim'),'k-','linewidth',1)
xlabel('t lag (s)');ylabel('EMG')
title('Pupil-EMG CCG')

subplot(3,2,3);
shadederror(WTBehavior.pwCCG.t_lag(:,:,1),...
    nanmean(WTBehavior.pwCCG.pupil.WhOn-min(min(WTBehavior.pwCCG.pupil.WhOn)),3)',...
    nanstd(WTBehavior.pwCCG.pupil.WhOn-min(min(WTBehavior.pwCCG.pupil.WhOn)),0,3)','k','none',1,'k');
shadederror(KOBehavior.pwCCG.t_lag(:,:,1),...
    nanmean(KOBehavior.pwCCG.pupil.WhOn-min(min(KOBehavior.pwCCG.pupil.WhOn)),3)',...
    nanstd(KOBehavior.pwCCG.pupil.WhOn-min(min(KOBehavior.pwCCG.pupil.WhOn)),0,3)','r','none',1,'r');
axis tight
plot([0 0],get(gca,'ylim'),'k-','linewidth',1)
xlim(lagwin)
xlabel('t lag (s)');ylabel('Pupil')

subplot(3,2,5);
shadederror(WTBehavior.pwCCG.t_lag(:,:,1),...
    nanmean(WTBehavior.pwCCG.EMG.WhOn,3)',...
    nanstd(WTBehavior.pwCCG.EMG.WhOn,0,3)','k','none',1,'k');
shadederror(KOBehavior.pwCCG.t_lag(:,:,1),...
    nanmean(KOBehavior.pwCCG.EMG.WhOn,3)',...
    nanstd(KOBehavior.pwCCG.EMG.WhOn,0,3)','r','none',1,'r');
axis tight
plot([0 0],get(gca,'ylim'),'k-','linewidth',1)
xlim(lagwin)
xlabel('t (s - aligned to WhOn)')
ylabel('EMG')

subplot(3,2,4);
shadederror(WTBehavior.wpCCG.t_lag(:,:,1),...
    nanmean(WTBehavior.wpCCG.pupil.PupOn,3)',...
    nanstd(WTBehavior.wpCCG.pupil.PupOn,0,3)','k','none',1,'k');
shadederror(KOBehavior.wpCCG.t_lag(:,:,1),...
    nanmean(KOBehavior.wpCCG.pupil.PupOn,3)',...
    nanstd(KOBehavior.wpCCG.pupil.PupOn,0,3)','r','none',1,'r');
axis tight
plot([0 0],get(gca,'ylim'),'k-','linewidth',1)
%xlim(lagwin)
ylabel('Pupil')

subplot(3,2,6);
shadederror(WTBehavior.wpCCG.t_lag(:,:,1),...
    nanmean(WTBehavior.wpCCG.EMG.PupOn,3)',...
    nanstd(WTBehavior.wpCCG.EMG.PupOn,0,3)','k','none',1,'k');
shadederror(KOBehavior.wpCCG.t_lag(:,:,1),...
    nanmean(KOBehavior.wpCCG.EMG.PupOn,3)',...
    nanstd(KOBehavior.wpCCG.EMG.PupOn,0,3)','r','none',1,'r');
axis tight
plot([0 0],get(gca,'ylim'),'k-','linewidth',1)
%xlim(lagwin)
xlabel('t (s - aligned to PupOn)')
ylabel('EMG')

NiceSave('Pupil_EMG_ACG_CCG',figfolder,baseName)

%% FIGURE:
rwbcolormap = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);
plotx = linspace(-pi,3*pi,100);

figure;
subplot(2,2,1);
shadederror(log10(WTBehavior.WPcoupling.freqs(:,:,1)),...
    nanmean(WTBehavior.WPcoupling.coupling,3),...
    nanstd(WTBehavior.WPcoupling.coupling,0,3),'k','none',1,'k');
shadederror(log10(KOBehavior.WPcoupling.freqs(:,:,1)),...
    nanmean(KOBehavior.WPcoupling.coupling,3),...
    nanstd(KOBehavior.WPcoupling.coupling,0,3),'r','none',1,'r');
LogScale('x',10);
xlabel('f (Hz)');
axis tight
title('Pupil-EMG phase coupling')

subplot(2,2,2);
imagesc(WTBehavior.PhaseAmpCoup.phasebins(:,:,1),WTBehavior.PhaseAmpCoup.ampbins(:,:,1),...
    nanmean(WTBehavior.PhaseAmpCoup.phaseamphist,3)); hold on;
imagesc(WTBehavior.PhaseAmpCoup.phasebins(:,:,1)+2*pi,WTBehavior.PhaseAmpCoup.ampbins(:,:,1),...
    nanmean(WTBehavior.PhaseAmpCoup.phaseamphist,3))
plot(circ_mean(WTBehavior.PhaseAmpCoup.sig2prefangle,[],3),WTBehavior.PhaseAmpCoup.ampbins(:,:,1),'.k')
plot(circ_mean(WTBehavior.PhaseAmpCoup.sig2prefangle,[],3)+2*pi,WTBehavior.PhaseAmpCoup.ampbins(:,:,1),'.k')
plot(plotx,cos(plotx),'k')
colormap(gca,rwbcolormap)
axis xy
axis tight
ColorbarWithAxis([-0.5 0.5],['Phase-EMG dist'])
caxis([-0.5 0.5])
xlim([-pi 3*pi]); 
xlabel('Pupil phase'); ylabel('Pupil (Z)')
title('WT')

subplot(2,2,4)
imagesc(KOBehavior.PhaseAmpCoup.phasebins(:,:,1),KOBehavior.PhaseAmpCoup.ampbins(:,:,1),...
    nanmean(KOBehavior.PhaseAmpCoup.phaseamphist,3)); hold on;
imagesc(KOBehavior.PhaseAmpCoup.phasebins(:,:,1)+2*pi,KOBehavior.PhaseAmpCoup.ampbins(:,:,1),...
    nanmean(KOBehavior.PhaseAmpCoup.phaseamphist,3))
plot(circ_mean(KOBehavior.PhaseAmpCoup.sig2prefangle,[],3),KOBehavior.PhaseAmpCoup.ampbins(:,:,1),'.k')
plot(circ_mean(KOBehavior.PhaseAmpCoup.sig2prefangle,[],3)+2*pi,KOBehavior.PhaseAmpCoup.ampbins(:,:,1),'.k')
plot(plotx,cos(plotx),'k')
colormap(gca,rwbcolormap)
axis xy
axis tight
ColorbarWithAxis([-0.5 0.5],['Phase-EMG dist'])
caxis([-0.5 0.5])
xlim([-pi 3*pi]); 
xlabel('Pupil phase'); ylabel('Pupil (Z)')
title('KO')

subplot(2,2,3)
shadederror(WTBehavior.PhaseAmpCoup.ampbins(:,:,1),...
    nanmean(WTBehavior.PhaseAmpCoup.sig2powerskew,3),...
    nanstd(WTBehavior.PhaseAmpCoup.sig2powerskew,0,3),'k','none',1,'k');
shadederror(KOBehavior.PhaseAmpCoup.ampbins(:,:,1),...
    nanmean(KOBehavior.PhaseAmpCoup.sig2powerskew,3),...
    nanstd(KOBehavior.PhaseAmpCoup.sig2powerskew,0,3),'r','none',1,'r');
xlabel('Pupil (Z)');
ylabel('Phase-EMG MI (mrl)')
axis tight
xlim([-2 2]);

NiceSave('Pupil_EMG_PhaseAmpCoup',figfolder,baseName)
