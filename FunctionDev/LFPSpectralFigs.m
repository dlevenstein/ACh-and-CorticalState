%%
basePath = pwd;
baseName = 'EM1M3';

figfolder = fullfile(basePath,'SummaryFigures');

load(fullfile(basePath,['WT_EM1M3.LFPSpectralAnalysis.mat']));
WTLFP = groupLFPSpectral;
WTDepth = groupDepth;
load(fullfile(basePath,['KO_EM1M3.LFPSpectralAnalysis.mat']));
KOLFP = groupLFPSpectral;
KODepth = groupDepth;

%%
numberchans = 100;

% Wh/NWh
[WTLayerSpec_NWh,normcolumn] = normCx(WTDepth.ndepth,WTLFP.cLayerSpec_NWh,numberchans);
[WTLayerSpec_Wh,normcolumn] = normCx(WTDepth.ndepth,WTLFP.cLayerSpec_Wh,numberchans);
[KOLayerSpec_NWh,normcolumn] = normCx(KODepth.ndepth,KOLFP.cLayerSpec_NWh,numberchans);
[KOLayerSpec_Wh,normcolumn] = normCx(KODepth.ndepth,KOLFP.cLayerSpec_Wh,numberchans);

% lo/hi Pupil
[WTLayerSpec_loP,normcolumn] = normCx(WTDepth.ndepth,WTLFP.cLayerSpec_loP,numberchans);
[WTLayerSpec_hiP,normcolumn] = normCx(WTDepth.ndepth,WTLFP.cLayerSpec_hiP,numberchans);
[KOLayerSpec_loP,normcolumn] = normCx(KODepth.ndepth,KOLFP.cLayerSpec_loP,numberchans);
[KOLayerSpec_hiP,normcolumn] = normCx(KODepth.ndepth,KOLFP.cLayerSpec_hiP,numberchans);

%% FIGURE:
freqs = logspace(log10(1),log10(128),100);
cmin = min([min(min(WTLayerSpec_Wh-WTLayerSpec_NWh))...
    min(min(KOLayerSpec_Wh-KOLayerSpec_NWh))]);
cmax = max([max(max(WTLayerSpec_Wh-WTLayerSpec_NWh))...
    max(max(KOLayerSpec_Wh-KOLayerSpec_NWh))]);

figure;
subplot(2,3,1);
imagesc(log10(freqs),normcolumn,(WTLayerSpec_Wh-WTLayerSpec_NWh)');
axis tight
LogScale('x',10)
colormap('jet');
%ColorbarWithAxis([cmin cmax],['power (dB)'])
caxis([cmax*-1 cmax]);
%set(gca,'YDir','reverse');
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('WT PSpec Wh-NWh diff');

subplot(2,3,4);
imagesc(log10(freqs),normcolumn,(KOLayerSpec_Wh-KOLayerSpec_NWh)');
axis tight
LogScale('x',10)
colormap('jet');
%ColorbarWithAxis([cmin cmax],['power (dB)'])
caxis([cmax*-1 cmax]);
%set(gca,'YDir','reverse');
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('KO PSpec Wh-NWh diff');

cmin = min([min(min(WTLayerSpec_Wh))...
    min(min(WTLayerSpec_NWh))...
    min(min(KOLayerSpec_Wh))...
    min(min(KOLayerSpec_NWh))]);
cmax = max([max(max(WTLayerSpec_Wh))...
    max(max(WTLayerSpec_NWh))...
    max(max(KOLayerSpec_Wh))...
    max(max(KOLayerSpec_NWh))]);

subplot(2,3,2);
imagesc(log10(freqs),normcolumn,WTLayerSpec_NWh');
axis tight
LogScale('x',10)
colormap('jet');
%ColorbarWithAxis([cmin cmax],['power (dB)'])
caxis([cmax*-1 cmax]);
%set(gca,'YDir','reverse');
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('Power spectra NWh');

subplot(2,3,3);
imagesc(log10(freqs),normcolumn,WTLayerSpec_Wh');
axis tight
LogScale('x',10)
colormap('jet');
%ColorbarWithAxis([cmin cmax],['power (dB)'])
caxis([cmax*-1 cmax]);
%set(gca,'YDir','reverse');
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('Power spectra Wh');

subplot(2,3,5);
imagesc(log10(freqs),normcolumn,KOLayerSpec_NWh');
axis tight
LogScale('x',10)
colormap('jet');
%ColorbarWithAxis([cmin cmax],['power (dB)'])
caxis([cmax*-1 cmax]);
%set(gca,'YDir','reverse');
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('Power spectra NWh');

subplot(2,3,6);
imagesc(log10(freqs),normcolumn,KOLayerSpec_Wh');
axis tight
LogScale('x',10)
colormap('jet');
%ColorbarWithAxis([cmin cmax],['power (dB)'])
caxis([cmax*-1 cmax]);
%set(gca,'YDir','reverse');
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('Power spectra Wh');

%NiceSave('LaminarPspec_Wh_NWh',figfolder,baseName)

%% FIGURE:
cmin = min([min(min(WTLayerSpec_hiP-WTLayerSpec_loP))...
    min(min(KOLayerSpec_hiP-KOLayerSpec_loP))]);
cmax = max([max(max(WTLayerSpec_hiP-WTLayerSpec_loP))...
    max(max(KOLayerSpec_hiP-KOLayerSpec_loP))]);

figure;
subplot(2,3,1);
imagesc(log10(freqs),normcolumn,(WTLayerSpec_hiP-WTLayerSpec_loP)');
axis tight
LogScale('x',10)
colormap('jet');
%ColorbarWithAxis([cmin cmax],['power (dB)'])
caxis([cmax*-1 cmax]);
%set(gca,'YDir','reverse');
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('WT PSpec Wh-NWh diff');

subplot(2,3,4);
imagesc(log10(freqs),normcolumn,(KOLayerSpec_hiP-KOLayerSpec_loP)');
axis tight
LogScale('x',10)
colormap('jet');
%ColorbarWithAxis([cmin cmax],['power (dB)'])
caxis([cmax*-1 cmax]);
%set(gca,'YDir','reverse');
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('KO PSpec Wh-NWh diff');

cmin = min([min(min(WTLayerSpec_hiP))...
    min(min(WTLayerSpec_loP))...
    min(min(KOLayerSpec_hiP))...
    min(min(KOLayerSpec_loP))]);
cmax = max([max(max(WTLayerSpec_hiP))...
    max(max(WTLayerSpec_loP))...
    max(max(KOLayerSpec_hiP))...
    max(max(KOLayerSpec_loP))]);

subplot(2,3,2);
imagesc(log10(freqs),normcolumn,WTLayerSpec_loP');
axis tight
LogScale('x',10)
colormap('jet');
%ColorbarWithAxis([cmin cmax],['power (dB)'])
caxis([cmax*-1 cmax]);
%set(gca,'YDir','reverse');
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('Power spectra NWh');

subplot(2,3,3);
imagesc(log10(freqs),normcolumn,WTLayerSpec_hiP');
axis tight
LogScale('x',10)
colormap('jet');
%ColorbarWithAxis([cmin cmax],['power (dB)'])
caxis([cmax*-1 cmax]);
%set(gca,'YDir','reverse');
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('Power spectra Wh');

subplot(2,3,5);
imagesc(log10(freqs),normcolumn,KOLayerSpec_loP');
axis tight
LogScale('x',10)
colormap('jet');
%ColorbarWithAxis([cmin cmax],['power (dB)'])
caxis([cmax*-1 cmax]);
%set(gca,'YDir','reverse');
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('Power spectra NWh');

subplot(2,3,6);
imagesc(log10(freqs),normcolumn,KOLayerSpec_hiP');
axis tight
LogScale('x',10)
colormap('jet');
%ColorbarWithAxis([cmin cmax],['power (dB)'])
caxis([cmax*-1 cmax]);
%set(gca,'YDir','reverse');
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);
title('Power spectra Wh');

%NiceSave('LaminarPspec_lo_hiPup',figfolder,baseName)

%%
wavespec.samplingRate = 250;
twin = [0.75 0.75].*wavespec.samplingRate;

WTL1eventSpec_Wh = nanmean(WTLFP.L1eventSpec_Wh,3);
WTL23eventSpec_Wh = nanmean(WTLFP.L23eventSpec_Wh,3);
WTL4eventSpec_Wh = nanmean(WTLFP.L4eventSpec_Wh,3);
WTL5aeventSpec_Wh = nanmean(WTLFP.L5aeventSpec_Wh,3);
WTL56eventSpec_Wh = nanmean(WTLFP.L56eventSpec_Wh,3);
WTL6eventSpec_Wh = nanmean(WTLFP.L6eventSpec_Wh,3);

KOL1eventSpec_Wh = nanmean(KOLFP.L1eventSpec_Wh,3);
KOL23eventSpec_Wh = nanmean(KOLFP.L23eventSpec_Wh,3);
KOL4eventSpec_Wh = nanmean(KOLFP.L4eventSpec_Wh,3);
KOL5aeventSpec_Wh = nanmean(KOLFP.L5aeventSpec_Wh,3);
KOL56eventSpec_Wh = nanmean(KOLFP.L56eventSpec_Wh,3);
KOL6eventSpec_Wh = nanmean(KOLFP.L6eventSpec_Wh,3);

%% FIGURE:
taxis = (-(twin(1)/wavespec.samplingRate):(1/wavespec.samplingRate):(twin(2)/wavespec.samplingRate))*1e3;

cmax = max([max(max(WTL1eventSpec_Wh)) max(max(WTL23eventSpec_Wh))...
    max(max(WTL4eventSpec_Wh)) max(max(WTL5aeventSpec_Wh))...
    max(max(WTL56eventSpec_Wh)) max(max(WTL6eventSpec_Wh))...
    max(max(KOL1eventSpec_Wh)) max(max(KOL23eventSpec_Wh))...
    max(max(KOL4eventSpec_Wh)) max(max(KOL5aeventSpec_Wh))...
    max(max(KOL56eventSpec_Wh)) max(max(KOL6eventSpec_Wh))]);

figure;
subplot(6,2,1);
imagesc(taxis,log10(freqs),WTL1eventSpec_Wh');hold on;
colormap jet;
LogScale('y',10);
caxis([-cmax cmax]);
axis xy
xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(freqs(1)) log10(freqs(end))],'--k');hold on;
title('WT L1');
    
subplot(6,2,2);
imagesc(taxis,log10(freqs),KOL1eventSpec_Wh');hold on;
colormap jet;
LogScale('y',10);
caxis([-cmax cmax]);
axis xy
xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(freqs(1)) log10(freqs(end))],'--k');hold on;
title('KO L1');

subplot(6,2,3);
imagesc(taxis,log10(freqs),WTL23eventSpec_Wh');hold on;
colormap jet;
LogScale('y',10);
caxis([-cmax cmax]);
axis xy
xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(freqs(1)) log10(freqs(end))],'--k');hold on;
title('L2/3');
    
subplot(6,2,4);
imagesc(taxis,log10(freqs),KOL23eventSpec_Wh');hold on;
colormap jet;
LogScale('y',10);
caxis([-cmax cmax]);
axis xy
xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(freqs(1)) log10(freqs(end))],'--k');hold on;
title('L2/3');

subplot(6,2,5);
imagesc(taxis,log10(freqs),WTL4eventSpec_Wh');hold on;
colormap jet;
LogScale('y',10);
caxis([-cmax cmax]);
axis xy
xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(freqs(1)) log10(freqs(end))],'--k');hold on;
title('L4')

subplot(6,2,6);
imagesc(taxis,log10(freqs),KOL4eventSpec_Wh');hold on;
colormap jet;
LogScale('y',10);
caxis([-cmax cmax]);
axis xy
xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(freqs(1)) log10(freqs(end))],'--k');hold on;
title('L4')

subplot(6,2,7);
imagesc(taxis,log10(freqs),WTL5aeventSpec_Wh');hold on;
colormap jet;
LogScale('y',10);
caxis([-cmax cmax]);
axis xy
xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(freqs(1)) log10(freqs(end))],'--k');hold on;
title('L5a')
    
subplot(6,2,8);
imagesc(taxis,log10(freqs),KOL5aeventSpec_Wh');hold on;
colormap jet;
LogScale('y',10);
caxis([-cmax cmax]);
axis xy
xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(freqs(1)) log10(freqs(end))],'--k');hold on;
title('L5a')

subplot(6,2,9);
imagesc(taxis,log10(freqs),WTL56eventSpec_Wh');hold on;
colormap jet;
LogScale('y',10);
caxis([-cmax cmax]);
axis xy
xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(freqs(1)) log10(freqs(end))],'--k');hold on;
title('L5/6')
    
subplot(6,2,10);
imagesc(taxis,log10(freqs),KOL56eventSpec_Wh');hold on;
colormap jet;
LogScale('y',10);
caxis([-cmax cmax]);
axis xy
xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(freqs(1)) log10(freqs(end))],'--k');hold on;
title('L5/6')

subplot(6,2,11);
imagesc(taxis,log10(freqs),WTL6eventSpec_Wh');hold on;
colormap jet;
LogScale('y',10);
caxis([-cmax cmax]);
axis xy
xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(freqs(1)) log10(freqs(end))],'--k');hold on;
title('L6')
  
subplot(6,2,12);
imagesc(taxis,log10(freqs),KOL6eventSpec_Wh');hold on;
colormap jet;
LogScale('y',10);
caxis([-cmax cmax]);
axis xy
xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(freqs(1)) log10(freqs(end))],'--k');hold on;
title('L6')

%NiceSave('Laminar_eventSpec_Wh',figfolder,baseName)

%%
L1comodcorrs = NaN(length(freqs),length(freqs),6);
L1comodcorrs(:,:,1) = nanmean(WTLFP.L1comodcorrs_all.L1,3);
L1comodcorrs(:,:,2) = nanmean(WTLFP.L1comodcorrs_all.L23,3);
L1comodcorrs(:,:,3) = nanmean(WTLFP.L1comodcorrs_all.L4,3);
L1comodcorrs(:,:,4) = nanmean(WTLFP.L1comodcorrs_all.L5a,3);
L1comodcorrs(:,:,5) = nanmean(WTLFP.L1comodcorrs_all.L56,3);
L1comodcorrs(:,:,6) = nanmean(WTLFP.L1comodcorrs_all.L6,3);

L23comodcorrs = NaN(length(freqs),length(freqs),6);
L23comodcorrs(:,:,1) = nanmean(WTLFP.L23comodcorrs_all.L1,3);
L23comodcorrs(:,:,2) = nanmean(WTLFP.L23comodcorrs_all.L23,3);
L23comodcorrs(:,:,3) = nanmean(WTLFP.L23comodcorrs_all.L4,3);
L23comodcorrs(:,:,4) = nanmean(WTLFP.L23comodcorrs_all.L5a,3);
L23comodcorrs(:,:,5) = nanmean(WTLFP.L23comodcorrs_all.L56,3);
L23comodcorrs(:,:,6) = nanmean(WTLFP.L23comodcorrs_all.L6,3);

L4comodcorrs = NaN(length(freqs),length(freqs),6);
L4comodcorrs(:,:,1) = nanmean(WTLFP.L4comodcorrs_all.L1,3);
L4comodcorrs(:,:,2) = nanmean(WTLFP.L4comodcorrs_all.L23,3);
L4comodcorrs(:,:,3) = nanmean(WTLFP.L4comodcorrs_all.L4,3);
L4comodcorrs(:,:,4) = nanmean(WTLFP.L4comodcorrs_all.L5a,3);
L4comodcorrs(:,:,5) = nanmean(WTLFP.L4comodcorrs_all.L56,3);
L4comodcorrs(:,:,6) = nanmean(WTLFP.L4comodcorrs_all.L6,3);

L5acomodcorrs = NaN(length(freqs),length(freqs),6);
L5acomodcorrs(:,:,1) = nanmean(WTLFP.L5acomodcorrs_all.L1,3);
L5acomodcorrs(:,:,2) = nanmean(WTLFP.L5acomodcorrs_all.L23,3);
L5acomodcorrs(:,:,3) = nanmean(WTLFP.L5acomodcorrs_all.L4,3);
L5acomodcorrs(:,:,4) = nanmean(WTLFP.L5acomodcorrs_all.L5a,3);
L5acomodcorrs(:,:,5) = nanmean(WTLFP.L5acomodcorrs_all.L56,3);
L5acomodcorrs(:,:,6) = nanmean(WTLFP.L5acomodcorrs_all.L6,3);

L56comodcorrs = NaN(length(freqs),length(freqs),6);
L56comodcorrs(:,:,1) = nanmean(WTLFP.L56comodcorrs_all.L1,3);
L56comodcorrs(:,:,2) = nanmean(WTLFP.L56comodcorrs_all.L23,3);
L56comodcorrs(:,:,3) = nanmean(WTLFP.L56comodcorrs_all.L4,3);
L56comodcorrs(:,:,4) = nanmean(WTLFP.L56comodcorrs_all.L5a,3);
L56comodcorrs(:,:,5) = nanmean(WTLFP.L56comodcorrs_all.L56,3);
L56comodcorrs(:,:,6) = nanmean(WTLFP.L56comodcorrs_all.L6,3);

L6comodcorrs = NaN(length(freqs),length(freqs),6);
L6comodcorrs(:,:,1) = nanmean(WTLFP.L6comodcorrs_all.L1,3);
L6comodcorrs(:,:,2) = nanmean(WTLFP.L6comodcorrs_all.L23,3);
L6comodcorrs(:,:,3) = nanmean(WTLFP.L6comodcorrs_all.L4,3);
L6comodcorrs(:,:,4) = nanmean(WTLFP.L6comodcorrs_all.L5a,3);
L6comodcorrs(:,:,5) = nanmean(WTLFP.L6comodcorrs_all.L56,3);
L6comodcorrs(:,:,6) = nanmean(WTLFP.L6comodcorrs_all.L6,3);

%% FIGURE:

figure;
for i = 1:size(L1comodcorrs,3)
    subplot(6,6,i*1);
    imagesc(log10(freqs),log10(freqs),L1comodcorrs(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(L23comodcorrs,3)
    subplot(6,6,i+6);
    imagesc(log10(freqs),log10(freqs),L23comodcorrs(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(L4comodcorrs,3)
    subplot(6,6,i+12);
    imagesc(log10(freqs),log10(freqs),L4comodcorrs(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(L5acomodcorrs,3)
    subplot(6,6,i+18);
    imagesc(log10(freqs),log10(freqs),L5acomodcorrs(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(L56comodcorrs,3)
    subplot(6,6,i+24);
    imagesc(log10(freqs),log10(freqs),L56comodcorrs(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

for i = 1:size(L6comodcorrs,3)
    subplot(6,6,i+30);
    imagesc(log10(freqs),log10(freqs),L6comodcorrs(:,:,i))
    axis xy
    xlabel('f (Hz)'); ylabel('f (Hz)');
    colormap(gca,'jet')
    %caxis([cmin cmax])
    LogScale('x',10); LogScale('y',10);
end

%NiceSave('LaminarCoMOD_all',figfolder,baseName)
