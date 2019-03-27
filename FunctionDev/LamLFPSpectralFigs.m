%%
basePath = pwd;
baseName = 'EM1M3';

figfolder = fullfile(basePath,'SummaryFigures');
if (~exist(figfolder,'dir'))
    mkdir(figfolder)
end

load(fullfile(basePath,['WT_EM1M3.LamLFPSpectralAnalysis.mat']));
WTLamLFPSpectral = groupLamLFPSpectral;

load(fullfile(basePath,['KO_EM1M3.LamLFPSpectralAnalysis.mat']));
KOLamLFPSpectral = groupLamLFPSpectral;

%% FIGURE:
twin = round([0.75 0.75].*800);
taxis = (-(twin(1)/800):(1/800):(twin(2)/800))*1e3;

cmin = min([min(min(L1eventSpec.Wh.spec)) min(min(L23eventSpec.Wh.spec))...
    min(min(L4eventSpec.Wh.spec)) min(min(L5aeventSpec.Wh.spec))...
    min(min(L56eventSpec.Wh.spec)) min(min(L6eventSpec.Wh.spec))]);
cmax = max([max(max(L1eventSpec.Wh.spec)) max(max(L23eventSpec.Wh.spec))...
    max(max(L4eventSpec.Wh.spec)) max(max(L5aeventSpec.Wh.spec))...
    max(max(L56eventSpec.Wh.spec)) max(max(L6eventSpec.Wh.spec))]);

figure;
subplot(6,3,1);
imagesc(taxis,log10(wavespec.freqs),L1eventSpec.Wh.spec');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmax*-1 cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.freqs(1)) log10(wavespec.freqs(end))],'--k');hold on;
title('L1 modZ-PSpec');

subplot(6,3,4);
imagesc(taxis,log10(wavespec.freqs),L23eventSpec.Wh.spec');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmax*-1 cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.freqs(1)) log10(wavespec.freqs(end))],'--k');hold on;
title('L2/3');

subplot(6,3,7);
imagesc(taxis,log10(wavespec.freqs),L4eventSpec.Wh.spec');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmax*-1 cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.freqs(1)) log10(wavespec.freqs(end))],'--k');hold on;
title('L4');

subplot(6,3,10);
imagesc(taxis,log10(wavespec.freqs),L5aeventSpec.Wh.spec');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmax*-1 cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.freqs(1)) log10(wavespec.freqs(end))],'--k');hold on;
title('L5a');

subplot(6,3,13);
imagesc(taxis,log10(wavespec.freqs),L56eventSpec.Wh.spec');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmax*-1 cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.freqs(1)) log10(wavespec.freqs(end))],'--k');hold on;
title('L5b');

subplot(6,3,16);
imagesc(taxis,log10(wavespec.freqs),L6eventSpec.Wh.spec');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmax*-1 cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.freqs(1)) log10(wavespec.freqs(end))],'--k');hold on;
title('L6');

%%
% For Fractal
cmin = min([min(min(L1eventSpec.Wh.frac)) min(min(L23eventSpec.Wh.frac))...
    min(min(L4eventSpec.Wh.frac)) min(min(L5aeventSpec.Wh.frac))...
    min(min(L56eventSpec.Wh.frac)) min(min(L6eventSpec.Wh.frac))]);
cmax = max([max(max(L1eventSpec.Wh.frac)) max(max(L23eventSpec.Wh.frac))...
    max(max(L4eventSpec.Wh.frac)) max(max(L5aeventSpec.Wh.frac))...
    max(max(L56eventSpec.Wh.frac)) max(max(L6eventSpec.Wh.frac))]);

subplot(6,3,2);
imagesc(taxis,log10(wavespec.validfreq),L1eventSpec.Wh.frac');hold on;
colormap jet;
LogScale('y',10);
LogScale('c',10);
caxis([cmax/4*-1 cmax/4]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L1 Fractal');

subplot(6,3,5);
imagesc(taxis,log10(wavespec.validfreq),L23eventSpec.Wh.frac');hold on;
colormap jet;
LogScale('y',10);
LogScale('c',10);
caxis([cmax/4*-1 cmax/4]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L2/3');

subplot(6,3,8);
imagesc(taxis,log10(wavespec.validfreq),L4eventSpec.Wh.frac');hold on;
colormap jet;
LogScale('y',10);
LogScale('c',10);
caxis([cmax/4*-1 cmax/4]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L4');

subplot(6,3,11);
imagesc(taxis,log10(wavespec.validfreq),L5aeventSpec.Wh.frac');hold on;
colormap jet;
LogScale('y',10);
LogScale('c',10);
caxis([cmax/4*-1 cmax/4]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L5a');

subplot(6,3,14);
imagesc(taxis,log10(wavespec.validfreq),L56eventSpec.Wh.frac');hold on;
colormap jet;
LogScale('y',10);
LogScale('c',10);
caxis([cmax/4*-1 cmax/4]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L5b');

subplot(6,3,17);
imagesc(taxis,log10(wavespec.validfreq),L6eventSpec.Wh.frac');hold on;
colormap jet;
LogScale('y',10);
LogScale('c',10);
caxis([cmax/4*-1 cmax/4]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L6');

% for Osci
cmin = min([min(min(L1eventSpec.Wh.osci)) min(min(L23eventSpec.Wh.osci))...
    min(min(L4eventSpec.Wh.osci)) min(min(L5aeventSpec.Wh.osci))...
    min(min(L56eventSpec.Wh.osci)) min(min(L6eventSpec.Wh.osci))]);
cmax = max([max(max(L1eventSpec.Wh.osci)) max(max(L23eventSpec.Wh.osci))...
    max(max(L4eventSpec.Wh.osci)) max(max(L5aeventSpec.Wh.osci))...
    max(max(L56eventSpec.Wh.osci)) max(max(L6eventSpec.Wh.osci))]);

subplot(6,3,3);
imagesc(taxis,log10(wavespec.validfreq),L1eventSpec.Wh.osci');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmin cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L1 Osci');

subplot(6,3,6);
imagesc(taxis,log10(wavespec.validfreq),L23eventSpec.Wh.osci');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmin cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L2/3');

subplot(6,3,9);
imagesc(taxis,log10(wavespec.validfreq),L4eventSpec.Wh.osci');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmin cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L4');

subplot(6,3,12);
imagesc(taxis,log10(wavespec.validfreq),L5aeventSpec.Wh.osci');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmin cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L5a');

subplot(6,3,15);
imagesc(taxis,log10(wavespec.validfreq),L56eventSpec.Wh.osci');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmin cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L5b');

subplot(6,3,18);
imagesc(taxis,log10(wavespec.validfreq),L6eventSpec.Wh.osci');hold on;
colormap jet;
LogScale('y',10);
%LogScale('c',10);
caxis([cmin cmax]);
axis xy
ylim(log10([1 100]));
set(gca,'Ytick',log10([1 5 10 25 50 100]));
set(gca,'Yticklabel',{'1','5','10','25','50','100'});
%xlabel('time (ms)'); ylabel('f (Hz)');
plot([0 0],[log10(wavespec.validfreq(1)) log10(wavespec.validfreq(end))],'--k');hold on;
title('L6');

%NiceSave('Laminar_WhSpec',figfolder,baseName)