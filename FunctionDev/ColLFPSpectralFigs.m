%%
basePath = pwd;
baseName = 'EM1M3';

figfolder = fullfile(basePath,'SummaryFigures');
if (~exist(figfolder,'dir'))
    mkdir(figfolder)
end

load(fullfile(basePath,['WT_EM1M3.ColLFPSpectralAnalysis.mat']));
WTLFP = groupColLFPSpectral;
WTDepth = groupDepth;

load(fullfile(basePath,['KO_EM1M3.ColLFPSpectralAnalysis.mat']));
KOLFP = groupColLFPSpectral;
KODepth = groupDepth;

%%
numberchans = 100;

[WTLFP.columnspec_all.db,normcolumn] = normCx(WTDepth.ndepth,WTLFP.columnspec_all.db,numberchans);
[WTLFP.columnspec_all.modz,~] = normCx(WTDepth.ndepth,WTLFP.columnspec_all.modz,numberchans);
[WTLFP.columnspec_all.frac,~] = normCx(WTDepth.ndepth,WTLFP.columnspec_all.frac,numberchans);
[WTLFP.columnspec_all.osci,~] = normCx(WTDepth.ndepth,WTLFP.columnspec_all.osci,numberchans);

[KOLFP.columnspec_all.db,~] = normCx(KODepth.ndepth,KOLFP.columnspec_all.db,numberchans);
[KOLFP.columnspec_all.modz,~] = normCx(KODepth.ndepth,KOLFP.columnspec_all.modz,numberchans);
[KOLFP.columnspec_all.frac,~] = normCx(KODepth.ndepth,KOLFP.columnspec_all.frac,numberchans);
[KOLFP.columnspec_all.osci,~] = normCx(KODepth.ndepth,KOLFP.columnspec_all.osci,numberchans);

[WTLFP.columnspec_sWh.db,normcolumn] = normCx(WTDepth.ndepth,WTLFP.columnspec_sWh.db,numberchans);
[WTLFP.columnspec_sWh.modz,~] = normCx(WTDepth.ndepth,WTLFP.columnspec_sWh.modz,numberchans);
[WTLFP.columnspec_sWh.frac,~] = normCx(WTDepth.ndepth,WTLFP.columnspec_sWh.frac,numberchans);
[WTLFP.columnspec_sWh.osci,~] = normCx(WTDepth.ndepth,WTLFP.columnspec_sWh.osci,numberchans);

[KOLFP.columnspec_sWh.db,~] = normCx(KODepth.ndepth,KOLFP.columnspec_sWh.db,numberchans);
[KOLFP.columnspec_sWh.modz,~] = normCx(KODepth.ndepth,KOLFP.columnspec_sWh.modz,numberchans);
[KOLFP.columnspec_sWh.frac,~] = normCx(KODepth.ndepth,KOLFP.columnspec_sWh.frac,numberchans);
[KOLFP.columnspec_sWh.osci,~] = normCx(KODepth.ndepth,KOLFP.columnspec_sWh.osci,numberchans);

[WTLFP.columnspec_lNWh.db,normcolumn] = normCx(WTDepth.ndepth,WTLFP.columnspec_lNWh.db,numberchans);
[WTLFP.columnspec_lNWh.modz,~] = normCx(WTDepth.ndepth,WTLFP.columnspec_lNWh.modz,numberchans);
[WTLFP.columnspec_lNWh.frac,~] = normCx(WTDepth.ndepth,WTLFP.columnspec_lNWh.frac,numberchans);
[WTLFP.columnspec_lNWh.osci,~] = normCx(WTDepth.ndepth,WTLFP.columnspec_lNWh.osci,numberchans);

[KOLFP.columnspec_lNWh.db,~] = normCx(KODepth.ndepth,KOLFP.columnspec_lNWh.db,numberchans);
[KOLFP.columnspec_lNWh.modz,~] = normCx(KODepth.ndepth,KOLFP.columnspec_lNWh.modz,numberchans);
[KOLFP.columnspec_lNWh.frac,~] = normCx(KODepth.ndepth,KOLFP.columnspec_lNWh.frac,numberchans);
[KOLFP.columnspec_lNWh.osci,~] = normCx(KODepth.ndepth,KOLFP.columnspec_lNWh.osci,numberchans);

[WTLFP.columnspec_lWh.db,normcolumn] = normCx(WTDepth.ndepth,WTLFP.columnspec_lWh.db,numberchans);
[WTLFP.columnspec_lWh.modz,~] = normCx(WTDepth.ndepth,WTLFP.columnspec_lWh.modz,numberchans);
[WTLFP.columnspec_lWh.frac,~] = normCx(WTDepth.ndepth,WTLFP.columnspec_lWh.frac,numberchans);
[WTLFP.columnspec_lWh.osci,~] = normCx(WTDepth.ndepth,WTLFP.columnspec_lWh.osci,numberchans);

[KOLFP.columnspec_lWh.db,~] = normCx(KODepth.ndepth,KOLFP.columnspec_lWh.db,numberchans);
[KOLFP.columnspec_lWh.modz,~] = normCx(KODepth.ndepth,KOLFP.columnspec_lWh.modz,numberchans);
[KOLFP.columnspec_lWh.frac,~] = normCx(KODepth.ndepth,KOLFP.columnspec_lWh.frac,numberchans);
[KOLFP.columnspec_lWh.osci,~] = normCx(KODepth.ndepth,KOLFP.columnspec_lWh.osci,numberchans);

%% FIGURE:
cmin = min([min(min(WTLFP.columnspec_all.db))...
    min(min(KOLFP.columnspec_all.db))]);
cmax = max([max(max(WTLFP.columnspec_all.db))...
    max(max(KOLFP.columnspec_all.db))]);

figure;
subplot(2,3,1);
imagesc(log10(WTLFP.freqs(:,:,1)),normcolumn,...
    WTLFP.columnspec_all.db');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);

subplot(2,3,4);
imagesc(log10(KOLFP.freqs(:,:,1)),normcolumn,...
    KOLFP.columnspec_all.db');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);

cmin = min([min(min(WTLFP.columnspec_all.frac))...
    min(min(KOLFP.columnspec_all.frac))]);
cmax = max([max(max(WTLFP.columnspec_all.frac))...
    max(max(KOLFP.columnspec_all.frac))]);

subplot(2,3,2);
imagesc(log10(WTLFP.validfreqs(:,:,1)),normcolumn,...
    WTLFP.columnspec_all.frac');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);

subplot(2,3,5);
imagesc(log10(KOLFP.validfreqs(:,:,1)),normcolumn,...
    KOLFP.columnspec_all.frac');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);

cmin = min([min(min(WTLFP.columnspec_all.osci))...
    min(min(KOLFP.columnspec_all.osci))]);
cmax = max([max(max(WTLFP.columnspec_all.osci))...
    max(max(KOLFP.columnspec_all.osci))]);

subplot(2,3,3);
imagesc(log10(WTLFP.validfreqs(:,:,1)),normcolumn,...
    WTLFP.columnspec_all.osci');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);

subplot(2,3,6);
imagesc(log10(KOLFP.validfreqs(:,:,1)),normcolumn,...
    KOLFP.columnspec_all.osci');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);

NiceSave('ColumnarPspec_all',figfolder,baseName)

%% FIGURE:
cmin = min([min(min(WTLFP.columnspec_lNWh.db))...
    min(min(KOLFP.columnspec_lNWh.db))]);
cmax = max([max(max(WTLFP.columnspec_lNWh.db))...
    max(max(KOLFP.columnspec_lNWh.db))]);

figure;
subplot(2,3,1);
imagesc(log10(WTLFP.freqs(:,:,1)),normcolumn,...
    WTLFP.columnspec_lNWh.db');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin*1.5 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);

subplot(2,3,4);
imagesc(log10(KOLFP.freqs(:,:,1)),normcolumn,...
    KOLFP.columnspec_lNWh.db');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin*1.5 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);

cmin = min([min(min(WTLFP.columnspec_lNWh.frac))...
    min(min(KOLFP.columnspec_lNWh.frac))]);
cmax = max([max(max(WTLFP.columnspec_lNWh.frac))...
    max(max(KOLFP.columnspec_lNWh.frac))]);

subplot(2,3,2);
imagesc(log10(WTLFP.validfreqs(:,:,1)),normcolumn,...
    WTLFP.columnspec_lNWh.frac');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin*1.5 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);

subplot(2,3,5);
imagesc(log10(KOLFP.validfreqs(:,:,1)),normcolumn,...
    KOLFP.columnspec_lNWh.frac');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin*1.5 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);

cmin = min([min(min(WTLFP.columnspec_lNWh.osci))...
    min(min(KOLFP.columnspec_lNWh.osci))]);
cmax = max([max(max(WTLFP.columnspec_lNWh.osci))...
    max(max(KOLFP.columnspec_lNWh.osci))]);

subplot(2,3,3);
imagesc(log10(WTLFP.validfreqs(:,:,1)),normcolumn,...
    WTLFP.columnspec_lNWh.osci');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin/4 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);

subplot(2,3,6);
imagesc(log10(KOLFP.validfreqs(:,:,1)),normcolumn,...
    KOLFP.columnspec_lNWh.osci');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin/4 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);

NiceSave('ColumnarPspec_lNWh',figfolder,baseName)

%% FIGURE:
cmin = min([min(min(WTLFP.columnspec_lWh.db))...
    min(min(KOLFP.columnspec_lWh.db))]);
cmax = max([max(max(WTLFP.columnspec_lWh.db))...
    max(max(KOLFP.columnspec_lWh.db))]);

figure;
subplot(2,3,1);
imagesc(log10(WTLFP.freqs(:,:,1)),normcolumn,...
    WTLFP.columnspec_lWh.db');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin*1.5 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);

subplot(2,3,4);
imagesc(log10(KOLFP.freqs(:,:,1)),normcolumn,...
    KOLFP.columnspec_lWh.db');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin*1.5 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);

cmin = min([min(min(WTLFP.columnspec_lWh.frac))...
    min(min(KOLFP.columnspec_lWh.frac))]);
cmax = max([max(max(WTLFP.columnspec_lWh.frac))...
    max(max(KOLFP.columnspec_lWh.frac))]);

subplot(2,3,2);
imagesc(log10(WTLFP.validfreqs(:,:,1)),normcolumn,...
    WTLFP.columnspec_lWh.frac');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin*1.5 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);

subplot(2,3,5);
imagesc(log10(KOLFP.validfreqs(:,:,1)),normcolumn,...
    KOLFP.columnspec_lWh.frac');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin*1.5 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);

cmin = min([min(min(WTLFP.columnspec_lWh.osci))...
    min(min(KOLFP.columnspec_lWh.osci))]);
cmax = max([max(max(WTLFP.columnspec_lWh.osci))...
    max(max(KOLFP.columnspec_lWh.osci))]);

subplot(2,3,3);
imagesc(log10(WTLFP.validfreqs(:,:,1)),normcolumn,...
    WTLFP.columnspec_lWh.osci');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin/4 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);

subplot(2,3,6);
imagesc(log10(KOLFP.validfreqs(:,:,1)),normcolumn,...
    KOLFP.columnspec_lWh.osci');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin/4 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);

NiceSave('ColumnarPspec_lWh',figfolder,baseName)

%% FIGURE:
cmin = min([min(min(WTLFP.columnspec_sWh.db))...
    min(min(KOLFP.columnspec_sWh.db))]);
cmax = max([max(max(WTLFP.columnspec_sWh.db))...
    max(max(KOLFP.columnspec_sWh.db))]);

figure;
subplot(2,3,1);
imagesc(log10(WTLFP.freqs(:,:,1)),normcolumn,...
    WTLFP.columnspec_sWh.db');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin*1.5 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);

subplot(2,3,4);
imagesc(log10(KOLFP.freqs(:,:,1)),normcolumn,...
    KOLFP.columnspec_sWh.db');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin*1.5 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);

cmin = min([min(min(WTLFP.columnspec_sWh.frac))...
    min(min(KOLFP.columnspec_sWh.frac))]);
cmax = max([max(max(WTLFP.columnspec_sWh.frac))...
    max(max(KOLFP.columnspec_sWh.frac))]);

subplot(2,3,2);
imagesc(log10(WTLFP.validfreqs(:,:,1)),normcolumn,...
    WTLFP.columnspec_sWh.frac');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin*1.5 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);

subplot(2,3,5);
imagesc(log10(KOLFP.validfreqs(:,:,1)),normcolumn,...
    KOLFP.columnspec_sWh.frac');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin*1.5 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);

cmin = min([min(min(WTLFP.columnspec_sWh.osci))...
    min(min(KOLFP.columnspec_sWh.osci))]);
cmax = max([max(max(WTLFP.columnspec_sWh.osci))...
    max(max(KOLFP.columnspec_sWh.osci))]);

subplot(2,3,3);
imagesc(log10(WTLFP.validfreqs(:,:,1)),normcolumn,...
    WTLFP.columnspec_sWh.osci');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin/4 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);

subplot(2,3,6);
imagesc(log10(KOLFP.validfreqs(:,:,1)),normcolumn,...
    KOLFP.columnspec_sWh.osci');
axis tight
LogScale('x',10)
LogScale('c',10)
colormap('jet');
caxis([cmin/4 cmax]);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
xtickangle(45);
xlabel('f (Hz)');
set(gca,'Ytick',[0.1 0.35 0.5 0.6 0.9]);
set(gca,'Yticklabel',{'L1/2','L3/4','L4/5a','L5b','L6'});
set(gca,'YGrid','on', 'GridColor','w','GridAlpha',0.45);

NiceSave('ColumnarPspec_sWh',figfolder,baseName)

%% FIGURE:
figure;
shadederror(log10(WTLFP.freqs(:,:,1)),...
    nanmean(WTLFP.columnspec_lNWh.db,2),...
    nanstd(WTLFP.columnspec_lNWh.db,0,2),'k','none',1,'k');
hold on;
shadederror(log10(KOLFP.freqs(:,:,1)),...
    nanmean(WTLFP.columnspec_lWh.db,2),...
    nanstd(WTLFP.columnspec_lWh.db,0,2),'r','none',1,'r');
xlabel('f (hz)');
ylabel('power (dB)')
axis square
LogScale('x',10);
xlim(log10([1 100]));
set(gca,'Xtick',log10([1 5 10 25 50 100]));
set(gca,'Xticklabel',{'1','5','10','25','50','100'});
xtickangle(45);
xlabel('f (Hz)');

NiceSave('Pspec_Wh_NWh',figfolder,baseName)
