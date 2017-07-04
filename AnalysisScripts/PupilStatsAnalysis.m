function [ puphist,pupACG ] = PupilStatsAnalysis( basePath,figfolder )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%
basePath= '/mnt/proraidDL/Database/WMProbeData/170421_Layers_LFP_Pupil_EMG_Emx1M1M3/170421_Layers_LFP_Pupil_EMG_170421_180427';
figfolder = '/mnt/data1/Dropbox/research/Current Projects/S1State/AnalysisScripts/figures/PupilStatsAnalysis';

baseName = bz_BasenameFromBasepath(basePath);
%%

bz_LoadBehavior('pupildiameter',basePath);

%% Compute the x-covariance

nantimes = isnan(pupildilation.data);
pupildilation.interpdata = interp1(pupildilation.timestamps(~nantimes),...
    pupildilation.data(~nantimes),pupildilation.timestamps);

acgwin = 200; %s
[pupACG.ACG,pupACG.tlag] = xcov(pupildilation.interpdata(~isnan(pupildilation.interpdata)),...
    round(acgwin.*pupildilation.samplingRate),'coeff');
pupACG.tlag = pupACG.tlag./pupildilation.samplingRate;

%% Histogram
puphist.bins = linspace(0,1,20);
puphist.counts = hist(pupildilation.interpdata,puphist.bins);
puphist.counts =puphist.counts./sum(puphist.counts);

%% D/Dt
smoothwin =2;%s
pupildt = diff(smooth(pupildilation.interpdata,smoothwin.*pupildilation.samplingRate,'moving')).*pupildilation.samplingRate;
pupildt = smooth(pupildt,smoothwin.*pupildilation.samplingRate,'moving');

%Histogram
pupdthist.bins={linspace(0,1,60),linspace(-0.15,0.15,60)};
pupdthist.counts = hist3([pupildilation.data(1:end-1),pupildt],pupdthist.bins);
pupdthist.counts = pupdthist.counts./sum(pupdthist.counts(:));

%% Figure
winsize = 300; %s
figure
subplot(4,1,1)
plot(pupildilation.timestamps,pupildilation.data,'k')
xlim(pupildilation.timestamps([1 end]))

subplot(4,1,2)
plot(pupildilation.timestamps,pupildilation.interpdata,'k','linewidth',2)
hold on
plot(get(gca,'xlim'),[0 0],'r-')
plot(pupildilation.timestamps(1:end-1),pupildt,'k','linewidth',1)
xlim(100+[0 winsize])
ylim([-0.15 1])

subplot(4,2,7)
plot(pupACG.tlag,pupACG.ACG,'k','linewidth',2)
hold on
axis tight
plot(get(gca,'xlim'),[0 0])
xlabel('t lag (s)')
ylabel('Pupil Autocov.')

subplot(4,2,5)
bar(puphist.bins,puphist.counts,'facecolor','k')
xlim([0 1])

subplot(2,2,4)
imagesc(pupdthist.bins{1},pupdthist.bins{2},pupdthist.counts')
hold on
plot(get(gca,'xlim'),[0 0],'r-')
colorbar

NiceSave('PupilStats',figfolder,baseName)

end

