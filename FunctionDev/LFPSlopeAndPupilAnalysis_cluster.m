function [slopepupilcorr,PSShist,pupcyclePSS] = LFPSlopeAndPupilAnalysis_cluster( basePath,figfolder )
%UNTITLED3 Summary of this function goes here
%
%
%
%%
baseName = bz_BasenameFromBasepath(basePath);
sessionInfo = bz_getSessionInfo(basePath,'noPrompts',true);

%%
pupildilation = bz_LoadBehavior(basePath,'pupildiameter');

smoothwin =2; %s
nantimes = isnan(pupildilation.data);
pupildilation.data = pupildilation.data(~isnan(pupildilation.data));
pupildilation.dpdt = diff(smooth(pupildilation.data,smoothwin.*pupildilation.samplingRate,'moving')).*pupildilation.samplingRate;
pupildilation.dpdt = smooth(pupildilation.dpdt,smoothwin.*pupildilation.samplingRate,'moving');
pupildilation.timestamps = pupildilation.timestamps(~nantimes);

%%
slopepupilcorr.pup = zeros(size(sessionInfo.AnatGrps.Channels));
slopepupilcorr.dpdt = zeros(size(sessionInfo.AnatGrps.Channels));
for cc = 1:length(sessionInfo.AnatGrps.Channels)
    cc
    channum = sessionInfo.AnatGrps.Channels(cc);
    %channum = 31;
    %%
    lfp = bz_GetLFP(channum,'basepath',basePath,'noPrompts',true);
    
    %%
    dt = 0.5;
    winsize = 2;
    [specslope] = bz_PowerSpectrumSlope(lfp,winsize,dt,'showfig',false);
    
    %%
    specslope.pupilsize = interp1(pupildilation.timestamps,pupildilation.data,...
        specslope.timestamps,'nearest');
    specslope.dpdt = interp1(pupildilation.timestamps(1:end-1),pupildilation.dpdt,...
        specslope.timestamps,'nearest');
    
    %%
    [slopepupilcorr.pup(cc),slopepupilcorr.pup_p(cc)] =...
        corr(log10(specslope.pupilsize),specslope.data,'type','spearman','rows','complete');
    [ slopepupilcorr.dpdt(cc),slopepupilcorr.dpdt_p(cc)] =...
        corr(specslope.dpdt,specslope.data,'type','spearman','rows','complete');
    clear lfp
    
end

%%
figure;
subplot(2,2,1);
plot(1:sessionInfo.nChannels,slopepupilcorr.pup,'k')
hold on
plot(1:sessionInfo.nChannels,slopepupilcorr.dpdt,'k--')
xlabel('Channel (Sup-->Deep)');ylabel('Corr')
axis tight
legend('p','dpdt','location','southeast')

%% Take a look at the channels with dpdt and p
[~, bestchans.pup] = max(slopepupilcorr.pup);
[~, bestchans.dpdt] = max(slopepupilcorr.dpdt);
repchans = [sessionInfo.AnatGrps.Channels(bestchans.pup) sessionInfo.AnatGrps.Channels(bestchans.dpdt)];
%repchans = 42;
lfp = bz_GetLFP(sessionInfo.AnatGrps.Channels(repchans(1)),'basepath',basePath,'noPrompts',true);

dt = 0.2;
winsize = 2;
[specslope,specgram] = bz_PowerSpectrumSlope(lfp,winsize,dt,'showfig',true,...
    'saveMat',basePath);

specslope.pupilsize = interp1(pupildilation.timestamps,pupildilation.data,...
    specslope.timestamps,'nearest');
specslope.dpdt = interp1(pupildilation.timestamps(1:end-1),pupildilation.dpdt,...
    specslope.timestamps,'nearest');

%Check residuals

%%
numbins = 15;
bins = linspace(-0.5,0.5,numbins+1);
pupcyclePSS.bincenters = bins(1:end-1)+0.5.*diff(bins([1 2]));
bins([1 end])=[-Inf Inf];
[N,~,~,BINX,BINY] = histcounts2(log10(specslope.pupilsize),specslope.dpdt,...
    bins,bins);
pupcyclePSS.meanPSS = zeros(size(N));

for xx = 1:numbins
    for yy = 1:numbins
        pupcyclePSS.meanPSS(xx,yy) = nanmean(specslope.data(BINX==xx & BINY==yy));
    end
end

nbinthresh = 10;  %Must have more than 10 time windows
pupcyclePSS.meanPSS(N<nbinthresh) = nan;

%% PSS and UP/DOWN
SlowWaves = bz_LoadEvents(basePath,'SlowWaves');

%%
updown = {'DOWN','UP'};
UDcolor = {'b','r'};
for ss = 1:2
    SlowWaves.dur.(updown{ss}) = diff(SlowWaves.ints.(updown{ss}),1,2);
    SlowWaves.midpoint.(updown{ss}) = mean(SlowWaves.ints.(updown{ss}),2);
    SlowWaves.PSS.(updown{ss}) = interp1(specslope.timestamps,specslope.data,SlowWaves.midpoint.(updown{ss}));
end

%%
numbins = 30;
PSShist.bins = linspace(-2,0,numbins);
PSShist.hist = hist(specslope.data,PSShist.bins);

%%
exwinsize = 300;
exwin = bz_RandomWindowInIntervals(specslope.timestamps([1 end])',exwinsize);

winsize = 8;
maxtime = pupildilation.timestamps(...
    (pupildilation.timestamps>exwin(1) & pupildilation.timestamps<exwin(2)) & pupildilation.data==...
    max(pupildilation.data(pupildilation.timestamps>exwin(1) & pupildilation.timestamps<exwin(2))));
subsamplewin = maxtime(1)+[-0.5 0.5].*winsize;

mintime = pupildilation.timestamps(...
    (pupildilation.timestamps>exwin(1) & pupildilation.timestamps<exwin(2)) & pupildilation.data==...
    min(pupildilation.data(pupildilation.timestamps>exwin(1) & pupildilation.timestamps<exwin(2))));
subsamplewin2 = mintime(1)+[-0.5 0.5].*winsize;

figure
subplot(4,1,1)
imagesc(specgram.timestamps,log2(specgram.freqs),specgram.amp)
axis xy
xlim(exwin)
ylabel({'Specgram','f (Hz)'})
LogScale('y',2)

subplot(6,1,3)
plot(specslope.timestamps,specslope.data,'k')
axis tight
xlim(exwin)
box off
xlabel('t (s)')
ylabel('PSS')

subplot(6,1,4)
plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2)
hold on
plot(pupildilation.timestamps(1:end-1),pupildilation.dpdt,'k')
axis tight
xlim(exwin)
box off
ylabel('Pupil')
plot(exwin,[0 0],'k--')
plot(subsamplewin,3.*ones(size(subsamplewin)),'r','linewidth',2)
plot(subsamplewin2,3.*ones(size(subsamplewin2)),'r','linewidth',2)

subplot(6,2,11)
bz_MultiLFPPlot(lfp,'timewin',subsamplewin)
hold on
% plot(slow.timestamps,10000.*slow.data,'g','linewidth',1)
ylabel('High Pupil')
xlim(subsamplewin)

subplot(6,2,12)
bz_MultiLFPPlot(lfp,'timewin',subsamplewin2)
hold on
% plot(slow.timestamps,10000.*slow.data,'g','linewidth',1)
ylabel('Low Pupil')
xlim(subsamplewin2)

subplot(6,2,9)
plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2)
hold on
plot(pupildilation.timestamps(1:end-1),pupildilation.dpdt,'k')
plot(specslope.timestamps,specslope.data)
axis tight
xlim(subsamplewin)
box off
ylabel('Pupil')
plot(subsamplewin,[0 0],'k--')

subplot(6,2,10)
plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2)
hold on
plot(pupildilation.timestamps(1:end-1),pupildilation.dpdt,'k')
plot(specslope.timestamps,specslope.data)
axis tight
xlim(subsamplewin2)
box off
ylabel('Pupil')
plot(subsamplewin2,[0 0],'k--')

NiceSave('PSSexample',figfolder,baseName)

%%
figure
subplot(3,3,1)
h = imagesc(pupcyclePSS.bincenters,pupcyclePSS.bincenters,pupcyclePSS.meanPSS');
set(h,'AlphaData',~isnan(pupcyclePSS.meanPSS'));
hold on
plot(pupcyclePSS.bincenters([1 end]),[0 0],'k--')
LogScale('x',10)
axis xy
colorbar
xlabel('Pupil Area (med^-^1)')
ylabel('dp/dt')
LogScale('x',10)
title('PSS')

subplot(6,3,12)
plot(PSShist.bins,PSShist.hist,'k','linewidth',2)
box off
xlabel('PSS')
axis tight
xlim([-1.75 -0.25])

subplot(3,2,2)
for ss = 1:2
    plot(SlowWaves.PSS.(updown{ss}),log10(SlowWaves.dur.(updown{ss})),'.','color',UDcolor{ss},'markersize',3)
    hold on
end
xlabel('PSS');ylabel('Dur (s)')
axis tight
box off
LogScale('y',10)
legend(updown{:},'location','eastoutside')

subplot(3,3,4)
plot(log10(specslope.pupilsize),specslope.data,'k.','markersize',1)
xlabel('Pupil Area (med^-^1)');ylabel('PSS')
box off
axis tight
subplot(3,3,5)
plot((specslope.dpdt),specslope.data,'k.','markersize',1)
xlabel('dpdt (med^-^1s^-^1)');ylabel('PSS')
box off
axis tight

NiceSave('PSSandUPDOWN',figfolder,baseName)

%%
samplewin = bz_RandomWindowInIntervals(lfp.timestamps([1 end])',300);
winsize = 10;
maxtime = pupildilation.timestamps(...
    (pupildilation.timestamps>samplewin(1) & pupildilation.timestamps<samplewin(2)) & pupildilation.data==...
    max(pupildilation.data(pupildilation.timestamps>samplewin(1) & pupildilation.timestamps<samplewin(2))));
subsamplewin = maxtime(1)+[-0.5 0.5].*winsize;

mintime = pupildilation.timestamps(...
    (pupildilation.timestamps>samplewin(1) & pupildilation.timestamps<samplewin(2)) & pupildilation.data==...
    min(pupildilation.data(pupildilation.timestamps>samplewin(1) & pupildilation.timestamps<samplewin(2))));
subsamplewin2 = mintime(1)+[-0.5 0.5].*winsize;

figure
% subplot(3,3,1)
%     plot(slopepupilcorr.dpdt,'k','linewidth',0.5)
%     hold on
%     plot(slopepupilcorr.p,'k','linewidth',2)
%     plot(repchans(1),slopepupilcorr.dpdt(repchans(1)),'o')
%     plot(repchans(2),slopepupilcorr.p(repchans(2)),'o')
%     xlabel('Channel (sorted by depth)')
%     axis tight
%     box off

for cc = 1%:2
    subplot(6,6,12+cc)
    plot(log10(specslope.pupilsize),specslope.data(:,cc),'k.','markersize',3)
    xlabel('Pupil Area (med^-^1)');ylabel('Spectrum Slope')
    box off
    axis tight
    subplot(6,6,18+cc)
    plot((specslope.dpdt),specslope.data(:,cc),'k.','markersize',1)
    xlabel('dpdt (med^-^1s^-^1)');ylabel('Spectrum Slope')
    box off
    axis tight
end

% subplot(6,1,4)
%     bz_MultiLFPPlot(lfp,'timewin',samplewin)
subplot(6,6,9:12)
plot(specslope.timestamps,specslope.data)
axis tight
xlim(samplewin)
box off
xlabel('t (s)')
ylabel('PSS')

subplot(6,6,3:6)
plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2)
hold on
plot(pupildilation.timestamps(1:end-1),pupildilation.dpdt,'k')
axis tight
xlim(samplewin)
box off
ylabel('Pupil')
plot(samplewin,[0 0],'k--')
plot(subsamplewin,3.*ones(size(subsamplewin)),'r','linewidth',2)
plot(subsamplewin2,3.*ones(size(subsamplewin2)),'r','linewidth',2)

subplot(6,3,11:12)
bz_MultiLFPPlot(lfp,'timewin',samplewin)
hold on
% plot(slow.timestamps,10000.*slow.data,'g','linewidth',1)
ylabel('High Pupil')
xlim(subsamplewin)

subplot(6,3,17:18)
bz_MultiLFPPlot(lfp,'timewin',subsamplewin2)
hold on
% plot(slow.timestamps,10000.*slow.data,'g','linewidth',1)
ylabel('Low Pupil')
xlim(subsamplewin2)

subplot(6,3,8:9)
plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2)
hold on
plot(pupildilation.timestamps(1:end-1),pupildilation.dpdt,'k')
plot(specslope.timestamps,specslope.data)
axis tight
xlim(subsamplewin)
box off
ylabel('Pupil')
plot(subsamplewin,[0 0],'k--')

subplot(6,3,14:15)
plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2)
hold on
plot(pupildilation.timestamps(1:end-1),pupildilation.dpdt,'k')
plot(specslope.timestamps,specslope.data)
axis tight
xlim(subsamplewin2)
box off
ylabel('Pupil')
plot(subsamplewin2,[0 0],'k--')

NiceSave('PSSandPupil',figfolder,baseName)

end