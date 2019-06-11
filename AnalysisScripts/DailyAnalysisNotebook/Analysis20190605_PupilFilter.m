function [ ] = AnalysisXXXXXXXX(basePath,figfolder)
% Date XX/XX/20XX
%
%Question: 
%
%Plots
%-
%-
%
%% Load Header
%Initiate Paths
reporoot = '/home/dlevenstein/ProjectRepos/ACh-and-CorticalState/';
%reporoot = '/Users/dlevenstein/Project Repos/ACh-and-CorticalState/';
basePath = '/mnt/proraidDL/Database/WMData/AChPupil/171209_WT_EM1M3/';
%basePath = '/mnt/proraidDL/Database/WMData/AChPupil/180706_WT_EM1M3/';
%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath,'noPrompts',true);



%% Loading behavior...
% Pupil diameter
pupildilation = bz_LoadBehavior(basePath,'pupildiameter');

smoothwin = 0.5; %s
pupildilation.data = smooth(pupildilation.data,smoothwin.*pupildilation.samplingRate,'moving');
nantimes = isnan(pupildilation.data);
pupildilation.data = pupildilation.data(~isnan(pupildilation.data));

if length(pupildilation.data) < 1
    warning('Not enough pupil data >)');
    return
end

smoothwin = 2; %s
pupildilation.dpdt = diff(smooth(pupildilation.data,smoothwin.*pupildilation.samplingRate,'moving')).*pupildilation.samplingRate;
pupildilation.dpdt = smooth(pupildilation.dpdt,smoothwin.*pupildilation.samplingRate,'moving');
pupildilation.timestamps = pupildilation.timestamps(~nantimes);



% EMG
EMGwhisk = bz_LoadBehavior(basePath,'EMGwhisk');

EMGwhisk.iswhisk = InIntervals(EMGwhisk.timestamps,EMGwhisk.ints.Wh);
%%
maxtimejump =1;
%Restricting SPONT UP/DOWNs
load(fullfile(basePath,[baseName,'.MergePoints.events.mat']),'MergePoints');
sidx = find(startsWith(MergePoints.foldernames,"Spont"));
sponttimes = [MergePoints.timestamps(sidx(1),1) MergePoints.timestamps(sidx(end),2)];

spontidx = find(EMGwhisk.timestamps < sponttimes(2));
EMGwhisk.timestamps = EMGwhisk.timestamps(spontidx);
EMGwhisk.EMGsm = EMGwhisk.EMGsm(spontidx);
EMGwhisk.EMGsm = NanPadJumps( EMGwhisk.timestamps,EMGwhisk.EMGsm,maxtimejump );


spontidx = find(pupildilation.timestamps < sponttimes(2));
pupildilation.timestamps = pupildilation.timestamps(spontidx);
pupildilation.data = pupildilation.data(spontidx);

pupildilation.EMG = interp1(EMGwhisk.timestamps,EMGwhisk.EMGsm,pupildilation.timestamps,'nearest');


%% Measure coupling with a range of lower/uppwer bounds, at a few different ncyc resolutions

lowerbounds = logspace(-2,-1,20);
upperbounds = logspace(-1.5,0,20);
orders = 1:4;

%For each combination of bands

%Filter the pupil
amp2 = pupildilation.EMG./nanmean(pupildilation.EMG);
coupling = nan(length(lowerbounds),length(upperbounds),length(orders));
for oo = orders
for ll = 1:length(lowerbounds)
    ll
    for uu = 1:length(upperbounds)
        uu
        if upperbounds(uu)<=lowerbounds(ll)
            coupling(ll,uu,oo)=nan;
            continue
        end
        pupilcycle = bz_Filter(pupildilation,'passband',[lowerbounds(ll) upperbounds(uu)],...
            'filter' ,'fir1','order',oo);

    %Normalize power to the median (mean?)
        amp1 = pupilcycle.amp./nanmean(pupilcycle.amp);
        
        coupling(ll,uu,oo) = abs(mean(amp1.*amp2.*exp(1i.*pupilcycle.phase)));
    %Calculate the coupling with EMG
    end
end
end

%%
trybounds = [0.02 0.33];
figure
for oo = orders
    subplot(2,2,oo)
imagesc(log10(lowerbounds),log10(upperbounds),coupling(:,:,oo)')
hold on
axis xy
plot(log10(trybounds(1)),log10(trybounds(2)),'r+')
LogScale('xy',10)
clim([0.25 0.75])
xlabel('Lower Bound (Hz)');ylabel('Upper Bounr (Hz)')
colorbar
title(num2str(oo))
end


%% Look at traces for each order...
pupilcycle = bz_Filter(pupildilation,'passband',trybounds,...
            'filter' ,'fir1','order',1);


% [thresh,cross,bihist,diptest,overthresh] = bz_BimodalThresh(log10(pupilcycle.amp),...
%     'Schmidt',true,'setthresh',pupilcycle.pupthresh);
        
maxtimejump = 1; %s
pupilcycle.amp = NanPadJumps( pupilcycle.timestamps,pupilcycle.amp,maxtimejump );
pupilcycle.phase = NanPadJumps( pupilcycle.timestamps,pupilcycle.phase,maxtimejump );

pupilcycle.pupthresh = -0.8;
pupilcycle.highpup = log10(pupilcycle.amp)>pupilcycle.pupthresh; 
pupilcycle.lowpup = log10(pupilcycle.amp)<=pupilcycle.pupthresh;

%%
EMGwhisk.pupphase = interp1(pupilcycle.timestamps,pupilcycle.phase,EMGwhisk.timestamps,'nearest');
EMGwhisk.pupamp = interp1(pupilcycle.timestamps,pupilcycle.amp,EMGwhisk.timestamps,'nearest');

[pupilphaseEMG.meanZ,pupilphaseEMG.N,pupilphaseEMG.Xbins,pupilphaseEMG.Ybins] = ...
    ConditionalHist3( EMGwhisk.pupphase(EMGwhisk.EMGsm~=0),...
    log10(EMGwhisk.pupamp(EMGwhisk.EMGsm~=0)),log10(EMGwhisk.EMGsm(EMGwhisk.EMGsm~=0)),...
    'minXY',500,'Xbounds',[-pi pi],'Ybounds',[-2 0],...
    'numXbins',30,'numYbins',30);


%%
figure
cosx = linspace(-pi,pi,100);
cospamp = [0.08 0.8];
figure


subplot(3,2,1)
a = imagesc(pupilphaseEMG.Xbins,pupilphaseEMG.Ybins,pupilphaseEMG.meanZ');
hold on
alpha(a,double(~isnan(pupilphaseEMG.meanZ')))

%imagesc(pupilphaseEMG.Xbins+2*pi,pupilphaseEMG.Ybins,pupilphaseEMG.meanZ')
crameri lapaz
plot([-pi 3*pi],pupilcycle.pupthresh.*[1 1],'w--')
plot(cosx,(cos(cosx)+1).*cospamp(pp)-2,'k')
axis xy
box off
%xlim([-pi 3*pi])
ColorbarWithAxis([-0.7 0.7],'Mean EMG')
LogScale('c',10)
LogScale('y',10)
xlabel('Pupil Phase');ylabel('Pupil Amplitude')
 bz_piTickLabel('x')
 ylim([-2 0.1])

%% Example Figure
windows(1,:) = [100 250];
windows(2,:) = bz_RandomWindowInIntervals(pupildilation.timestamps([1 end]),150);
windows(3,:) = bz_RandomWindowInIntervals(pupildilation.timestamps([1 end]),150);

figure;
for ww = 1:3
subplot(3,1,ww);
%plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2); hold on;
scatter(pupildilation.timestamps,pupildilation.data,3,pupilcycle.phase,'filled')
hold on
plot(pupilcycle.timestamps,pupilcycle.amp.*2)
%scatter(pupilcycle.timestamps,pupilcycle.amp,3,pupilcycle.phase,'filled')
plot(EMGwhisk.timestamps,EMGwhisk.EMG./max(EMGwhisk.EMG),...
    'color',[0.5 0.5 0.5],'linewidth',0.5);
plot(EMGwhisk.timestamps,5*EMGwhisk.EMGsm./(max(EMGwhisk.EMG)),...
    'color','k','linewidth',0.5);
plot(pupilcycle.timestamps(pupilcycle.highpup),3*ones(sum(pupilcycle.highpup),1),'r.')
plot(pupilcycle.timestamps(~pupilcycle.highpup),3*ones(sum(~pupilcycle.highpup),1),'k.')
plot(EMGwhisk.timestamps(EMGwhisk.iswhisk),2.8*ones(sum(EMGwhisk.iswhisk),1),'b.')
%plot(highpupildata.timestamps,highpupildata.data+nanmean(pupildilation.data),'r')
colormap(gca,hsv)
ColorbarWithAxis([min(pupilcycle.phase) max(pupilcycle.phase)],['pupil phase'])
%h1 = plot(get(gca,'xlim'),[0 0],'k-');
xlim(windows(ww,:)); ylim([-0.05 3]);
bz_ScaleBar('s')
ylabel('Pupil');
%set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%legend({'diameter','phase','dPdt'},'location','northeast');

end