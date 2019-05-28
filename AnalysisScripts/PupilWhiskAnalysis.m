function [puphist,pupACG,pupPSD,pupdthist,...
    EMGhist,Whdurhist,...
    pupilEMGdist,pupildynamicsEMG,...
    pupilEMGcorr,pwCCG,...
    pupilphaseEMG,WPcoupling] = PupilWhiskAnalysis(basePath,figfolder)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% DEV
%basePath = '/mnt/proraidDL/Database/WMProbeData/170421_Layers_LFP_Pupil_EMG_Emx1M1M3/170421_Layers_LFP_Pupil_EMG_170421_180427/';
%basePath = '/mnt/proraidDL/Database/WMProbeData/170606_Layers_LFP_Pupil_EMG_Emx1M2M4/170606_Layers_LFP_Pupil_EMG_170606_204730/';
basePath= '/mnt/proraidDL/Database/WMProbeData/170421_Layers_LFP_Pupil_EMG_Emx1M1M3/170421_Layers_LFP_Pupil_EMG_170421_180427';
%basePath = '/home/dlevenstein/ProjectRepos/ACh-and-CorticalState/Dataset/180605_WT_M1M3_LFP_Layers_Pupil_EMG_180605_121846';
%basePath = 'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180706_WT_EM1M3';
%basePath = 'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180525_KO_EM1M3\Spont2';
%reporoot = 'H:\Personal\GitHub\ACh-and-CorticalState';
%figfolder = [reporoot,'/AnalysisScripts/AnalysisFigs/PupilWhiskAnalysis'];

%basePath = pwd;
%figfolder = '/home/dlevenstein/ProjectRepos/ACh-and-CorticalState/AnalysisScripts/AnalysisFigs/PupilWhiskAnalysis';



baseName = bz_BasenameFromBasepath(basePath);
%figfolder = '/mnt/data1/Dropbox/research/Current Projects/S1State/AnalysisScripts/figures/PupilWhiskAnalysis';

%% PUPIL %%
pupildilation = bz_LoadBehavior(basePath,'pupildiameter');

smoothwin =2;%s
pupildilation.dpdt = diff(smooth(pupildilation.data,smoothwin.*pupildilation.samplingRate,'moving')).*pupildilation.samplingRate;
pupildilation.dpdt = smooth(pupildilation.dpdt,smoothwin.*pupildilation.samplingRate,'moving');

nantimes = isnan(pupildilation.data);
pupildilation.interpdata = interp1(pupildilation.timestamps(~nantimes),...
    pupildilation.data(~nantimes),pupildilation.timestamps);

%% Compute the x-covariance


acgwin = 200; %s
[pupACG.ACG,pupACG.tlag] = xcov(pupildilation.interpdata(~isnan(pupildilation.interpdata)),...
    round(acgwin.*pupildilation.samplingRate),'coeff');
pupACG.tlag = pupACG.tlag./pupildilation.samplingRate;

%% Histogram
puphist.bins = linspace(0,5,20);
puphist.counts = hist(pupildilation.interpdata,puphist.bins);
puphist.counts =puphist.counts./sum(puphist.counts);

%% D/Dt

%Histogram
pupdthist.bins={linspace(-0.5,0.5,25),linspace(-0.5,0.5,25)};
pupdthist.counts = hist3([log10(pupildilation.interpdata(1:end-1)),pupildilation.dpdt],pupdthist.bins);
pupdthist.counts = pupdthist.counts./sum(pupdthist.counts(:));


%% Calculate FFT of the pupil
frange = [0.001 1];
pupilspec = bz_WaveSpec(pupildilation.interpdata(~isnan(pupildilation.interpdata)),...
    'frange',frange,'nfreqs',200,'ncyc',3,'samplingRate',pupildilation.samplingRate);
freqs = pupilspec.freqs;
spec = pupilspec.data;
spec = (abs(spec));
clear pupilspec
%%
lowfilter = [0.01 0.1];
highfilter = [0.3 0.8];

pupil4filter = pupildilation;
pupil4filter.data = pupildilation.interpdata(~isnan(pupildilation.interpdata));
%pupil4filter.t
pupil4filter.timestamps = pupil4filter.timestamps(~isnan(pupildilation.interpdata));
lowpupildata = bz_Filter(pupil4filter,'passband',lowfilter,'filter' ,'fir1','order',3);
%highpupildata = bz_Filter(pupil4filter,'passband',highfilter,'filter' ,'fir1');

pupPSD.freqs = freqs;
pupPSD.psd = mean(log10(spec),1);
clear spec
%% Figure: PUPIL
winsize = 300; %s
viewwin = bz_RandomWindowInIntervals(pupildilation.timestamps([1 end]),winsize);
figure
    subplot(6,4,1:3)
        plot(pupildilation.timestamps,pupildilation.data,'k')
        hold on
        xlim(pupildilation.timestamps([1 end]))
        ylabel('Pupil: All')
        box off
    subplot(6,1,2:3)
        plot(pupildilation.timestamps,pupildilation.interpdata,'k','linewidth',2)
        hold on
        box off
      %  scatter(lowpupildata.timestamps,lowpupildata.data+0.5.*nanmean(pupildilation.data),4,lowpupildata.phase)
      %  plot(highpupildata.timestamps,highpupildata.data+nanmean(pupildilation.data),'r')
       % colormap(gca,hsv)
        colorbar
        plot(get(gca,'xlim'),[0 0],'r-')
        plot(pupildilation.timestamps(1:end-1),pupildilation.dpdt,'k','linewidth',1)
        xlim(viewwin)
        ylim([-1 4])
        ylabel('Pupil: Zoom')
        xlabel('t (s)')
        
    subplot(6,4,4)
        bar(puphist.bins,puphist.counts,'facecolor','k')
        xlim([0 5])
        box off
        axis tight
        xlabel('Pupil Diameter');title('Pupil Diameter Distribution')
        
    subplot(4,2,5)
        plot(pupACG.tlag,pupACG.ACG,'k','linewidth',2)
        hold on
        axis tight
        box off
        plot(get(gca,'xlim'),[0 0])
        xlabel('t lag (s)')
        ylabel('Pupil Autocov.')

    subplot(2,2,4)
        imagesc(pupdthist.bins{1},pupdthist.bins{2},pupdthist.counts')
            colormap(gca,[1 1 1;colormap])

        hold on
        plot(get(gca,'xlim'),[0 0],'r-')
        colorbar
        xlabel('Pupil Diameter');ylabel('d/dt')
        title('Pupil Dynamics Histogram')
        LogScale('x',10)
        %caxis([-0.01 max(pupdthist.counts(:))])

        
        
% subplot(4,2,7)
% plot(log10(freqs),mean(spec,2),'k','linewidth',2)
% hold on
% plot(log10(lowfilter),[0 0 ],'r')
% %plot(log10(highfilter),[0 0 ],'r')
% LogScale('x',10)
% xlabel('f (Hz)');ylabel('power')
% axis tight
 
    subplot(4,2,7)
        plot(log10(pupPSD.freqs),pupPSD.psd,'k','linewidth',2)
        hold on
        plot(log10(lowfilter),[-1 -1],'r')
        %plot(log10(highfilter),[0 0 ],'r')
        LogScale('x',10)
        xlabel('f (Hz)');ylabel('power (dB)')
        axis tight
 
% subplot(2,2,3)
% plot(log10(freqs),(mean(spec,2)./(1./freqs)'))
% hold on
% plot(log10(lowfilter),[0 0 ],'r')
% plot(log10(highfilter),[0 0 ],'r')
% LogScale('x',10)
% xlabel('f (Hz)')
% axis tight
        
        
NiceSave('PupilStats',figfolder,baseName)


%% Pupil Space figure 
viewwin = bz_RandomWindowInIntervals(pupildilation.timestamps([1 end]),winsize);


figure
    subplot(6,1,1:2)
        plot(pupildilation.timestamps,pupildilation.interpdata,'k','linewidth',2)
        hold on
        box off
        scatter(lowpupildata.timestamps,lowpupildata.data,4,lowpupildata.phase)
      %  plot(highpupildata.timestamps,highpupildata.data+nanmean(pupildilation.data),'r')
        colormap(gca,hsv)
        colorbar
        %plot(get(gca,'xlim'),[0 0],'r-')
        xlim(viewwin)
        ylim([-1 4])
        ylabel('Pupil: Zoom')
        xlabel('t (s)')
    subplot(2,3,4)
        colormap(gca,hsv)
        scatter(log10(pupildilation.data(1:end-1)),pupildilation.dpdt,0.2,lowpupildata.phase(1:end-1))
        %plot(pupildilation.data(1:end-1),pupildilation.dpdt,'k.','markersize',2)
        ylim([-0.7 0.7]);xlim([-0.5 0.5])
        LogScale('x',10)
        %
        %colorbar
        
        ColorbarWithAxis([-pi pi],'Pupil Phase','location','northoutside')
        xlabel('Pupil Diameter');ylabel('d/dt')
        
    subplot(2,3,5)
        %colormap(gca,distcolor)
        scatter(log10(pupildilation.data(1:end-1)),pupildilation.dpdt,0.2,lowpupildata.amp(1:end-1))
        %plot(pupildilation.data(1:end-1),pupildilation.dpdt,'k.','markersize',2)
        ylim([-0.7 0.7]);xlim([-0.5 0.5])
        LogScale('x',10)
        %
        %colorbar
        
        ColorbarWithAxis([0 1.2],'Pupil Power','location','northoutside')
        xlabel('Pupil Diameter');ylabel('d/dt')
        
%     subplot(2,2,3)
%         hist3([lowpupildata.phase,log10(lowpupildata.amp)])
%         
        
NiceSave('PupilSpace',figfolder,baseName)
%% EMG %%
EMGwhisk = bz_LoadBehavior( basePath,'EMGwhisk');
if isempty(EMGwhisk)
    EMGwhisk = GetWhiskFromEMG(basePath);
end

%%
EMGwhisk.pupiltime = interp1(EMGwhisk.timestamps,EMGwhisk.EMGsm,pupildilation.timestamps);

%% Histogram
EMGhist.bins = linspace(0,20,20);
EMGhist.counts = hist(EMGwhisk.EMGsm,EMGhist.bins);
EMGhist.counts =EMGhist.counts./sum(EMGhist.counts);

EMGhist.logbins = linspace(-1.5,1.5,20);
EMGhist.logcounts = hist(log10(EMGwhisk.EMGsm),EMGhist.logbins);
EMGhist.logcounts =EMGhist.logcounts./sum(EMGhist.logcounts);

EMGhist.threshold = EMGwhisk.detectorparms.Whthreshold;
%% Whisk Durations
EMGwhisk.Whdurs = diff(EMGwhisk.ints.Wh,1,2);
EMGwhisk.InterWhdurs = EMGwhisk.ints.Wh(2:end,1)-EMGwhisk.ints.Wh(1:end-1,2);

Whdurhist.bins = linspace(-1,2,20);
Whdurhist.Whdurs = hist(log10(EMGwhisk.Whdurs),Whdurhist.bins);
Whdurhist.Whdurs = Whdurhist.Whdurs./diff(EMGwhisk.timestamps([1 end])); %Units: whisks/s
Whdurhist.InterWhdurs = hist(log10(EMGwhisk.InterWhdurs),Whdurhist.bins);
Whdurhist.InterWhdurs = Whdurhist.InterWhdurs./diff(EMGwhisk.timestamps([1 end])); %Units: whisks/s



%% Figure: EMG/WHISK
figure
    subplot(4,1,1)
        plot(EMGwhisk.timestamps,EMGwhisk.EMG,'color',[0.5 0.5 0.5],'linewidth',0.5)
        hold on
        plot(pupildilation.timestamps,EMGwhisk.pupiltime,'b','linewidth',2)
        plot(EMGwhisk.ints.Wh',...
            ones(size(EMGwhisk.ints.Wh))',...
            'g-','linewidth',1)
%             EMGwhisk.detectorparms.Whthreshold.*ones(size(EMGwhisk.ints.Wh))',...
%             'g-','linewidth',1)
        %plot(EMGwhisk.timestamps,EMGwhisk.EMGsm,'b')
        ylabel('EMG: All')
        axis tight
        xlabel('t (s)')
        box off
        xlim(pupildilation.timestamps([1 end]))
        ylim([-20 40])
    subplot(4,1,2)
        plot(EMGwhisk.timestamps,EMGwhisk.EMG,'color',[0.5 0.5 0.5],'linewidth',0.5)
        hold on
        plot(pupildilation.timestamps,EMGwhisk.pupiltime,'b','linewidth',2)
        plot(EMGwhisk.ints.Wh',...
            ones(size(EMGwhisk.ints.Wh))',...
            'g-','linewidth',1)
%             EMGwhisk.detectorparms.Whthreshold.*ones(size(EMGwhisk.ints.Wh))',...
%             'g-','linewidth',1)
        %plot(EMGwhisk.timestamps,EMGwhisk.EMGsm,'b')
        axis tight
        ylabel('EMG: zoom')
        xlabel('t (s)')
        box off
        xlim(100+[0 winsize]);ylim([-20 40])
    subplot(4,2,5)
        bar(EMGhist.bins,EMGhist.counts,'facecolor','b')
        axis tight
        xlabel('EMG');title('EMG Distribution(linear)')
        
    subplot(4,2,7)
        bar(EMGhist.logbins,EMGhist.logcounts,'facecolor','b')
        axis tight
        hold on
        %plot(log10(EMGwhisk.detectorparms.NWhthreshold).*[1 1],get(gca,'ylim'),'r--')
        plot(log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],get(gca,'ylim'),'g--')
        %xlim([-1 1.5])
        LogScale('x',10)
        xlabel('EMG');title('EMG Distribution(log)')
        
    subplot(4,2,8)
        plot(Whdurhist.bins,Whdurhist.Whdurs,'g','linewidth',2)
        hold on
        plot(Whdurhist.bins,Whdurhist.InterWhdurs,'r','linewidth',2)
        LogScale('x',10)
        legend('Wh','NWh')
        xlabel('Duration');title('Whisking/nonwhisking durations')
NiceSave('EMGStats',figfolder,baseName)

%% Relate EMG/Whisk with Pupil %%

%Codistribution
pupilEMGdist.bins = {linspace(-0.5,0.5,70),linspace(-1.5,1,70)};
[pupilEMGdist.counts,pupilEMGdist.bins] = hist3([log10(pupildilation.data),log10(EMGwhisk.pupiltime)],pupilEMGdist.bins);
pupilEMGdist.counts = pupilEMGdist.counts./sum(pupilEMGdist.counts(:));

%XCovariance
[pupilEMGcorr.xcorr,pupilEMGcorr.corrlags] = xcov(pupildilation.interpdata(~isnan(pupildilation.interpdata)),EMGwhisk.pupiltime(~isnan(pupildilation.interpdata)),'unbiased');
pupilEMGcorr.corrlags = pupilEMGcorr.corrlags.*(1./pupildilation.samplingRate);

%Pupil/EMG at whisking onset
[pwCCG.pupil.WhOn,t_lag,~,alltrans.pupil.WhOn,skippedWh] = EventVsContinousCCG(pupildilation.data,pupildilation.timestamps,EMGwhisk.ints.Wh(:,1),20);
[pwCCG.pupilxy.WhOn(:,1),t_lag,~,~] = EventVsContinousCCG(pupildilation.pupilxy(:,1),pupildilation.timestamps,EMGwhisk.ints.Wh(:,1),20);
[pwCCG.pupilxy.WhOn(:,2),t_lag,~,~] = EventVsContinousCCG(pupildilation.pupilxy(:,2),pupildilation.timestamps,EMGwhisk.ints.Wh(:,1),20);
[pwCCG.EMG.WhOn,t_lag,~,alltrans.EMG.WhOn] = EventVsContinousCCG(EMGwhisk.pupiltime,pupildilation.timestamps,EMGwhisk.ints.Wh(:,1),20);
pwCCG.t_lag = t_lag;
%% Histogram: Whisking by Pupil Dynamics

[pupildynamicsEMG,pupildynamicsEMG.bins]=PairMatHist(log10(EMGwhisk.pupiltime(1:end-1)),...
    [log10(pupildilation.data(1:end-1)),pupildilation.dpdt],15,[-0.5 0.5]);

[~,pupildilation.iswhisk] = RestrictInts(pupildilation.timestamps,EMGwhisk.ints.Wh);
[pWhisk]=PairMatHist(single(pupildilation.iswhisk(1:end-1)),[log10(pupildilation.data(1:end-1)),pupildilation.dpdt],pupildynamicsEMG.binedges);
pupildynamicsEMG.pWhisk = pWhisk.mean;


[pupilphaseEMG,pupilphaseEMG.bins]=PairMatHist(log10(EMGwhisk.pupiltime),...
    [log10(lowpupildata.amp),lowpupildata.phase],20,[-pi pi]);
[pWhisk]=PairMatHist(single(pupildilation.iswhisk),[log10(lowpupildata.amp),lowpupildata.phase],pupilphaseEMG.binedges);
pupilphaseEMG.pWhisk = pWhisk.mean;


%Wh onsets/offsets in pupil space
whints_pupil = interp1(pupildilation.timestamps,pupildilation.data,EMGwhisk.ints.Wh);
whints_pupildt = interp1(pupildilation.timestamps(1:end-1),pupildilation.dpdt,EMGwhisk.ints.Wh);
whints_pupilphase = interp1(lowpupildata.timestamps,lowpupildata.phase,EMGwhisk.ints.Wh);
whints_pupilamp = interp1(lowpupildata.timestamps,lowpupildata.amp,EMGwhisk.ints.Wh);

[pWhiskStart]=PairMatHist(EMGwhisk.Whdurs(:,1),[log10(whints_pupil(:,1)),whints_pupildt(:,1)],pupildynamicsEMG.binedges);
pupildynamicsEMG.meanWhdur = pWhiskStart.mean;
pupildynamicsEMG.numWhstarts = pWhiskStart.num;
pupildynamicsEMG.pWhstarts = pupildynamicsEMG.numWhstarts./sum(pupildynamicsEMG.numWhstarts(:));
pupildynamicsEMG.occupancy = pupildynamicsEMG.num./pupildilation.samplingRate;
pupildynamicsEMG.Whstartrate = pupildynamicsEMG.numWhstarts./pupildynamicsEMG.occupancy;


[pWhiskStart]=PairMatHist(EMGwhisk.Whdurs(:,1),[log10(whints_pupilamp(:,1)),whints_pupilphase(:,1)],pupilphaseEMG.binedges);
pupilphaseEMG.meanWhdur = pWhiskStart.mean;
pupilphaseEMG.numWhstarts = pWhiskStart.num;
pupilphaseEMG.pWhstarts = pupilphaseEMG.numWhstarts./sum(pupilphaseEMG.numWhstarts(:));
pupilphaseEMG.occupancy = pupilphaseEMG.num./pupildilation.samplingRate;
pupilphaseEMG.Whstartrate = pupilphaseEMG.numWhstarts./pupilphaseEMG.occupancy;


occupancythresh = 2; %s (don't count bins with less than 1s total time)
whthresh = 2;
pupildynamicsEMG.mean(pupildynamicsEMG.occupancy<=occupancythresh)=nan;
pupilphaseEMG.mean(pupilphaseEMG.occupancy<=occupancythresh)=nan;
pupildynamicsEMG.Whstartrate(pupildynamicsEMG.occupancy<=occupancythresh)=nan;
pupilphaseEMG.Whstartrate(pupilphaseEMG.occupancy<=occupancythresh)=nan;

pupildynamicsEMG.meanWhdur(pupildynamicsEMG.numWhstarts<=whthresh)=nan;
pupilphaseEMG.meanWhdur(pupilphaseEMG.numWhstarts<=whthresh)=nan;
%%
figure
subplot(2,2,1)
scatter(log10(whints_pupil(:,1)),whints_pupildt(:,1),1,log10(EMGwhisk.Whdurs))
xlabel('Pupil Diameter');ylabel('d/dt')
LogScale('x',10)

subplot(2,2,2)
scatter(whints_pupilphase(:,1),log10(whints_pupilamp(:,1)),1,log10(EMGwhisk.Whdurs))
hold on
scatter(whints_pupilphase(:,1)+2*pi,log10(whints_pupilamp(:,1)),1,log10(EMGwhisk.Whdurs))
xlabel('Pupil Phase');ylabel('Pupil Amplitude')
LogScale('y',10)

NiceSave('PupilandWhiskDur',figfolder,baseName)
%%

PhaseAmpCouplingByAmp( lowpupildata.phase,log10(lowpupildata.amp),...
    log10(EMGwhisk.pupiltime),10 );

%% Figure WHisk in pupil space
emgcolor = [1 1 1;makeColorMap([0.5 0.5 0.5],[0 0 0.8])];


figure

    subplot(4,2,1)
        
        imagesc(pupildynamicsEMG.bins,pupildynamicsEMG.bins,(pupildynamicsEMG.mean)');
        colormap(gca,emgcolor)
        axis xy
        xlim([-0.5 0.5]);ylim([-0.5 0.5])
        %caxis([-0.01 1])
        %colormap(gca)
        colorbar
        %ColorbarWithAxis([-0.1 1],'EMG Envelope')
        xlabel('Pupil Diameter');ylabel('d/dt')
        title('EMG')
        hold on
        caxis([-0.75 0.5])
        LogScale('c',10)
        %colormap(x,emgcolor)
        %colormap(y,'jet')
        plot(log10(whints_pupil(:,1)),whints_pupildt(:,1),'g.','markersize',0.1)
        %plot(whints_pupil(:,2),whints_pupildt(:,2),'r.')
        plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color',0.5.*[1 1 1])
    subplot(4,2,3)
        
        imagesc(pupildynamicsEMG.bins,pupildynamicsEMG.bins,pupildynamicsEMG.pWhisk');
        colormap(gca,emgcolor)
        axis xy
        xlim([-0.5 0.5]);ylim([-0.5 0.5])
        caxis([-0.01 1])
        colorbar
        %ColorbarWithAxis([-0.1 1],'EMG Envelope')
        xlabel('Pupil Diameter');ylabel('d/dt')
        title('Probability of Whisking')
        hold on
        %colormap(x,emgcolor)
        %colormap(y,'jet')
        plot(log10(whints_pupil(:,1)),whints_pupildt(:,1),'g.','markersize',0.1)
        %plot(whints_pupil(:,2),whints_pupildt(:,2),'r.')
        plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color',0.5.*[1 1 1])
        
    subplot(4,2,5)
        
        imagesc(pupildynamicsEMG.bins,pupildynamicsEMG.bins,(pupildynamicsEMG.Whstartrate)');
        colormap(gca,emgcolor)
        axis xy
        xlim([-0.5 0.5]);ylim([-0.5 0.5])
        %caxis([-0.01 1])
        %colormap(gca)
        colorbar
        %ColorbarWithAxis([-0.1 1],'EMG Envelope')
        xlabel('Pupil Diameter');ylabel('d/dt')
        title('p(WhiskStart)')
        hold on
        caxis([-0.1 1])
       % LogScale('c',10)
        %colormap(x,emgcolor)
        %colormap(y,'jet')
        plot(log10(whints_pupil(:,1)),whints_pupildt(:,1),'g.','markersize',0.1)
        %plot(whints_pupil(:,2),whints_pupildt(:,2),'r.')
        plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color',0.5.*[1 1 1])
        
    subplot(4,2,7)
        
        imagesc(pupildynamicsEMG.bins,pupildynamicsEMG.bins,log10(pupildynamicsEMG.meanWhdur)');
        %colormap(gca,emgcolor)
        axis xy
        xlim([-0.5 0.5]);ylim([-0.5 0.5])
        %caxis([-0.01 1])
        %colormap(gca)
        colorbar
        %ColorbarWithAxis([-0.1 1],'EMG Envelope')
        xlabel('Pupil Diameter');ylabel('d/dt')
        title('Whisk Duration')
        hold on
        caxis([-1 0.5])
        LogScale('c',10)
        %colormap(x,emgcolor)
        %colormap(y,'jet')
        plot(log10(whints_pupil(:,1)),whints_pupildt(:,1),'g.','markersize',0.1)
        %plot(whints_pupil(:,2),whints_pupildt(:,2),'r.')
        plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color',0.5.*[1 1 1])
        
        
        
    subplot(4,2,2)
        
        imagesc(pupilphaseEMG.bins,pupilphaseEMG.bins,(pupilphaseEMG.mean));
        hold on
        imagesc(pupilphaseEMG.bins+2*pi,pupilphaseEMG.bins,(pupilphaseEMG.mean));

        colormap(gca,emgcolor)
        axis xy
        ylim([-2 1]);xlim([-pi 3*pi])
        %caxis([-0.01 1])
        %colormap(gca)
        colorbar
        %ColorbarWithAxis([-0.1 1],'EMG Envelope')
        xlabel('Pupil Phase');ylabel('Power')
        title('EMG')
        hold on
        caxis([-0.75 0.5])
        LogScale('c',10)
        %colormap(x,emgcolor)
        %colormap(y,'jet')
        %plot(whints_pupil(:,2),whints_pupildt(:,2),'r.')
        plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color',0.5.*[1 1 1])
        
        
    subplot(4,2,4)
        
        imagesc(pupilphaseEMG.bins,pupilphaseEMG.bins,(pupilphaseEMG.pWhisk));
        hold on
        imagesc(pupilphaseEMG.bins+2*pi,pupilphaseEMG.bins,(pupilphaseEMG.pWhisk));

        colormap(gca,emgcolor)
        axis xy
        ylim([-2 1]);xlim([-pi 3*pi])
        %caxis([-0.01 1])
        %colormap(gca)
        colorbar
        %ColorbarWithAxis([-0.1 1],'EMG Envelope')
        xlabel('Pupil Phase');ylabel('Power')
        title('Probability of Whisking')
        hold on
        caxis([-0.01 1])

        %LogScale('c',10)
        %colormap(x,emgcolor)
        %colormap(y,'jet')
        %plot(whints_pupil(:,2),whints_pupildt(:,2),'r.')
        plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color',0.5.*[1 1 1])
        
    subplot(4,2,6)
        
        imagesc(pupilphaseEMG.bins,pupilphaseEMG.bins,(pupilphaseEMG.Whstartrate));
        hold on
        imagesc(pupilphaseEMG.bins+2*pi,pupilphaseEMG.bins,(pupilphaseEMG.Whstartrate));
        colormap(gca,emgcolor)
        axis xy
        ylim([-2 1]);xlim([-pi 3*pi])
        %caxis([-0.01 1])
        %colormap(gca)
        colorbar
        %ColorbarWithAxis([-0.1 1],'EMG Envelope')
        xlabel('Pupil Diameter');ylabel('d/dt')
        title('p(WhiskStart)')
        hold on
        caxis([-0.1 1])
       % LogScale('c',10)
        %colormap(x,emgcolor)
        %colormap(y,'jet')
        %plot(whints_pupil(:,2),whints_pupildt(:,2),'r.')
        plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color',0.5.*[1 1 1])
        
    subplot(4,2,8)
        
        imagesc(pupilphaseEMG.bins,pupilphaseEMG.bins,log10(pupilphaseEMG.meanWhdur));
        hold on
        imagesc(pupilphaseEMG.bins+2*pi,pupilphaseEMG.bins,log10(pupilphaseEMG.meanWhdur));
        %colormap(gca,emgcolor)
        axis xy
        ylim([-2 1]);xlim([-pi 3*pi])
        %caxis([-0.01 1])
        %colormap(gca)
        colorbar
        %ColorbarWithAxis([-0.1 1],'EMG Envelope')
        xlabel('Pupil Diameter');ylabel('d/dt')
        title('Whisk Duration')
        hold on
        caxis([-1 0.5])
        LogScale('c',10)
        %colormap(x,emgcolor)
        %colormap(y,'jet')
        %plot(whints_pupil(:,2),whints_pupildt(:,2),'r.')
        plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color',0.5.*[1 1 1])
        
NiceSave('EMGinPupilSpace',figfolder,baseName)
%% Figure: Whisk/EMG and Pupil
lagwin = [-5 5];
distcolor = [1 1 1; makeColorMap([0.7 0.7 0.7],[0 0.5 0],[0.7 0.6,0])];
windur = 200;
viewrange = bz_RandomWindowInIntervals(pupildilation.timestamps([1 end]),windur);
%viewrange = [400 600];
figure
    subplot(6,3,9)
        plot(pupilEMGcorr.corrlags,pupilEMGcorr.xcorr','k')
        hold on
        box off
        axis tight
        xlim([-40 40])
        plot([0 0],get(gca,'ylim'),'b-')
        xlabel('t lag (s)');ylabel('EMG')
        title('Pupil-EMG CCG')
        
    subplot(6,3,1:2)
        plot(pupildilation.timestamps,pupildilation.interpdata,'k','linewidth',2)
        xlim(viewrange);ylim([0 4])
        %xlabel('t (s)')
        box off
         ylabel('Pupil')
        set(gca,'xticklabel',[])

    subplot(6,3,4:5)
        plot(EMGwhisk.timestamps,EMGwhisk.EMG,'color',[0.5 0.5 0.5],'linewidth',0.5)
        hold on
        plot(pupildilation.timestamps,EMGwhisk.pupiltime,'b','linewidth',2)
        plot(EMGwhisk.ints.Wh',...
            ones(size(EMGwhisk.ints.Wh))',...
            'g-','linewidth',1)
        axis tight
        xlim(viewrange);
        xlabel('t (s)')
        ylabel('EMG')
        box off

    subplot(3,3,7)
        %colormap(gca,distcolor)
        colormap(gca,[1 1 1;colormap])
        imagesc(pupilEMGdist.bins{1},pupilEMGdist.bins{2},...
            (pupilEMGdist.counts)'./max(pupilEMGdist.counts(:)))
        hold on
        plot(get(gca,'xlim'),log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],'g--')
        
        axis xy
        colorbar
        LogScale('y',10)
        caxis([0 0.6])
        %ylim([-1 1.5])
        %xlim([0 1])
        xlabel('Pupil Diameter (0-1)');ylabel('EMG Envelope')

    subplot(6,3,3)
        plot(t_lag,pwCCG.pupil.WhOn,'k','linewidth',2)
        axis tight
        box off
        hold on
        %plot(t_lag,Wh_on_mean+Wh_on_std,'k--')
        plot([0 0],get(gca,'ylim'),'g')
        xlim(lagwin)
        ylabel('Pupil')
        set(gca,'xticklabel',[])

    subplot(6,3,6)
        plot(t_lag,pwCCG.EMG.WhOn,'b','linewidth',2)
        axis tight
        box off
        hold on
        %plot(t_lag,Wh_on_mean+Wh_on_std,'k--')
        plot([0 0],get(gca,'ylim'),'g')
        xlim(lagwin)
        xlabel('t (s - aligned to WhOn)')
        ylabel('EMG')
        
        
    subplot(3,3,8)
        
        imagesc(pupildynamicsEMG.bins,pupildynamicsEMG.bins,pupildynamicsEMG.pWhisk');
        colormap(gca,emgcolor)
        axis xy
        xlim([-0.5 0.5]);ylim([-0.5 0.5])
        caxis([-0.01 1])
        colorbar
        %ColorbarWithAxis([-0.1 1],'EMG Envelope')
        xlabel('Pupil Diameter');ylabel('d/dt')
        title('Probability of Whisking')
        hold on
        %colormap(x,emgcolor)
        %colormap(y,'jet')
        plot(log10(whints_pupil(:,1)),whints_pupildt(:,1),'g.','markersize',0.1)
        %plot(whints_pupil(:,2),whints_pupildt(:,2),'r.')
        plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color',0.5.*[1 1 1])

    subplot(6,4,10)
        %colormap(gca,emgcolor)
        scatter(log10(pupildilation.data(1:end-1)),pupildilation.dpdt,0.2,log10(EMGwhisk.pupiltime(1:end-1)))
        %plot(pupildilation.data(1:end-1),pupildilation.dpdt,'k.','markersize',2)
        ylim([-0.5 0.5]);xlim([-0.5 0.5])
        LogScale('x',10)
        %
        %colorbar
        ColorbarWithAxis([-1.5 1.8],'EMG Envelope')
        xlabel('Pupil Diameter');ylabel('d/dt')
        
%     subplot(6,4,10)
%         colormap(gca,emgcolor)
%         imagesc(pupildynamicsEMG.bins,pupildynamicsEMG.bins,pupildynamicsEMG.mean')
%         axis xy
%         xlim([0 1]);ylim([-0.2 0.2])
%         %colorbar
%         ColorbarWithAxis([-1.5 1.5],'EMG Envelope')
%         xlabel('Pupil Diameter');ylabel('d/dt')


subplot(4,3,9)
%scatter(pupildilation.data(~isnan(pupildilation.interpdata)),pupildilation.dpdt(~isnan(pupildilation.interpdata)),0.2,lowpupildata.phase)
scatter(lowpupildata.phase,log10(lowpupildata.amp),0.2,log10(EMGwhisk.pupiltime(~isnan(pupildilation.interpdata))))
hold on
scatter(lowpupildata.phase+2*pi,log10(lowpupildata.amp),0.2,log10(EMGwhisk.pupiltime(~isnan(pupildilation.interpdata))))
axis tight 
caxis([-1.5 1.8])
xlabel('Phase');ylabel('Amplitude')

subplot(4,3,12)
colormap(gca,emgcolor(2:end,:))
%scatter(pupildilation.data(~isnan(pupildilation.interpdata)),pupildilation.dpdt(~isnan(pupildilation.interpdata)),0.2,lowpupildata.phase)
scatter(lowpupildata.phase,log10(lowpupildata.amp),0.2,(pupildilation.iswhisk(~isnan(pupildilation.interpdata))))
hold on
scatter(lowpupildata.phase+2*pi,log10(lowpupildata.amp),0.2,(pupildilation.iswhisk(~isnan(pupildilation.interpdata))))
axis tight 
xlabel('Phase');ylabel('Amplitude')

        
NiceSave('WhiskAndPupil',figfolder,baseName)



%% Pupil-Whisk Phase coupling

frange = [0.001 1];
[wavespec] = bz_WaveSpec(pupildilation.interpdata,...
    'frange',frange,'nfreqs',100,'ncyc',4,'samplingRate',pupildilation.samplingRate);
wavespec.timestamps=pupildilation.timestamps;

%%
whpow = EMGwhisk.pupiltime./mean(EMGwhisk.pupiltime);
powbins = 10;
%sig2powerskew = zeros(powbins,wavespec.nfreqs);
coupling = zeros(wavespec.nfreqs,1);
for ff = 1:wavespec.nfreqs
%ff
% [ ampbins,phasebins,sig2powerskew(:,ff),~,~ ] = ...
%     PhaseAmpCouplingByAmp( angle(wavespec.data(:,ff)),log10(abs(wavespec.data(:,ff))),...
%     log10(EMGwhisk.pupiltime),powbins );
%pow = abs(wavespec.data(:,ff))./mean(abs(wavespec.data(:,ff)));
coupling(ff) = abs(mean(whpow.*exp(1i.*angle(wavespec.data(:,ff)))));

%close all
end
WPcoupling.coupling = coupling;
WPcoupling.freqs = wavespec.freqs;
%%
figure
plot(log10(wavespec.freqs),coupling)
LogScale('x',10)
xlabel('f (Hz)')
ylabel('Pupil Phase - Whisk Coupling')
NiceSave('PupilWhiskCoupling',figfolder,baseName)
%%
% figure
% imagesc(log10(wavespec.freqs),ampbins,sig2powerskew)
% LogScale('x',10)
% axis xy
% colorbar
%%
% 
% %Find frames closest to WHon, NWHon transitions
% Wh_ints = interp1(pupildilation.timestamps,pupildilation.timestamps,EMGwhisk.ints.Wh,'nearest');
% NWh_ints = interp1(pupildilation.timestamps,pupildilation.timestamps,EMGwhisk.ints.NWh,'nearest');
% %Get average pupil, EMG aligned to transitions (20s lag)
% 
% t_lag = 20; %s
% t_lag_frames = round(t_lag./(1./pupildilation.samplingRate));
% t_lag = [-t_lag_frames:t_lag_frames].*(1./pupildilation.samplingRate);
% 
% Wh_onints = find(ismember(pupildilation.timestamps,Wh_ints(:,1)));
% Wh_onints = bsxfun(@(A,B) A+B,Wh_onints,t_lag_frames.*[-1 1]);
% 
% [Wh_onepochs,~] = IsolateEpochs2(pupildilation.data,Wh_onints,0,1)
% Wh_onepochs_pupil = cat(2,Wh_onepochs{:});
% Wh_on_mean_pupil = mean(Wh_onepochs_pupil,2);
% %Wh_on_std = std(Wh_onepochs,[],2);
% 
% [Wh_onepochs,~] = IsolateEpochs2(EMGwhisk.pupiltime,Wh_onints,0,1)
% Wh_onepochs_EMG = cat(2,Wh_onepochs{:});
% Wh_on_mean_EMG = mean(Wh_onepochs_EMG,2);
% %Wh_on_std = std(Wh_onepochs,[],2);
% 
% 
% %%
% 
% [pwCCG.pupil.WhOff,t_lag,~,alltrans.pupil.WhOff] = EventVsContinousCCG(pupildilation.data,pupildilation.timestamps,EMGwhisk.ints.Wh(:,2),20);
% [pwCCG.EMG.WhOff,t_lag,~,alltrans.EMG.WhOff] = EventVsContinousCCG(EMGwhisk.pupiltime,pupildilation.timestamps,EMGwhisk.ints.Wh(:,2),20);
% 
% [pwCCG.pupil.NWhOn,t_lag,~,alltrans.pupil.NWhOn] = EventVsContinousCCG(pupildilation.data,pupildilation.timestamps,EMGwhisk.ints.NWh(:,1),20);
% [pwCCG.EMG.NWhOn,t_lag,~,alltrans.EMG.NWhOn] = EventVsContinousCCG(EMGwhisk.pupiltime,pupildilation.timestamps,EMGwhisk.ints.NWh(:,1),20);
% 
% [pwCCG.pupil.NWhOff,t_lag,~,alltrans.pupil.NWhOff] = EventVsContinousCCG(pupildilation.data,pupildilation.timestamps,EMGwhisk.ints.NWh(:,2),20);
% [pwCCG.EMG.NWhOff,t_lag,~,alltrans.EMG.NWhOff] = EventVsContinousCCG(EMGwhisk.pupiltime,pupildilation.timestamps,EMGwhisk.ints.NWh(:,2),20);
% 
% %%
% 
% figure
% subplot(5,2,1)
% plot(t_lag,pwCCG.pupil.WhOn,'k','linewidth',2)
% axis tight
% box off
% hold on
% %plot(t_lag,Wh_on_mean+Wh_on_std,'k--')
% plot([0 0],get(gca,'ylim'),'g')
% xlim(lagwin)
% 
% subplot(5,2,3)
% plot(t_lag,pwCCG.EMG.WhOn,'b','linewidth',2)
% axis tight
% box off
% hold on
% %plot(t_lag,Wh_on_mean+Wh_on_std,'k--')
% plot([0 0],get(gca,'ylim'),'g')
% xlim(lagwin)
% 
% subplot(5,2,2)
% plot(t_lag,pwCCG.pupil.WhOff,'k','linewidth',2)
% axis tight
% box off
% hold on
% %plot(t_lag,Wh_on_mean+Wh_on_std,'k--')
% plot([0 0],get(gca,'ylim'),'g')
% xlim(lagwin)
% 
% subplot(5,2,4)
% plot(t_lag,pwCCG.EMG.WhOff,'b','linewidth',2)
% axis tight
% box off
% hold on
% %plot(t_lag,Wh_on_mean+Wh_on_std,'k--')
% plot([0 0],get(gca,'ylim'),'g')
% xlim(lagwin)
% 
% 
% subplot(5,2,5)
% plot(t_lag,pwCCG.pupil.NWhOn,'k','linewidth',2)
% axis tight
% box off
% hold on
% %plot(t_lag,Wh_on_mean+Wh_on_std,'k--')
% plot([0 0],get(gca,'ylim'),'r')
% xlim(lagwin)
% 
% subplot(5,2,7)
% plot(t_lag,pwCCG.EMG.NWhOn,'b','linewidth',2)
% axis tight
% box off
% hold on
% %plot(t_lag,Wh_on_mean+Wh_on_std,'k--')
% plot([0 0],get(gca,'ylim'),'r')
% xlim(lagwin)
% 
% subplot(5,2,6)
% plot(t_lag,pwCCG.pupil.NWhOff,'k','linewidth',2)
% axis tight
% box off
% hold on
% %plot(t_lag,Wh_on_mean+Wh_on_std,'k--')
% plot([0 0],get(gca,'ylim'),'r')
% xlim(lagwin)
% 
% subplot(5,2,8)
% plot(t_lag,pwCCG.EMG.NWhOff,'b','linewidth',2)
% axis tight
% box off
% hold on
% %plot(t_lag,Wh_on_mean+Wh_on_std,'k--')
% plot([0 0],get(gca,'ylim'),'r')
% xlim(lagwin)
% 
% subplot(5,2,9)
% plotyy(t_lag,pwCCG.pupilxy.WhOn(:,1),t_lag,pwCCG.pupilxy.WhOn(:,2))
% %hold on
% %plot(t_lag,pwCCG.pupilxy.WhOn(:,2),'k','linewidth',2)
% %axis tight
% box off
% hold on
% %plot(t_lag,Wh_on_mean+Wh_on_std,'k--')
% plot([0 0],get(gca,'ylim'),'r')
% xlim(lagwin)
% 
% NiceSave('PupilEMG_Whisk',figfolder,baseName)
% 
% %%
% Whdur = diff(EMGwhisk.ints.Wh,1,2);
% [dursorted,sortwhdur] = sort(Whdur(~skippedWh));
% 
% %%
% pupilWhnorm = NormToInt(alltrans.pupil.WhOn,'modZ',[-5 0]-t_lag(1),1./(1./pupildilation.samplingRate));
% %%
% figure
% subplot(3,1,1)
%     imagesc(t_lag,[1 length(sortwhdur)],alltrans.EMG.WhOn(:,sortwhdur)')
%     hold on
%     plot(dursorted,1:length(sortwhdur),'r')
%     plot([0 0],[1 length(sortwhdur)],'r')
%     ColorbarWithAxis([0 20],'EMG Envelope')
%     caxis([0 20])
%     xlim(lagwin)
% subplot(3,1,2)
%     imagesc(t_lag,[1 length(sortwhdur)],(alltrans.pupil.WhOn(:,sortwhdur))')
%     hold on
%     plot(dursorted,1:length(sortwhdur),'r')
%     plot([0 0],[1 length(sortwhdur)],'r')
%     ColorbarWithAxis([0 1],{'Pupil Diameter', '(0-1 normalized to whole recording)'})
%     xlim(lagwin)
%     
% subplot(3,1,3)
%     imagesc(t_lag,[1 length(sortwhdur)],pupilWhnorm(:,sortwhdur)')
%     hold on
%     plot(dursorted,1:length(sortwhdur),'r')
%     plot([0 0],[1 length(sortwhdur)],'r')
%     colorbar
%     xlim(lagwin)
%     ColorbarWithAxis([0 5],{'Pupil Diameter', '(z-normalized to 5s before Whisking)'})
%     xlabel('t (s, aligned to whisking onset)')
%     NiceSave('WhiskEpochs',figfolder,baseName)
% %%
% 
% %%
% 
% 
% 
% 
% 
% 
% 
% 
% figure
% 
% 
% subplot(2,2,2)
% colormap(gca,distcolor)
% imagesc(pupildynamicsEMG.bins,pupildynamicsEMG.bins,(pupildynamicsEMG.num)'./max(pupildynamicsEMG.num(:)))
% axis xy
% xlim([0 1]);ylim([-0.2 0.2])
% colorbar
% %ColorbarWithAxis([-1.5 1.8],'EMG Envelope')
% xlabel('Pupil Diameter');ylabel('d/dt')
% 
% subplot(2,2,3)
% imagesc(pupildynamicsEMG.bins,pupildynamicsEMG.bins,(pupildynamicsEMG.std)')
% axis xy
% colorbar
% caxis([0 1])
% xlim([0 1]);ylim([-0.2 0.2])
% %ColorbarWithAxis([-1.5 1.8],'EMG Envelope')
% xlabel('Pupil Diameter');ylabel('d/dt')
% 
% 
% 
% NiceSave('PupilEMG_dddt',figfolder,baseName)
% %%
% emgcolor = makeColorMap([0.5 0.5 0.5],[0 0 0],[0 0 1]);
% figure
% subplot(4,1,1)
%     plot(pupildilation.timestamps,pupildilation.data,'k')
% subplot(4,1,2)
%     plot(pupildilation.timestamps(1:end-1),pupildilation.dpdt,'k')
%     ylim([-1 1])
% 
% subplot(2,2,3)
% %colormap(gca,emgcolor)
% scatter(pupildilation.data(1:end-1),pupildilation.dpdt,0.2,log10(EMGwhisk.pupiltime(1:end-1)))
% %plot(pupildilation.data(1:end-1),pupildilation.dpdt,'k.','markersize',2)
% ylim([-0.5 0.5])
% %
% %colorbar
% ColorbarWithAxis([-1.5 1.8],'EMG Envelope')
% xlabel('Pupil Diameter');ylabel('d/dt')


end

