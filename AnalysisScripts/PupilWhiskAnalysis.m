function [puphist,pupACG,pupPSD,...
    EMGhist,Whdurhist,...
    pupilEMGdist,pupildynamicsEMG,...
    pupilEMGcorr,pwCCG] = PupilWhiskAnalysis(basePath,figfolder)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%%
%basePath = '/mnt/proraidDL/Database/WMProbeData/170421_Layers_LFP_Pupil_EMG_Emx1M1M3/170421_Layers_LFP_Pupil_EMG_170421_180427/';
%basePath = '/mnt/proraidDL/Database/WMProbeData/170606_Layers_LFP_Pupil_EMG_Emx1M2M4/170606_Layers_LFP_Pupil_EMG_170606_204730/';
%basePath= '/mnt/proraidDL/Database/WMProbeData/170421_Layers_LFP_Pupil_EMG_Emx1M1M3/170421_Layers_LFP_Pupil_EMG_170421_180427';

baseName = bz_BasenameFromBasepath(basePath);
%figfolder = '/mnt/data1/Dropbox/research/Current Projects/S1State/AnalysisScripts/figures/PupilWhiskAnalysis';

%% PUPIL %%
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
pupdthist.counts = hist3([pupildilation.interpdata(1:end-1),pupildt],pupdthist.bins);
pupdthist.counts = pupdthist.counts./sum(pupdthist.counts(:));


%% Calculate FFT of the pupil
frange = [0.005 1];
[freqs,t,spec] = WaveSpec(pupildilation.interpdata(~isnan(pupildilation.interpdata)),frange,200,3,1/pupildilation.samplingRate,'log');
spec = (abs(spec));
%%
lowfilter = [0.015 0.08];
highfilter = [0.3 0.8];

pupil4filter = pupildilation;
pupil4filter.data = pupildilation.interpdata(~isnan(pupildilation.interpdata));
%pupil4filter.t
pupil4filter.timestamps = pupil4filter.timestamps(~isnan(pupildilation.interpdata));
lowpupildata = bz_Filter(pupil4filter,'passband',lowfilter,'filter' ,'fir1','order',3);
%highpupildata = bz_Filter(pupil4filter,'passband',highfilter,'filter' ,'fir1');

pupPSD.freqs = freqs;
pupPSD.psd = log10(mean(spec,2));
%% Figure: PUPIL
winsize = 300; %s
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
        plot(lowpupildata.timestamps,lowpupildata.data+nanmean(pupildilation.data),'g')
      %  plot(highpupildata.timestamps,highpupildata.data+nanmean(pupildilation.data),'r')
        plot(get(gca,'xlim'),[0 0],'r-')
        plot(pupildilation.timestamps(1:end-1),pupildt,'k','linewidth',1)
        xlim(100+[0 winsize])
        ylim([-0.15 1])
        ylabel('Pupil: Zoom')
        xlabel('t (s)')
        
    subplot(6,4,4)
        bar(puphist.bins,puphist.counts,'facecolor','k')
        xlim([0 1])
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
    colormap(gca,[1 1 1;colormap])
        imagesc(pupdthist.bins{1},pupdthist.bins{2},pupdthist.counts')
        hold on
        plot(get(gca,'xlim'),[0 0],'r-')
        colorbar
        xlabel('Pupil Diameter');ylabel('d/dt')
        title('Pupil Dynamics Histogram')
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
        plot(log10(freqs),log10(mean(spec,2)),'k','linewidth',2)
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

%% EMG %%
EMGwhisk = bz_LoadStates( basePath,'EMGwhisk');
if isempty(EMGwhisk)
    EMGwhisk = GetWhiskFromEMG(basePath);
end

%%
EMGwhisk.pupiltime = interp1(EMGwhisk.t,EMGwhisk.EMGsm,pupildilation.timestamps);

%% Histogram
EMGhist.bins = linspace(0,20,20);
EMGhist.counts = hist(EMGwhisk.EMGsm,EMGhist.bins);
EMGhist.counts =EMGhist.counts./sum(EMGhist.counts);

EMGhist.logbins = linspace(-1,1.5,20);
EMGhist.logcounts = hist(log10(EMGwhisk.EMGsm),EMGhist.logbins);
EMGhist.logcounts =EMGhist.logcounts./sum(EMGhist.logcounts);

%% Whisk Durations
EMGwhisk.Whdurs = diff(EMGwhisk.ints.Wh,1,2);
EMGwhisk.InterWhdurs = EMGwhisk.ints.Wh(2:end,1)-EMGwhisk.ints.Wh(1:end-1,2);

Whdurhist.bins = linspace(-1,2,20);
Whdurhist.Whdurs = hist(log10(EMGwhisk.Whdurs),Whdurhist.bins);
Whdurhist.Whdurs = Whdurhist.Whdurs./sum(Whdurhist.Whdurs);
Whdurhist.InterWhdurs = hist(log10(EMGwhisk.InterWhdurs),Whdurhist.bins);
Whdurhist.InterWhdurs = Whdurhist.InterWhdurs./sum(Whdurhist.InterWhdurs);

%% Figure: EMG/WHISK
figure
    subplot(4,1,1)
        plot(EMGwhisk.t,EMGwhisk.EMG,'color',[0.5 0.5 0.5],'linewidth',0.5)
        hold on
        plot(pupildilation.timestamps,EMGwhisk.pupiltime,'b','linewidth',2)
        plot(EMGwhisk.ints.Wh',...
            EMGwhisk.detectorparms.Whthreshold.*ones(size(EMGwhisk.ints.Wh))',...
            'g-','linewidth',1)
        %plot(EMGwhisk.t,EMGwhisk.EMGsm,'b')
        ylabel('EMG: All')
        axis tight
        xlabel('t (s)')
        box off
        xlim(pupildilation.timestamps([1 end]))
        ylim([-20 40])
    subplot(4,1,2)
        plot(EMGwhisk.t,EMGwhisk.EMG,'color',[0.5 0.5 0.5],'linewidth',0.5)
        hold on
        plot(pupildilation.timestamps,EMGwhisk.pupiltime,'b','linewidth',2)
        plot(EMGwhisk.ints.Wh',...
            EMGwhisk.detectorparms.Whthreshold.*ones(size(EMGwhisk.ints.Wh))',...
            'g-','linewidth',1)
        %plot(EMGwhisk.t,EMGwhisk.EMGsm,'b')
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
pupilEMGdist.bins = {linspace(0,1,50),linspace(-1,1.5,150)};
[pupilEMGdist.counts,pupilEMGdist.bins] = hist3([pupildilation.data,log10(EMGwhisk.pupiltime)],pupilEMGdist.bins);

%XCovariance
[pupilEMGcorr.xcorr,pupilEMGcorr.corrlags] = xcov(pupildilation.interpdata(~isnan(pupildilation.interpdata)),EMGwhisk.pupiltime(~isnan(pupildilation.interpdata)),'unbiased');
pupilEMGcorr.corrlags = pupilEMGcorr.corrlags.*(1./pupildilation.samplingRate);

%Pupil/EMG at whisking onset
[pwCCG.pupil.WhOn,t_lag,~,alltrans.pupil.WhOn,skippedWh] = EventVsContinousCCG(pupildilation.data,pupildilation.timestamps,EMGwhisk.ints.Wh(:,1),20);
[pwCCG.pupilxy.WhOn(:,1),t_lag,~,~] = EventVsContinousCCG(pupildilation.pupilxy(:,1),pupildilation.timestamps,EMGwhisk.ints.Wh(:,1),20);
[pwCCG.pupilxy.WhOn(:,2),t_lag,~,~] = EventVsContinousCCG(pupildilation.pupilxy(:,2),pupildilation.timestamps,EMGwhisk.ints.Wh(:,1),20);
[pwCCG.EMG.WhOn,t_lag,~,alltrans.EMG.WhOn] = EventVsContinousCCG(EMGwhisk.pupiltime,pupildilation.timestamps,EMGwhisk.ints.Wh(:,1),20);

%% Histogram: Whisking by Pupil Dynamics
[pupildynamicsEMG,pupildynamicsEMG.bins]=PairMatHist(log10(EMGwhisk.pupiltime(1:end-1)),[pupildilation.data(1:end-1),pupildt],100,[-1 1]);

[~,EMGwhisk.IDXvec] = RestrictInts(pupildilation.timestamps,EMGwhisk.ints.Wh);
[pWhisk]=PairMatHist(single(EMGwhisk.IDXvec(1:end-1)),[pupildilation.data(1:end-1),pupildt],100,[-1 1]);
pupildynamicsEMG.pWhisk = pWhisk.mean;

%Wh onsets/offsets in pupil space
whints_pupil = interp1(pupildilation.timestamps,pupildilation.data,EMGwhisk.ints.Wh);
whints_pupildt = interp1(pupildilation.timestamps(1:end-1),pupildt,EMGwhisk.ints.Wh);

%% Figure: Whisk/EMG and Pupil
lagwin = [-5 5];
distcolor = [1 1 1; makeColorMap([0.7 0.7 0.7],[0 0.5 0],[0.7 0.6,0])];
emgcolor = [1 1 1;makeColorMap([0.5 0.5 0.5],[0 0 0.8])];
viewrange = [400 600];
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
        xlim(viewrange);ylim([0 1])
        %xlabel('t (s)')
        box off
         ylabel('Pupil')
        set(gca,'xticklabel',[])

    subplot(6,3,4:5)
        plot(EMGwhisk.t,EMGwhisk.EMG,'color',[0.5 0.5 0.5],'linewidth',0.5)
        hold on
        plot(pupildilation.timestamps,EMGwhisk.pupiltime,'b','linewidth',2)
        plot(EMGwhisk.ints.Wh',...
            EMGwhisk.detectorparms.Whthreshold.*ones(size(EMGwhisk.ints.Wh))',...
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
        colormap(gca,emgcolor)
        imagesc(pupildynamicsEMG.bins,pupildynamicsEMG.bins,pupildynamicsEMG.pWhisk')
        axis xy
        xlim([0 1]);ylim([-0.2 0.2])
        caxis([-0.01 1])
        colorbar
        %ColorbarWithAxis([-0.1 1],'EMG Envelope')
        xlabel('Pupil Diameter');ylabel('d/dt')
        title('Probability of Whisking')
        hold on
        %plot(whints_pupil(:,1),whints_pupildt(:,1),'g.','markersize',0.1)
        %plot(whints_pupil(:,2),whints_pupildt(:,2),'r.')
        plot(get(gca,'xlim'),[0 0],'--','linewidth',0.5,'color',0.5.*[1 1 1])

    subplot(6,4,10)
        %colormap(gca,emgcolor)
        scatter(pupildilation.data(1:end-1),pupildt,0.2,log10(EMGwhisk.pupiltime(1:end-1)))
        %plot(pupildilation.data(1:end-1),pupildt,'k.','markersize',2)
        ylim([-0.15 0.15])
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
%scatter(pupildilation.data(~isnan(pupildilation.interpdata)),pupildt(~isnan(pupildilation.interpdata)),0.2,lowpupildata.phase)
scatter(lowpupildata.phase,lowpupildata.amp,0.2,log10(EMGwhisk.pupiltime(~isnan(pupildilation.interpdata))))
hold on
scatter(lowpupildata.phase+2*pi,lowpupildata.amp,0.2,log10(EMGwhisk.pupiltime(~isnan(pupildilation.interpdata))))
axis tight 
caxis([-1.5 1.8])
xlabel('Phase');ylabel('Amplitude')

subplot(4,3,12)
colormap(gca,emgcolor(2:end,:))
%scatter(pupildilation.data(~isnan(pupildilation.interpdata)),pupildt(~isnan(pupildilation.interpdata)),0.2,lowpupildata.phase)
scatter(lowpupildata.phase,lowpupildata.amp,0.2,(EMGwhisk.IDXvec(~isnan(pupildilation.interpdata))))
hold on
scatter(lowpupildata.phase+2*pi,lowpupildata.amp,0.2,(EMGwhisk.IDXvec(~isnan(pupildilation.interpdata))))
axis tight 
xlabel('Phase');ylabel('Amplitude')

        
NiceSave('WhiskAndPupil',figfolder,baseName)
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
%     plot(pupildilation.timestamps(1:end-1),pupildt,'k')
%     ylim([-1 1])
% 
% subplot(2,2,3)
% %colormap(gca,emgcolor)
% scatter(pupildilation.data(1:end-1),pupildt,0.2,log10(EMGwhisk.pupiltime(1:end-1)))
% %plot(pupildilation.data(1:end-1),pupildt,'k.','markersize',2)
% ylim([-0.5 0.5])
% %
% %colorbar
% ColorbarWithAxis([-1.5 1.8],'EMG Envelope')
% xlabel('Pupil Diameter');ylabel('d/dt')


end

