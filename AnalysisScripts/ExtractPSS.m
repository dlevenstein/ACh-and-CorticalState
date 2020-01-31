function [  ] = ExtractPSS(basePath,figfolder)
%basePath = '/mnt/proraidDL/Database/WMData/AChPupil/EM1M3/180209_WT_EM1M3';

[lfp] = bz_GetLFP('all','basepath',basePath);
%%
winsize = 1;
dt = 0.1;
frange = [2.5 100];
[~] = bz_PowerSpectrumSlope(lfp,winsize,dt,'nfreqs',100,...
    'IRASA',true,'Redetect',true,'frange',frange,'saveMat',basePath);

winsize = 8;
dt = 0.02;
frange = [2 128];
[~] = bz_PowerSpectrumSlope(lfp,winsize,dt,...
    'IRASA',true,'Redetect',true,'frange',frange,'saveMat',basePath,'saveName','_wav',...
    'spectype','wavelet','nfreqs',100);

%% PSS example figure
% if strcmp(baseName,'171209_WT_EM1M3')
%     PSS.spec = PSpecSlope.Shortwin.SPEC(:,:,PSS.CTXchans);
%     PSS.frac = PSpecSlope.Shortwin.FRAC(:,:,PSS.CTXchans);
% %exwins = [600 775];
% %exwins = [lfp.timestamps(1) lfp.timestamps(end)];
% %exwins = [600 666];
% %exwins = [100 250];
% exwins = [1175 1300;1242 1270;1260 1266;1244 1250];
% for ww = 1:length(exwins)
% %for ww = 2
%     timewin = exwins(ww,:);
% 
% lfpinwin = InIntervals(lfp.timestamps,timewin);
% %pssinwin = InIntervals(PSS.timestamps,timewin);
% exchan_spec = strcmp(spec.chanlayers,'L5b6');
% exchan_PSS = find(spec.channels(exchan_spec)==PSS.chan);
% 
% speclim = [0.65 1.35]; %Med norm
% PSSrange = [-2.4 -1.2];
% %PSSrange = [-2.6 -0.8];
% 
% figure
% subplot(8,1,1:2)
% %plot(pupildilation.timestamps,pupildilation.data_raw,'k--','Linewidth',2)
% hold on
% plot(pupildilation.timestamps,pupildilation.data,'r','Linewidth',2)
% hold on
% plot(EMGwhisk.timestamps,EMGwhisk.EMG./50+0.4,'color',0.5.*[1 1 1])
% plot(EMGwhisk.timestamps,EMGwhisk.EMGsm./9+0.4,'k')
% ylim([0.35 2.5])
% if ww>7
%    plot(mean(timewin).*[1 1], ylim(gca),'k')
% end
% 
% xlim(timewin)
% %colorbar
% box off
% bz_ScaleBar('s')
% 
% subplot(8,1,3:4)
% imagesc(spec.timestamps,log10(spec.freqs),spec.data(:,:,exchan_spec)')
% hold on
% ylim([-1 2.5])
% %plot(PSS.timestamps,bz_NormToRange(PSS.data(:,exchan_PSS),[0 2.5]),'k','LineWidth',2)
% plot(lfp.timestamps(lfpinwin),bz_NormToRange(single(lfp.data(lfpinwin,exchan_spec)),[-1 -0.1]),'k','LineWidth',0.1)
% axis xy
% if ww>7
%    plot(mean(timewin).*[1 1], ylim(gca),'k')
% end
% LogScale('y',10)
%         %ColorbarWithAxis(speclim,'Power (med^-^1)')
%         clim(speclim);
%        box off
%        crameri vik
% 
% yyaxis right
% plot(PSS.timestamps,PSS.data(:,exchan_PSS),'k','LineWidth',2)
% ylim([-3.5 -0.85])
% xlim(timewin)
%     bz_ScaleBar('s')
% 
% depthscalefact=0.5e5;
% if ww>2
%     lfpscalefact = 1.5;
% else
%     lfpscalefact = 1;
% end
% 
% 
% subplot(2,1,2)
% imagesc(PSS.timestamps,PSS.interpdepth*depthscalefact,PSS.depthinterp')
% axis xy
% hold on
% bz_MultiLFPPlot(lfp,'timewin',timewin,'LFPlabels',lfp.chanlayers,...
%     'LFPmidpoints',lfp.chandepths*depthscalefact,'lfpcolor','k',...
%     'lfpwidth',0.1,'scaleLFP',lfpscalefact);
% %ColorbarWithAxis(PSSrange,'PSS')
% clim(PSSrange)
% ylim([-1.05 -0.025]*depthscalefact)
% 
%     %crameri bamako
% if ww>7
%    plot(mean(timewin).*[1 1], ylim(gca),'k')
% end    
% xlim(timewin)
%     bz_ScaleBar('s')   
% 
% NiceSave(['IllustrateExample',num2str(ww)],figfolder,baseName)
% 
% if ww == 2
% figure
%     subplot(8,1,7:8)
%         imagesc(spec.timestamps,log10(spec.freqs),spec.data(:,:,exchan_spec)')
%         hold on
%         ylim([-1 2.5])
%         %plot(PSS.timestamps,bz_NormToRange(PSS.data(:,exchan_PSS),[0 2.5]),'k','LineWidth',2)
%         plot(lfp.timestamps(lfpinwin),bz_NormToRange(single(lfp.data(lfpinwin,exchan_spec)),[-1 -0.1]),...
%             'k','LineWidth',0.1)
%         axis xy
%         if ww>7
%            plot(mean(timewin).*[1 1], ylim(gca),'k')
%         end
%         LogScale('y',10)
%                 %ColorbarWithAxis(speclim,'Power (med^-^1)')
%                 clim(speclim);
%                box off
%                crameri vik
%         xlim(timewin)
%         bz_ScaleBar('s')
%         yyaxis right
%         plot(PSS.timestamps,PSS.data(:,exchan_PSS),'k','LineWidth',2)
%         ylim([-3.5 -0.85])
%     
%         
%         times = (0.5+timewin(1)):2:timewin(2);
%         usefreqs = PSS.freqs>=1 & PSS.freqs<=100;
%         fracfreqs = PSS.freqs>=2 & PSS.freqs<=100;
%        
%     subplot(4,1,1)
%     for tt = 1:length(times)
%         hold on
%         [~,timepoint] = min(abs(PSS.timestamps-(times(tt))));
%         
%         linefit.validfreq = PSS.freqs; linefit.frac = PSS.frac(:,timepoint,exchan_PSS)';
%         [linefit] = WaveIRASA_plawfit( linefit, [2.5 100] );
%         
%         plot(log10(PSS.freqs(usefreqs))+times(tt)-0.5,log10(PSS.spec(usefreqs,timepoint,exchan_PSS)),'k','linewidth',1)
%         plot(log10(linefit.fitFreq)+times(tt)-0.5,linefit.Beta*log10(linefit.fitFreq)+linefit.Cons,'r','linewidth',2)
%     
%     end
%     %clim(PSSrange)
%         axis tight
%         box off
%         xlabel('f (Hz)');ylabel('Power (dB)')
%         %LogScale('x',10)
%         xlim(timewin)
%         %xlim(timewin)
%             NiceSave(['IllustrateExamplePSS',num2str(ww)],figfolder,baseName)   
% 
% end  
% end
% end

end