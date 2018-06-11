function [specvarcorr,specbyvar,specvarcorr_meannorm,specbyvar_meannorm,...
    pupilLFPcorr,LFPbypupil_meannorm,LFPbypupil] = LFPSpecByPupilAnalysis(basePath,figfolder)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%%
%basePath = pwd;
baseName = bz_BasenameFromBasepath(basePath);
figfolder = '/mnt/data1/Dropbox/research/Current Projects/S1State/AnalysisScripts/figures/LFPSpecByPupilAnalysis';


%%
recparms = bz_getSessionInfo(basePath,'noPrompts',true);



%%
pupildilation = bz_LoadBehavior(basePath,'pupildiameter');
% pupildilation.timestamps = pupildilation.t_interp;
% pupildilation.maxnormdata = pupildilation.puparea_pxl./max(pupildilation.puparea_pxl);
% pupildilation.sf = 1./(diff(pupildilation.t_pulse([1 2])));

%%

nantimes = isnan(pupildilation.data);
pupildilation.interpdata = interp1(pupildilation.timestamps(~nantimes),...
    pupildilation.data(~nantimes),pupildilation.timestamps);
%% Pupil ACG
acgwin = 150; %s
[pupilACG,tlag] = xcorr(pupildilation.interpdata,...
    round(acgwin.*pupildilation.samplingRate),'coeff');
tlag = tlag./pupildilation.samplingRate;
%
figure
plot(tlag,pupilACG,'k','linewidth',2)
xlabel('t lag (s)')
ylabel('Pupil Diameter Autocorrelation')
%%
figure
plot(pupildilation.data,'k')
hold on
%plot(pupildilation.maxnormdata,'r')

%%
specparms.frange = [1 128];
specparms.nfreqs = 100;
specparms.ncyc = 5;
specparms.space = 'log';
specparms.numvarbins = 10;
specparms.winsize = 5;
specparms.noverlap = 4.9;
specparms.type = 'wavelet';
specparms.varnorm = 'percentile';
specparms.specnorm = 'mean';

figparms.figfolder = figfolder;
figparms.baseName = baseName;

%Loop the LFP, sorted by depth
numprobechannels = length(recparms.SpkGrps.Channels);
for ll = 1:numprobechannels
%for ll = 1:2
    channum = recparms.SpkGrps.Channels(ll);
    specparms_loop = specparms;
    specparms_loop.specnorm = 'mean';
    figparms_loop = figparms;
    figparms_loop.plotname = ['MeanNorm_Depth',num2str(ll),'_Channum',num2str(channum)];
    lfp = bz_GetLFP(channum,'basepath',basePath,'noPrompts',true);
    [specvarcorr_meannorm(ll),specbyvar_meannorm(ll)] = bz_LFPSpecToExternalVar(lfp,pupildilation,specparms_loop,figparms_loop);
    close all
    
    specparms_loop.specnorm = 'log';
    figparms_loop.plotname = ['_Depth',num2str(ll),'_Channum',num2str(channum)];
    [specvarcorr(ll),specbyvar(ll)] = bz_LFPSpecToExternalVar(lfp,pupildilation,specparms_loop,figparms_loop);
    close all
    
end

save(fullfile(figfolder,[baseName,'_SpecByPupil']))
%%
pupilLFPcorr = CollapseStruct(specvarcorr_meannorm,1);
LFPbypupil_meannorm = CollapseStruct(specbyvar_meannorm,3);
LFPbypupil = CollapseStruct(specbyvar,3);

%%
%exchanidx = [2,17,29]; %By depth (not channum)
exchanidx = [2,32,60]; %By depth (not channum)


pupilcolormap = makeColorMap([0 0 1],[0 0 0],[1 0 0],specparms.numvarbins);

corrmap = [makeColorMap([1 1 1],[0.4 0.4 0.8],[0 0 0.8]);
    makeColorMap([0 0 0.8],[0 0 0],[0.8 0 0])];
%corrmap = makeColorMap([0 0 0.8],[0 0 0],[0.8 0 0]);
figure
suptitle(baseName)
    subplot(5,5,[1 2 6 7])
        imagesc(log2(specvarcorr(1).freqs),[1 numprobechannels],pupilLFPcorr.corr)
        hold on
        plot(log2(specvarcorr(1).freqs(end)).*ones(size(exchanidx)),exchanidx,'w*')
        ylabel('Channel (by Depth)')
        xlabel('f (Hz)')
        colormap(gca,corrmap)
        ColorbarWithAxis([-0.3 0.1],'Correlation (rho)')
        LogScale('x',2)
        title('Pupil-Power Correlation')
        
    for pp = 1:length(specbyvar(1).varbins)
        
    subplot(4,5,pp+10)
    colormap(gca,'jet')
        imagesc(log2(specvarcorr(1).freqs),[1 numprobechannels],squeeze(LFPbypupil_meannorm.mean(:,pp,:))')
        LogScale('x',2)
        caxis([0.7 1.5])
        title(['Pupil Diameter: ',num2str(specbyvar(1).varbins(pp))])
        xlabel('f (Hz)')
        %colorbar
        %ColorbarWithAxis([0.6 1.4],'Power (mean^1)');
    end
    
    for ee = 1:length(exchanidx)
        subplot(6,5,4+(ee-1).*5)
            set(gca,'colororder',pupilcolormap)
            hold all
            set(gca,'colororder',pupilcolormap)
            plot(log2(specvarcorr(1).freqs),squeeze(LFPbypupil.mean(:,:,exchanidx(ee)))','linewidth',0.5)
            set(gca,'colororder',pupilcolormap)
            LogScale('x',2)
            xlim(log2(specvarcorr(1).freqs([1 end])));
            ylim([min(LFPbypupil.mean(:)) max(LFPbypupil.mean(:))])
            
            
        subplot(6,5,3+(ee-1).*5)
            %set(gca,'colororder',pupilcolormap)
            %hold all
           % set(gca,'colororder',pupilcolormap)
            plot(log2(specvarcorr(1).freqs),pupilLFPcorr.corr(exchanidx(ee),:))
            hold on
            plot(get(gca,'xlim'),[0 0],'r--')
            plot(log2(specvarcorr(1).freqs),pupilLFPcorr.corr(exchanidx(ee),:),'k','linewidth',2)
            LogScale('x',2)
           % set(gca,'colororder',pupilcolormap)
            xlim(log2(specvarcorr(1).freqs([1 end])));
           ylim([min(pupilLFPcorr.corr(:)) max(pupilLFPcorr.corr(:))])
           
           
        subplot(6,5,5+(ee-1).*5)
            colormap(gca,'jet')

            imagesc(log2(specvarcorr(1).freqs),specbyvar(1).varbins(pp),squeeze(LFPbypupil.mean(:,:,exchanidx(ee)))')
            %xlim(log2(specvarcorr(1).freqs([1 end])));
            axis xy
            caxis([0.7 1.4])
            %colorbar
            LogScale('x',2)
           caxis([min(LFPbypupil.mean(:)) max(LFPbypupil.mean(:))])
    end
    
    NiceSave('LFPPupilbyDepth',figfolder,baseName)
    
    %%
    figure
        for pp = 1:length(specbyvar(1).varbins)
        
    subplot(4,5,pp+10)
    colormap(gca,'jet')
        imagesc(log2(specvarcorr(1).freqs),[1 numprobechannels],squeeze(LFPbypupil.mean(:,pp,:))')
        LogScale('x',2)
       % caxis([0.4 1.6])
        title(['Pupil Diameter: ',num2str(specbyvar(1).varbins(pp))])
        xlabel('f (Hz)')
        %colorbar
        %ColorbarWithAxis([0.6 1.4],'Power (mean^1)');
    end
end    
% %%
% lfp = bz_GetLFP(recparms.SpkGrps.Channels);
% 
% %% Wavelet for Example channels
% 
% %exchanidx = [2,17,29]; %By depth (not channum)
% 
% %%
% for ee = 1:length(exchanidx)
%         [freqs,~,spec{ee}] = WaveSpec(single(lfp.data(:,exchanidx(ee))),...
%             specparms.frange,specparms.nfreqs,specparms.ncyc,...
%             1/lfp.samplingRate,specparms.space);
%     %spec{ee} = log(abs(spec{ee}));
%     spec{ee} = (abs(spec{ee}));
%     spec_meannorm{ee} = bsxfun(@(X,Y) X./Y,spec{ee},mean(spec{ee},2));
%     spec{ee} = log(spec{ee});
%     
% %         %Calculate the FFT spectrogram parameters - covert from s to sf
% %         winsize = specparms.winsize*lfp.samplingRate;
% %         noverlap = specparms.noverlap*lfp.samplingRate; %Might want to calaulte overlap based on pupil...?
% %         %Calculate the FFT spectrogram
% %         [FFTspec{ee},~,spectimestamps] = spectrogram(single(lfp.data(:,exchanidx(ee))),...
% %             winsize,noverlap,freqs,lfp.samplingRate);
% %         FFTspec{ee} = log(abs(FFTspec{ee}));
%     
% end
% 
% %%
% figure
% plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2)
% %%
% winsize = 50;
% timewin = [635 655];
% %timewin = [600 700];
% timewin = [330 350];
% timewin = [1445 1495];
% timewin = 770+[0 winsize];
% % figure
% % subplot(3,1,1)
% % bz_MultiLFPPlot(lfp,'timewin',timewin,'channels',lfp.channels(1:2:end-2))
% % subplot(6,1,3)
% % plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2)
% % xlim(timewin)
% % 
% % %%
% depthnames = {'Superficial','Intermediate','Deep'};
% %timewin = [1445 1495];
% figure
% for ee = 1:length(exchanidx)
% subplot(10,1,(ee-1).*2+[1 2])
% colormap('jet')
% imagesc(lfp.timestamps,log2(freqs),(spec{ee}))
% xlim(timewin)
% axis xy
% LogScale('y',2)
% ylabel({depthnames{ee},'f (Hz)'});
% %SpecColorRange(spec{2},[1.4 1.5]);
% %SpecColorRange(spec{2},[1.4 1.7]);
% SpecColorRange(spec{2},[0.8 1.4]);
% %caxis([-3 3])
% set(gca,'xtick',[])
% %box off
% end
% subplot(12,1,8:10)
%     bz_MultiLFPPlot(lfp,'timewin',timewin,'channels',lfp.channels(exchanidx))
%         
%         set(gca,'xtick',[])
%         box off
%         ylabel({'LFP','by Depth'})
%         xlabel('')
% subplot(12,1,11:12)
%     plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2)
%     xlim(timewin)
%     xlabel('t (s)')
%     box off
%     ylabel({'Pupil', 'Diameter'})
%   NiceSave('ExampleDilation_Spec',figfolder,baseName,'tiff')  
%     NiceSave('ExampleDilation_Spec',figfolder,baseName,'pdf')  
% 
% %% Sample zooms
% winsize=5;
% examplewins = [690+[0 winsize];712+[0 winsize];717+[0 winsize]]
% %examplewins = [1448 1453; 1469.5 1474.5; 1485 1490];
% for xx = 1:3
%     timewin = examplewins(xx,:);
% figure
% for ee = 1:length(exchanidx)
% %subplot(12,1,(ee-1).*2+[1 2])
% subplot(12,2,(ee-1).*4+[1 3])
% colormap('jet')
% imagesc(lfp.timestamps,log2(freqs),(spec{ee}))
% xlim(timewin)
% axis xy
% LogScale('y',2)
% 
% %SpecColorRange(spec{2},[1.9 1.9]);
% SpecColorRange(spec{2},[0.8 1.4]);
% %caxis([-3 3])
% set(gca,'xtick',[])
% %box off
% end
% subplot(12,2,[13:2:19])
%     bz_MultiLFPPlot(lfp,'timewin',timewin,'channels',lfp.channels(1:end-2))
%         
%     %    set(gca,'xtick',[])
%         xlabel('')
%         box off
% % subplot(12,1,12)
% %     plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2)
% %     xlim(timewin)
% %     xlabel('t (s)')
% %     box off
% 
% subplot(12,2,[21 23])
%     plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2)
%     xlim(timewin)
%     xlabel('t (s)')
%     box off
%     ylabel({'Pupil', 'Diameter'})
%   NiceSave(['ExampleDilation_Spec',num2str(xx)],figfolder,baseName,'tiff')    
%     NiceSave(['ExampleDilation_Spec',num2str(xx)],figfolder,baseName,'pdf')     
% end
%     
% %%    
% figure
% subplot(3,1,1)
% imagesc(lfp.timestamps,[1 30], lfp.data')
% hold on
% colorbar
% caxis([-0.1e4 0.1e4])
% %plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2)
% xlim(timewin)
%     NiceSave('ExampleDilation',figfolder,baseName)    
% 
% %%
% 
%     
% %end
