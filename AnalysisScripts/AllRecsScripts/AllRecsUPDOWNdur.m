analysisfolder = '/Users/dlevenstein/Project Repos/ACh-and-CorticalState/AnalysisScripts/AnalysisFigs/UPDOWNDurationAnalysis';
UPDOWNDurAll = GetMatResults(analysisfolder,'UPDOWNDurationAnalysis');
%genotype = {PupilWhiskAll.genotype};
%[genotypes,~,genotypeidx] = unique(genotype);
%%
UPDOWNDur = bz_CollapseStruct(UPDOWNDurAll);

%%
animalname = extractBefore(UPDOWNDur.name, '_');
genotype = extractBetween(UPDOWNDur.name,'_','.');
%genotype = extractBetween(UPDOWNDur.name,'_','_');
%genotype = genotype(:,:,1);
WTKOtype = extractBefore(genotype,'_');
KOtype = extractAfter(genotype,'_');
[genotypes,~,genotypeidx] = unique(genotype);
%genotypes = char(genotypes);
genotypes{7} = 'AllWT';
KOgeneotypes = find(strncmpi(genotypes,'KO',2));
WTgeneotypes = find(strncmpi(genotypes,'WT',2));
KOgeneotypes = [7,KOgeneotypes]; %add allWT to KO comparison
%%
for gg = 1:6
    ConditionalUPDOWN.(genotypes{gg}) = bz_CollapseStruct(UPDOWNDur.ConditionalUPDOWN(genotypeidx==gg),3,'mean',true);
end
ConditionalUPDOWN.AllWT = bz_CollapseStruct(UPDOWNDur.ConditionalUPDOWN(strcmp(WTKOtype,'WT')),3,'mean',true);
%%
UPDOWN = {'UP','DOWN'};
UDcolor = {'r','b'};
WHNWH = {'wh','nwh'};
HILO = {'hipup','lopup'};
UDcmap = {makeColorMap([1 1 1],[0.8 0 0]),makeColorMap([1 1 1],[0 0 0.8])};

cosx = linspace(-pi,3*pi,100);
cospamp = [0.6 0.05];

groups = {KOgeneotypes,WTgeneotypes};
for ff = 1:2
for uu = 1:2
figure
for gi=1:length(groups{ff})
    gg = groups{ff}(gi);
    for ww = 1:2; for pp = 1:2
    subplot(4,4,4*(gi-1)+pp+2*(ww-1))
    colormap(gca,UDcmap{uu})
        imagesc(ConditionalUPDOWN.(genotypes{gg}).(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).Xbins,...
            ConditionalUPDOWN.(genotypes{gg}).(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).Ybins,...
            ConditionalUPDOWN.(genotypes{gg}).(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).pYX')
        hold on; axis xy; box off
        imagesc(ConditionalUPDOWN.(genotypes{gg}).(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).Xbins+2*pi,...
            ConditionalUPDOWN.(genotypes{gg}).(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).Ybins,...
            ConditionalUPDOWN.(genotypes{gg}).(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).pYX')
        plot(cosx,(cos(cosx)+1).*cospamp(pp)-1.5,'k')
        LogScale('y',10)
         xlim([-pi 3*pi])
         caxis([0 0.15])
         if uu==2
            xlabel('Pup Phase')
         end
        if ww == 1 &pp==1
            ylabel({genotypes{gg},[(UPDOWN{uu}),' Dur (s)']})
        end
        if gi==1
            title([(WHNWH{ww}),' ',(HILO{pp})])
        end
    
end;end;end
end
end




%%
pupcycleUPDOWN.ALL = bz_CollapseStruct(UPDOWNDur.pupcycleUPDOWN,3);
pupcycleUPDOWN.bins = UPDOWNDur.pupcycleUPDOWN.bincenters;

pupphaseUD.ALL = bz_CollapseStruct(UPDOWNDur.pupphaseUD,1,'justcat',true);
pupphaseUD.freqs = UPDOWNDur.pupphaseUD.freqs;

islowhist.ALL = bz_CollapseStruct(UPDOWNDur.islowhist,3);
islowhist.phasebins = UPDOWNDur.islowhist.phasebins;
islowhist.ampbins = UPDOWNDur.islowhist.ampbins;

pupbyislow.ALL = bz_CollapseStruct(UPDOWNDur.pupbyislow,3);
pupbyislow.phasebins = UPDOWNDur.pupbyislow.phasebins;
pupbyislow.ampbins = UPDOWNDur.pupbyislow.ampbins;

%%
for gg = 1:length(genotypes)
    pupcycleUPDOWN.pDOWN.(genotypes{gg}) = nanmean(pupcycleUPDOWN.ALL.pDOWN(:,:,genotypeidx==gg),3);
    pupcycleUPDOWN.logUPdur.(genotypes{gg}) = nanmean(pupcycleUPDOWN.ALL.logUPdur(:,:,genotypeidx==gg),3);
    
    pupphaseUD.phaseDN.(genotypes{gg}) = nanmean(pupphaseUD.ALL.phaseDN.mag(genotypeidx==gg,:),1);
    pupphaseUD.phaseUPdur.(genotypes{gg}) = nanmean(pupphaseUD.ALL.phaseUPdur.corr(genotypeidx==gg,:),1);
    
    islowhist.pDOWNPUP.(genotypes{gg}) = nanmean(islowhist.ALL.pDOWNPUP(:,:,genotypeidx==gg),3);
    islowhist.pDOWNPUP_marg.(genotypes{gg}) = nanmean(islowhist.ALL.pDOWNPUP_marg(:,:,genotypeidx==gg),3);
    
    pupbyislow.meanpupbyphase.(genotypes{gg}) = nanmean(pupbyislow.ALL.meanpupbyphase(:,:,genotypeidx==gg),3);
end

%% All Recordings: UP/DOWN and dpdt
figure
for gg = 1:length(genotypes)
    recs = find(genotypeidx==gg);
    for rr = 1:length(recs)
        subplot(4,3,rr+(gg-1)*6)
            h = imagesc(pupcycleUPDOWN.bins,pupcycleUPDOWN.bins,pupcycleUPDOWN.ALL.pDOWN(:,:,recs(rr))');
            set(h,'AlphaData',~isnan(pupcycleUPDOWN.ALL.pDOWN(:,:,recs(rr))'))
            colorbar
            hold on
            plot(pupcycleUPDOWN.bins([1 end]),[0 0],'k--')
            axis xy
            caxis([0 3])
            title(['DN rate: ',animalname{recs(rr)}])
            %xlim([-0.6 0.5])
            LogScale('x',10)
            xlabel('Pupil Area');ylabel('dp/dt')
    end
end
NiceSave('PupandDOWN_all',analysisfolder,'_')

figure
for gg = 1:length(genotypes)
    recs = find(genotypeidx==gg);
    for rr = 1:length(recs)
        subplot(4,3,rr+(gg-1)*6)
            h = imagesc(pupcycleUPDOWN.bins,pupcycleUPDOWN.bins,pupcycleUPDOWN.ALL.logUPdur(:,:,recs(rr))');
            set(h,'AlphaData',~isnan(pupcycleUPDOWN.ALL.logUPdur(:,:,recs(rr))'))
            colorbar
            hold on
            plot(pupcycleUPDOWN.bins([1 end]),[0 0],'k--')
            axis xy
            title(['UP dur: ',animalname{recs(rr)}])
            LogScale('x',10)
            caxis([-1 1])
            LogScale('c',10)
            LogScale('x',10)
            xlabel('Pupil Area');ylabel('dp/dt')
    end
end

NiceSave('PupandUP_all',analysisfolder,'_')


%% UP/DOWN and dpdt
figure
for gg = 1:length(genotypes)
    subplot(3,3,gg)
        h = imagesc(pupcycleUPDOWN.bins,pupcycleUPDOWN.bins,pupcycleUPDOWN.pDOWN.(genotypes{gg})');
        set(h,'AlphaData',~isnan(pupcycleUPDOWN.pDOWN.(genotypes{gg})'))
        colorbar
        hold on
        plot(pupcycleUPDOWN.bins([1 end]),[0 0],'k--')
        axis xy
        caxis([0 3])
        title(['DOWN rate: ',genotypes{gg}])
        %xlim([-0.6 0.5])
        LogScale('x',10)
        xlabel('Pupil Area');ylabel('dp/dt')
        
        
    subplot(3,3,gg+3)
        h = imagesc(pupcycleUPDOWN.bins,pupcycleUPDOWN.bins,pupcycleUPDOWN.logUPdur.(genotypes{gg})');
        set(h,'AlphaData',~isnan(pupcycleUPDOWN.logUPdur.(genotypes{gg})'))
        colorbar
        hold on
        plot(pupcycleUPDOWN.bins([1 end]),[0 0],'k--')
        axis xy
        %caxis([0 3])
        title(['mean UP duration: ',genotypes{gg}])
        LogScale('x',10)
        caxis([-1 1])
        LogScale('c',10)
        xlabel('Pupil Area');ylabel('dp/dt')
end


NiceSave('PupandUPDOWN_mean',analysisfolder,'_')

%% Phase Coupling: DOWN incidence

figure
    subplot(2,2,1)
        plot(log10(pupphaseUD.freqs),pupphaseUD.phaseDN.WT,'k','linewidth',2)
        hold on
        plot(log10(pupphaseUD.freqs),pupphaseUD.phaseDN.EMX,'r','linewidth',2)
        LogScale('x',10)
        xlabel('Pupil Freq (Hz)');ylabel('DOWN-Pupil Coupling')
    

        
    for gg = 1:length(genotypes)
        subplot(4,2,2.*gg)

        imagesc(islowhist.phasebins,islowhist.ampbins,islowhist.pDOWNPUP.(genotypes{gg})')
        hold on
        imagesc(islowhist.phasebins+2*pi,islowhist.ampbins,islowhist.pDOWNPUP.(genotypes{gg})')
        axis xy
        xlim(pi.*[-1 3])
        clim([0.3 2])
        colorbar('location','east')
        ylabel('iSlow Amp')
        title([genotypes{gg},': DN rate'])
    end
    
        subplot(4,2,6)
            bar(islowhist.phasebins,islowhist.pDOWNPUP_marg.EMX,'r')
            hold on
            bar(islowhist.phasebins+2*pi,islowhist.pDOWNPUP_marg.EMX,'r')
            bar(islowhist.phasebins,islowhist.pDOWNPUP_marg.WT,'k')
            bar(islowhist.phasebins+2*pi,islowhist.pDOWNPUP_marg.WT,'k')
            box off
            axis tight
            %colorbar
            ylabel('Mean DOWN rate')
            
        subplot(4,2,8)
            plot(pupbyislow.phasebins,pupbyislow.meanpupbyphase.WT,'k')
            hold on
            plot(pupbyislow.phasebins+2*pi,pupbyislow.meanpupbyphase.WT,'k')
            plot(pupbyislow.phasebins,pupbyislow.meanpupbyphase.EMX,'r')
            plot(pupbyislow.phasebins+2*pi,pupbyislow.meanpupbyphase.EMX,'r')
            axis tight
            xlabel('Pupil: iSlow Phase (0.005-0.05Hz)')
            ylabel('Mean Pupil Area')

NiceSave('PupPhaseandDOWN_mean',analysisfolder,'_')
 %% all recordings...
 
 figure
 for gg = 1:length(genotypes)
    recs = find(genotypeidx==gg);
    for rr = 1:length(recs)
        subplot(4,3,rr+(gg-1)*6)
            imagesc(islowhist.phasebins,islowhist.ampbins,islowhist.ALL.pDOWNPUP(:,:,recs(rr))')
            hold on
            imagesc(islowhist.phasebins+2*pi,islowhist.ampbins,islowhist.ALL.pDOWNPUP(:,:,recs(rr))')
            axis xy
            xlim(pi.*[-1 3])
            clim([0 2.5])
            colorbar('location','east')
            ylabel('iSlow Amp')
            title(['DN rate: ',animalname{recs(rr)}])
    end
 end
 NiceSave('DNbypupislow_all',analysisfolder,'_')
 
 figure
  for gg = 1:length(genotypes)
    recs = find(genotypeidx==gg);
    for rr = 1:length(recs)
        subplot(4,3,rr+(gg-1)*6)
            plot(log10(pupphaseUD.freqs),pupphaseUD.ALL.phaseDN.mag(recs(rr),:),'k','linewidth',1)
            hold on
            plot(log10(pupphaseUD.freqs(pupphaseUD.ALL.phaseDN.sig(recs(rr),:)>3)),...
                pupphaseUD.ALL.phaseDN.mag(recs(rr),(pupphaseUD.ALL.phaseDN.sig(recs(rr),:)>3)),'k.','markersize',10)
            LogScale('x',10)
            xlabel('Pupil Freq (Hz)');ylabel('DOWN-Pupil Coupling')
            ylim([0 0.25])
            title(['Pup-DN coupling: ',animalname{recs(rr)}])
    end
  end
  NiceSave('DNbypupphase_all',analysisfolder,'_')
%% Phase Coupling: UP duration
figure
    subplot(2,2,2)
        plot(log10(pupphaseUD.freqs),pupphaseUD.phaseUPdur.WT,'k','linewidth',2)
        hold on
        plot(log10(pupphaseUD.freqs),pupphaseUD.phaseUPdur.EMX,'r','linewidth',2)
        LogScale('x',10)
        xlabel('Pupil Freq (Hz)');ylabel('UPdur-Pupil Corr')
        
NiceSave('PupPhaseandUP_mean',analysisfolder,'_')