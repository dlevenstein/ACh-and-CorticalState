analysisfolder = '/Users/dlevenstein/Project Repos/ACh-and-CorticalState/AnalysisScripts/AnalysisFigs/UPDOWNandPupilAnalysis';
UPDOWNPupilAll = GetMatResults(analysisfolder,'UPDOWNandPupilAnalysis');
%genotype = {PupilWhiskAll.genotype};
%[genotypes,~,genotypeidx] = unique(genotype);
%%
UPDOWNpupil = CollapseStruct(UPDOWNPupilAll);

%%
animalname = extractBefore(UPDOWNpupil.name, '_');
genotype = extractBetween(UPDOWNpupil.name,'_','_');
genotype = genotype(:,:,1);
[genotypes,~,genotypeidx] = unique(genotype);
%genotypes = char(genotypes);
%%
pupcycleUPDOWN.ALL = CollapseStruct(UPDOWNpupil.pupcycleUPDOWN,3);
pupcycleUPDOWN.bins = UPDOWNpupil.pupcycleUPDOWN.bincenters;

pupphaseUD.ALL = CollapseStruct(UPDOWNpupil.pupphaseUD,1);
pupphaseUD.freqs = UPDOWNpupil.pupphaseUD.freqs;

%%
for gg = 1:length(genotypes)
    pupcycleUPDOWN.pDOWN.(genotypes{gg}) = nanmean(pupcycleUPDOWN.ALL.pDOWN(:,:,genotypeidx==gg),3);
    pupcycleUPDOWN.logUPdur.(genotypes{gg}) = nanmean(pupcycleUPDOWN.ALL.logUPdur(:,:,genotypeidx==gg),3);
    
    pupphaseUD.phaseDN = nanmean(pupphaseUD.ALL.phaseDN.mag,1);
    pupphaseUD.phaseUPdur = nanmean(pupphaseUD.ALL.phaseUPdur.corr,1);
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




%% UP/DOWN and dpdt
figure
for gg = 1:length(genotypes)
    subplot(2,2,gg)
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
        
        
    subplot(2,2,gg+2)
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


%% Phase Coupling

