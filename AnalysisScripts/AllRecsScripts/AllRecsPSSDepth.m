analysisfolder = '/Users/dlevenstein/Project Repos/ACh-and-CorticalState/AnalysisScripts/AnalysisFigs/LFPSlopebyDepthAnalysis';
SlopeDepthAll = GetMatResults(analysisfolder,'LFPSlopebyDepthAnalysis');
%genotype = {PupilWhiskAll.genotype};
%[genotypes,~,genotypeidx] = unique(genotype);
%%
SlopeDepth = bz_CollapseStruct(SlopeDepthAll);

%%
animalname = extractBefore(SlopeDepth.name, '_');
genotype = extractBetween(SlopeDepth.name,'_','.');
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

groups = {KOgeneotypes,WTgeneotypes};
groupnames = {'KOs','WTctr'};
%%
for gg = 1:6
    PSScomponents.(genotypes{gg}) = bz_CollapseStruct(SlopeDepth.PSScomponents(genotypeidx==gg),3,'mean',true);
    PSSdepth.(genotypes{gg}) = bz_CollapseStruct(SlopeDepth.PSSdepth(genotypeidx==gg),3,'mean',true);
    PSSdist.(genotypes{gg}) = bz_CollapseStruct(SlopeDepth.PSSdist(genotypeidx==gg),3,'mean',true);
end

PSScomponents.AllWT = bz_CollapseStruct(SlopeDepth.PSScomponents(strcmp(WTKOtype,'WT')),3,'mean',true);
PSSdepth.AllWT = bz_CollapseStruct(SlopeDepth.PSSdepth(strcmp(WTKOtype,'WT')),3,'mean',true);
PSSdist.AllWT = bz_CollapseStruct(SlopeDepth.PSSdist(strcmp(WTKOtype,'WT')),3,'mean',true);

%%

WHNWH = {'Wh','NWh'};
HILO = {'lopup','hipup'};
LAYERS = {'L1','L23','L4','L5a','L5b6','WM'};
depthinfo.boundaries = [0 0.1 0.35 0.5 1];
%%
cosx = linspace(-pi,pi,100);
cospamp = [0.025 0.225];
colorrange = [-2.4 -1.2];

for ff = 1:2
figure

for gi=1:length(groups{ff})
    gg = groups{ff}(gi);
for ww = 1:2
    subplot(5,4,(ww-1)*4+gi)
    
        for pp = 1:2
        imagesc( PSSdepth.(genotypes{gg}).(HILO{pp}).(WHNWH{ww}).varbins+2*pi*(pp-1),...
            PSSdepth.(genotypes{gg}).depth,...
            PSSdepth.(genotypes{gg}).(HILO{pp}).(WHNWH{ww}).mean_interp)
        hold on; axis xy; box off
        plot([-pi 3*pi],-depthinfo.boundaries'*[1 1],'w')
        plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp)-1,'k')
        end   
        
        
        %ColorbarWithAxis(colorrange,'Mean PSS')
        clim(colorrange)
        %colorbar
        xlim([-pi 3*pi])
        plot(pi*[1 1],get(gca,'ylim'),'k--')
        bz_piTickLabel('x')
        
        if ww == 2
            xlabel('Pupil Phase');
        end
        
        if gi==1
            ylabel({WHNWH{ww},'Depth'})
        end

        if ww == 1
            title(genotypes{gg})
        end
        
    subplot(5,4,(ww-1)*4+gi+8)
        imagesc( PSSdepth.(genotypes{gg}).pup.(WHNWH{ww}).varbins,...
             PSSdepth.(genotypes{gg}).depth,...
            PSSdepth.(genotypes{gg}).pup.(WHNWH{ww}).mean_interp)
        hold on; axis xy; box off
        %ColorbarWithAxis(colorrange,'Mean PSS')
        clim(colorrange)
        plot([-pi 3*pi],-depthinfo.boundaries'*[1 1],'w')
        if ww == 2
            xlabel('Pupil Size');
        end
        
        if gi==1
            ylabel({WHNWH{ww},'Depth'})
        end
end

subplot(5,4,gi+16)
        imagesc( PSSdepth.(genotypes{gg}).WhOn.all.varbins,...
            PSSdepth.(genotypes{gg}).depth,...
            PSSdepth.(genotypes{gg}).WhOn.all.mean_interp)
        hold on; axis xy; box off
        plot(PSSdepth.(genotypes{gg}).WhOn.all.varbins([1 end]),-depthinfo.boundaries'*[1 1],'w')
        plot([0 0],[-1 0],'w')
        %ColorbarWithAxis(colorrange,'Mean PSS')
        clim(colorrange)
        xlabel('t (s - relative to WhOn');ylabel('Depth')



end
NiceSave('DepthPSSandBeh',analysisfolder,groupnames{ff},'figtype','pdf')
end

%%
for ff = 1:2
    
figure
for gi=1:length(groups{ff})
    gg = groups{ff}(gi);
    subplot(4,4,gi)
    imagesc(PSScomponents.(genotypes{gg}).depth,PSScomponents.(genotypes{gg}).depth,...
        PSScomponents.(genotypes{gg}).corr)
    hold on
    plot(PSScomponents.(genotypes{gg}).depth([1 end]),-depthinfo.boundaries'*[1 1],'w')
    plot(-depthinfo.boundaries'*[1 1],PSScomponents.(genotypes{gg}).depth([1 end]),'w')

    colorbar
    clim([0.6 1])
    
    
    
    axis xy
    xlabel('Depth');ylabel('Depth')
    title('PSS Corr')
    
subplot(3,3,8)
    plot(log10(PSScomponents.(genotypes{gg}).EV),'o-')
    hold on
    xlim([1 6]);ylim([-1 2])
    xlabel('PC');ylabel('% EV')
    LogScale('y',10)
    %legend()

subplot(3,4,4+gi)
    plot(PSScomponents.(genotypes{gg}).depth,PSScomponents.(genotypes{gg}).PCAcoeff(:,1:3))
    hold on
    plot(PSScomponents.(genotypes{gg}).depth([1 end]),[0 0],'k--')
    legend('PC1','PC2','PC3','location','southoutside')
    xlabel('Deptah');ylabel('Weight')
    ylim([-0.25 0.4])
end

NiceSave('PSSDepthComponents',analysisfolder,groupnames{ff},'figtype','pdf')
end


%%

%%
cosx = linspace(-pi,pi,100);
cospamp = [0.05 0.5];


for ff = 1:2
    for ww = 1:2
figure
for gi=1:length(groups{ff})
    gg = groups{ff}(gi);
    
for ll = 1:length(LAYERS)

    subplot(6,4,(ll-1)*4+gi)
        for pp = 1:2
        imagesc( PSSdist.(genotypes{gg}).(LAYERS{ll}).(HILO{pp}).(WHNWH{ww}).Xbins+2*pi*(pp-1),...
            PSSdist.(genotypes{gg}).(LAYERS{ll}).(HILO{pp}).(WHNWH{ww}).Ybins,...
            PSSdist.(genotypes{gg}).(LAYERS{ll}).(HILO{pp}).(WHNWH{ww}).pYX')
        hold on; axis xy; box off

        plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp)-3,'k')
        end   
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        xlim([-pi 3*pi])
        plot(pi*[1 1],get(gca,'ylim'),'k--')
        bz_piTickLabel('x')
        %title(genotypes{gg})
    crameri bilbao
        if ll == 6
        xlabel('Pupil Phase');
        end
        if gi == 1
            ylabel({LAYERS{ll},'PSS'})
        end
        if ll == 1
        title({WHNWH{ww},genotypes{gg}})
        end
end
end

NiceSave(['PSSDistPupCycle_',WHNWH{ww}],analysisfolder,groupnames{ff},'figtype','pdf')

    end
end

%%
for ff = 1:2
    for ww = 1:2
figure
for gi=1:length(groups{ff})
    gg = groups{ff}(gi);
    
for ll = 1:length(LAYERS)
    subplot(6,4,(ll-1)*4+gi)
        imagesc( PSSdist.(genotypes{gg}).(LAYERS{ll}).pup.(WHNWH{ww}).Xbins,...
            PSSdist.(genotypes{gg}).(LAYERS{ll}).(HILO{pp}).(WHNWH{ww}).Ybins,...
            PSSdist.(genotypes{gg}).(LAYERS{ll}).pup.(WHNWH{ww}).pYX')
        hold on; axis xy; box off
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        if ll == 6
        xlabel('Pupil Size');
        end
        if gi == 1
            ylabel({LAYERS{ll},'PSS'})
        end
        if ll == 1
        title({WHNWH{ww},genotypes{gg}})
        end
   crameri bilbao
      
end
end

NiceSave(['PSSDistPupSize_',WHNWH{ww}],analysisfolder,groupnames{ff},'figtype','pdf')

    end
end
%%
        
for ff = 1:2

figure
for gi=1:length(groups{ff})
    gg = groups{ff}(gi);
    
for ll = 1:length(LAYERS)
    subplot(6,4,(ll-1)*4+gi)
        imagesc( PSSdist.(genotypes{gg}).(LAYERS{ll}).EMG.Xbins,...
            PSSdist.(genotypes{gg}).(LAYERS{ll}).EMG.Ybins,...
            PSSdist.(genotypes{gg}).(LAYERS{ll}).EMG.pYX')
        hold on; axis xy; box off
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        %colorbar
        xlim([-2 1])
        LogScale('x',10)
        if ll == 6
        xlabel('EMG');
        end
        if gi == 1
            ylabel({LAYERS{ll},'PSS'})
        end
        if ll == 1
        title(genotypes{gg})
        end
        
    crameri bilbao
end


end
NiceSave('PSSDistEMG',analysisfolder,groupnames{ff},'figtype','pdf')

%NiceSave('LayerPSSandBehavior',analysisfolder,groupnames{ff})
end