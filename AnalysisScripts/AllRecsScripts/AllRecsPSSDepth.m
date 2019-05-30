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
end

PSScomponents.AllWT = bz_CollapseStruct(SlopeDepth.PSScomponents(strcmp(WTKOtype,'WT')),3,'mean',true);
PSSdepth.AllWT = bz_CollapseStruct(SlopeDepth.PSSdepth(strcmp(WTKOtype,'WT')),3,'mean',true);

%%

WHNWH = {'Wh','NWh'};
HILO = {'lopup','hipup'};

%%
cosx = linspace(-pi,pi,100);
cospamp = [0.025 0.3];
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
        plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp)-1,'k')
        end   
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        clim(colorrange)
        xlim([-pi 3*pi])
        
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
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        clim(colorrange)
        if ww == 2
            xlabel('Pupil Size');
        end
        
        if gi==1
            ylabel({WHNWH{ww},'Depth'})
        end
end

subplot(5,4,gi+16)
        imagesc( PSSdepth.(genotypes{gg}).whOn.varbins,...
            PSSdepth.(genotypes{gg}).depth,...
            PSSdepth.(genotypes{gg}).whOn.mean_interp)
        hold on; axis xy; box off
        plot([0 0],[-1 0],'w')
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        clim(colorrange)
        xlabel('t (s - relative to WhOn');ylabel('Depth')



end
NiceSave('DepthPSSandBeh',analysisfolder,groupnames{ff})
end

%%
for ff = 1:2
figure
for gi=1:length(groups{ff})
    gg = groups{ff}(gi);
    subplot(4,4,gi)
    imagesc(PSScomponents.(genotypes{gg}).depth,PSScomponents.(genotypes{gg}).depth,...
        PSScomponents.(genotypes{gg}).corr)
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

NiceSave('PSSDepthComponents',analysisfolder,groupnames{ff})
end

%%


for ff = 1:2
    figure
    for uu = 1:2; for gi=1:length(groups{ff})
    gg = groups{ff}(gi);
    for ww = 1:2 
    subplot(4,4,(ww-1)*4+gi+(uu-1)*8)
        colormap(gca,UDcmap{uu})
    
        for pp = 1:2
        imagesc(ConditionalUPDOWN.(genotypes{gg}).(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).Xbins+2*pi*(pp-1),...
            ConditionalUPDOWN.(genotypes{gg}).(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).Ybins,...
            ConditionalUPDOWN.(genotypes{gg}).(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).pYX')
        hold on; axis xy; box off
        plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp)-1.5,'k')
        end
        LogScale('y',10)
         xlim([-pi 3*pi])
         %colorbar
         
         caxis([0 0.08])
        if ww == 1 & uu ==1
            title(genotypes{gg})
        end
        if ww == 2
            xlabel('Pup Phase')
        end
        if gi==1
            ylabel({(WHNWH{ww}),[(UPDOWN{uu}),' Dur (s)']})
        end
    
end;end;end

NiceSave('UPDOWNDurbyPupCycle',analysisfolder,groupnames{ff})

end
end


%% Figure: return maps
for ff = 1:2

figure
for gi=1:length(groups{ff})
    gg = groups{ff}(gi);
    for ww = 1:2; for uu = 1:2
    %subplot(4,4,9+(uu-1)*4+(ww-1)*2)
    subplot(4,4,(ww-1)*8+(uu-1)*4+gi)
    colormap(gca,UDcmap{uu})
    imagesc(ConditionalUPDOWN.(genotypes{gg}).(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).Ybins,...
        ConditionalUPDOWN.(genotypes{gg}).(UPDOWN{uu}).(HILO{pp}).(WHNWH{ww}).Ybins,...
        ReturnHist.(genotypes{gg}).(UPDOWN{uu}).allpup.(WHNWH{ww}))
    axis xy
    LogScale('xy',10)
    if uu==2 & ww == 2
        xlabel('Dur_n (s)')
    end
    if gi == 1
        ylabel({(WHNWH{ww}),'Dur_n_+_1 (s)'})
    end
        if ww == 1 & uu==1
            %ylabel({genotypes{gg},[(UPDOWN{uu}),' Dur (s)']})
            title(genotypes{gg})
        end
    
            end;end
end

NiceSave(['Return Maps'],analysisfolder,groupnames{ff})

end


%%
for ff = 1:2

figure
for gi=1:length(groups{ff})
    gg = groups{ff}(gi);
    for ww = 1:2; for uu = 1:2
    subplot(4,4,(ww-1)*8+(uu-1)*4+gi)
    colormap(gca,UDcmap{uu})
    imagesc(ConditionalUPDOWN.(genotypes{gg}).(UPDOWN{uu}).PSS.(WHNWH{ww}).Xbins,...
        ConditionalUPDOWN.(genotypes{gg}).(UPDOWN{uu}).PSS.(WHNWH{ww}).Ybins,...
        ConditionalUPDOWN.(genotypes{gg}).(UPDOWN{uu}).PSS.(WHNWH{ww}).pYX')
    hold on
%     plot(ConditionalUPDOWN.(genotypes{gg}).(UPDOWN{uu}).PSS.(WHNWH{ww}).Xbins([1 end]),...
%         log10(PSS.winsize*[1 1]),'k--')
    plot(PSShist.(genotypes{gg}).bins,PSShist.(genotypes{gg}).allpup.(WHNWH{ww})*8+...
        ConditionalUPDOWN.(genotypes{gg}).(UPDOWN{uu}).PSS.(WHNWH{ww}).Ybins(1)-0.5,'k','linewidth',1)
    axis xy; box off;axis tight
   % xlim([-1.5 1])
    LogScale('y',10)
    %colorbar
    caxis([0 0.1])
    if uu==2 & ww == 2
        xlabel('PSS')
    end
    if gi == 1
        ylabel({(WHNWH{ww}),'Dur (s)'})
    end
        if ww == 1 & uu==1
            %ylabel({genotypes{gg},[(UPDOWN{uu}),' Dur (s)']})
            title(genotypes{gg})
        end
       
            end;end
end

NiceSave('DurByPSS',analysisfolder,groupnames{ff})

end


%% UD by pupil size

for ff = 1:2

figure
for gi=1:length(groups{ff})
    gg = groups{ff}(gi);
    for ww = 1:2; for uu = 1:2
    subplot(4,4,(ww-1)*8+(uu-1)*4+gi)
    colormap(gca,UDcmap{uu})
    imagesc(ConditionalUPDOWN.(genotypes{gg}).(UPDOWN{uu}).pupmag.(WHNWH{ww}).Xbins,...
        ConditionalUPDOWN.(genotypes{gg}).(UPDOWN{uu}).pupmag.(WHNWH{ww}).Ybins,...
        ConditionalUPDOWN.(genotypes{gg}).(UPDOWN{uu}).pupmag.(WHNWH{ww}).pYX')
    axis xy
   % xlim([-1.5 1])
    LogScale('y',10)
    caxis([0 0.08])
    if uu==2 & ww == 2
        xlabel('Pupil (log)')
    end
    if gi == 1
        ylabel({(WHNWH{ww}),'Dur (s)'})
    end
        if ww == 1 & uu==1
            %ylabel({genotypes{gg},[(UPDOWN{uu}),' Dur (s)']})
            title(genotypes{gg})
        end
       
            end;end
end

NiceSave('DurByPupSize',analysisfolder,groupnames{ff})

end