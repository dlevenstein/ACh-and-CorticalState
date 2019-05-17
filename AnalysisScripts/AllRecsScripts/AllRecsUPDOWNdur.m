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

groups = {KOgeneotypes,WTgeneotypes};
groupnames = {'KOs','WTctr'};
%%
for gg = 1:6
    ConditionalUPDOWN.(genotypes{gg}) = bz_CollapseStruct(UPDOWNDur.ConditionalUPDOWN(genotypeidx==gg),3,'mean',true);
    ReturnHist.(genotypes{gg}) = bz_CollapseStruct(UPDOWNDur.ReturnHist(genotypeidx==gg),3,'mean',true);
    PSShist.(genotypes{gg}) = bz_CollapseStruct(UPDOWNDur.PSShist(genotypeidx==gg),3,'mean',true);
    PSShist.std.(genotypes{gg}) = bz_CollapseStruct(UPDOWNDur.PSShist(genotypeidx==gg),3,'std',true);
end
ConditionalUPDOWN.AllWT = bz_CollapseStruct(UPDOWNDur.ConditionalUPDOWN(strcmp(WTKOtype,'WT')),3,'mean',true);
ReturnHist.AllWT = bz_CollapseStruct(UPDOWNDur.ReturnHist(strcmp(WTKOtype,'WT')),3,'mean',true);
PSShist.AllWT = bz_CollapseStruct(UPDOWNDur.PSShist(strcmp(WTKOtype,'WT')),3,'mean',true);
PSShist.std.AllWT = bz_CollapseStruct(UPDOWNDur.PSShist(strcmp(WTKOtype,'WT')),3,'std',true);
%%
UPDOWN = {'UP','DOWN'};
UDcolor = {'r','b'};
WHNWH = {'wh','nwh'};
HILO = {'hipup','lopup'};
UDcmap = {makeColorMap([1 1 1],[0.8 0 0]),makeColorMap([1 1 1],[0 0 0.8])};

cosx = linspace(-pi,3*pi,100);
cospamp = [0.6 0.05];


for ff = 1:2
for uu = 1:2
figure
for gi=1:length(groups{ff})
    gg = groups{ff}(gi);
    for ww = 1:2; for pp = 1:2
    %subplot(4,4,4*(gi-1)+pp+2*(ww-1))
    subplot(4,4,(ww-1)*8+(pp-1)*4+gi)
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
         %colorbar
         
         caxis([0 0.08])
         if uu==2
            %xlabel('Pup Phase')
         end
        if ww == 1 & pp==1
            %ylabel({genotypes{gg},[(UPDOWN{uu}),' Dur (s)']})
            title(genotypes{gg})
        end
        if ww == 2 & pp==2
            %ylabel({genotypes{gg},[(UPDOWN{uu}),' Dur (s)']})
            xlabel('Pup Phase')
        end
        if gi==1
            %title([(WHNWH{ww}),' ',(HILO{pp})]
            ylabel({(WHNWH{ww}),[(UPDOWN{uu}),' Dur (s)']})
        end
    
end;end;end

NiceSave([(UPDOWN{uu}),'DurbyPupCycle'],analysisfolder,groupnames{ff})

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