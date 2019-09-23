analysisfolder = '/Users/dlevenstein/Project Repos/ACh-and-CorticalState/AnalysisScripts/AnalysisFigs/LFPWavSpecbyDepthAnalysis';
%analysisfolder = '/home/dlevenstein/ProjectRepos/ACh-and-CorticalState/AnalysisScripts/AnalysisFigs/LFPWavSpecbyDepthAnalysis';
SpecDepthAll = GetMatResults(analysisfolder,'LFPWavSpecbyDepthAnalysis');
%genotype = {PupilWhiskAll.genotype};
%[genotypes,~,genotypeidx] = unique(genotype);
%%
SpecDepth = bz_CollapseStruct(SpecDepthAll);

%%
animalname = extractBefore(SpecDepth.name, '_');
genotype = extractBetween(SpecDepth.name,'_','.');
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
    SPECdepth.(genotypes{gg}) = bz_CollapseStruct(SpecDepth.SPECdepth(genotypeidx==gg),3,'mean',true);
    OSCdepth.(genotypes{gg}) = bz_CollapseStruct(SpecDepth.OSCdepth(genotypeidx==gg),3,'mean',true);
    SPECdepth_std.(genotypes{gg}) = bz_CollapseStruct(SpecDepth.SPECdepth(genotypeidx==gg),3,'std',true);

    %speccorr.(genotypes{gg}) = bz_CollapseStruct(SpecDepth.speccorr(genotypeidx==gg),3,'mean',true);
end

SPECdepth.AllWT = bz_CollapseStruct(SpecDepth.SPECdepth(strcmp(WTKOtype,'WT')),3,'mean',true);
SPECdepth_std.AllWT = bz_CollapseStruct(SpecDepth.SPECdepth(strcmp(WTKOtype,'WT')),3,'std',true);

OSCdepth.AllWT = bz_CollapseStruct(SpecDepth.OSCdepth(strcmp(WTKOtype,'WT')),3,'mean',true);
%speccorr.AllWT = bz_CollapseStruct(SpecDepth.speccorr(strcmp(WTKOtype,'WT')),3,'mean',true);

%%

ONOFF = {'WhOn','WhOFF'};
WHNWH = {'Wh','NWh'};
HILO = {'lopup','hipup'};
LONGSHORT = {'long','short'};
LAYERS = {'L1','L23','L4','L5a','L5b6','L6'};
depthinfo.boundaries = [0 0.1 0.35 0.5 0.6 0.9 1];



%%
if exist('speccorr','var')
%band = [0.5 2];
for ff = 1:2
figure
for gi=1:length(groups{ff})
    gg = groups{ff}(gi);
subplot(4,4,gi)
imagesc(log10(speccorr.(genotypes{gg}).freqs),speccorr.(genotypes{gg}).interpdepth,speccorr.(genotypes{gg}).pupinterp')
hold on
plot(log10(speccorr.(genotypes{gg}).freqs([1 end])),-depthinfo.boundaries'*[1 1],'k')

LogScale('x',10)
axis xy
ColorbarWithAxis([-0.4 0.3],'Pupil Corr.')
crameri('vik','pivot',0)
%clim([-0.35 0.35])
xlabel('f (Hz)');ylabel('Depth')
title(genotypes{gg})

subplot(4,4,gi+4)
imagesc(log10(speccorr.(genotypes{gg}).freqs),speccorr.(genotypes{gg}).interpdepth,speccorr.(genotypes{gg}).EMGinterp')
hold on
plot(log10(speccorr.(genotypes{gg}).freqs([1 end])),-depthinfo.boundaries'*[1 1],'k')
%plot(log10([band; band]),[ylim(gca)', ylim(gca)'],'k')
LogScale('x',10)
axis xy
ColorbarWithAxis([-0.4 0.3],'EMG Corr.')
crameri('vik','pivot',0)
xlabel('f (Hz)');ylabel('Depth')


end
NiceSave('DepthFreqCorr',analysisfolder,groupnames{ff})
end
end




  %%
  
  cosx = linspace(-pi,pi,100);
cospamp = [0.025 0.3]*2;

%speclim = [1.75 3.75]; %Wavelets?
speclim = [0.9 1.1]; %Med norm

for ff = 1:2
    for ww = 1:2
figure
for gi=1:length(groups{ff})
    gg = groups{ff}(gi);
for dd = 1:6

    
    subplot(6,4,(dd-1)*4+gi)
        for pp = 1:2
        imagesc( SPECdepth.(genotypes{gg}).(LAYERS{dd}).(HILO{pp}).(WHNWH{ww}).varbins+2*pi*(pp-1),...
            log10(SPECdepth.(genotypes{gg}).freqs),...
            SPECdepth.(genotypes{gg}).(LAYERS{dd}).(HILO{pp}).(WHNWH{ww}).mean)
        hold on; axis xy; box off
        plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp),'w')
        end   
        LogScale('y',10)
        %crameri batlow
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        ColorbarWithAxis(speclim,'Power (med^-^1)')
       ylim([0 2.5])
       crameri vik
        xlim([-pi 3*pi])
                plot(pi*[1 1],get(gca,'ylim'),'k--')
        bz_piTickLabel('x')
        if dd == 6
        xlabel('Pupil Phase');
        end
        if gi == 1
            ylabel({LAYERS{dd},'Freq'})
        end
        if dd == 1
        title({WHNWH{ww},genotypes{gg}})
        end

          

end
end
NiceSave(['DepthSPECandPupPhase',WHNWH{ww}],analysisfolder,groupnames{ff})

figure
for gi=1:length(groups{ff})
    gg = groups{ff}(gi);
for dd = 1:6
    subplot(6,4,(dd-1)*4+gi)
        imagesc( SPECdepth.(genotypes{gg}).(LAYERS{dd}).pup.(WHNWH{ww}).varbins,...
            log10(SPECdepth.(genotypes{gg}).freqs),...
            SPECdepth.(genotypes{gg}).(LAYERS{dd}).pup.(WHNWH{ww}).mean)
        hold on; axis xy; box off
        LogScale('y',10)
        ColorbarWithAxis(speclim,'Power (med^-^1)')
        %clim(speclim)
        ylim([0 2.5])
                plot(pi*[1 1],get(gca,'ylim'),'k--')
        bz_piTickLabel('x')
        crameri vik
        if dd == 6
        xlabel('Pupil Size');
        end
        if gi == 1
            ylabel({LAYERS{dd},'Freq'})
        end
        if dd == 1
        title({WHNWH{ww},genotypes{gg}})
        end
        
end
end
NiceSave(['DepthSPECandPup',WHNWH{ww}],analysisfolder,groupnames{ff})
end

end



%%
for ff = 1:2
    for ww = 1:2
figure
for gi=1:length(groups{ff})
    gg = groups{ff}(gi);
for dd = 1:6

    
    subplot(6,4,(dd-1)*4+gi)
        for pp = 1:2
        imagesc( OSCdepth.(genotypes{gg}).(LAYERS{dd}).(HILO{pp}).(WHNWH{ww}).varbins+2*pi*(pp-1),...
            log10(OSCdepth.(genotypes{gg}).freqs),...
            OSCdepth.(genotypes{gg}).(LAYERS{dd}).(HILO{pp}).(WHNWH{ww}).mean)
        hold on; axis xy; box off
        plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp),'w')
        end 
        
        LogScale('y',10)
        %crameri batlow
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        clim([0 0.25])
        %crameri('vik','pivot',0)
        xlim([-pi 3*pi])
        %ylim([0 2.5])
                plot(pi*[1 1],get(gca,'ylim'),'k--')
        bz_piTickLabel('x')
        if dd == 6
        xlabel('Pupil Phase');
        end
        if gi == 1
            ylabel({LAYERS{dd},'Freq'})
        end
        if dd == 1
        title({WHNWH{ww},genotypes{gg}})
        end

          

end
end
NiceSave(['DepthOSCandPupPhase',WHNWH{ww}],analysisfolder,groupnames{ff})


figure
for gi=1:length(groups{ff})
    gg = groups{ff}(gi);
for dd = 1:6
    subplot(6,4,(dd-1)*4+gi)
        imagesc( OSCdepth.(genotypes{gg}).(LAYERS{dd}).pup.(WHNWH{ww}).varbins,...
            log10(OSCdepth.(genotypes{gg}).freqs),...
            OSCdepth.(genotypes{gg}).(LAYERS{dd}).pup.(WHNWH{ww}).mean)
        hold on; axis xy; box off
        LogScale('y',10)
        %ylim([0 2.5])
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        clim([0 0.25])
        %crameri('vik','pivot',0)
                plot(pi*[1 1],get(gca,'ylim'),'k--')
        bz_piTickLabel('x')
        if dd == 6
        xlabel('Pupil Size');
        end
        if gi == 1
            ylabel({LAYERS{dd},'Freq'})
        end
        if dd == 1
        title({WHNWH{ww},genotypes{gg}})
        end
        
end
end
NiceSave(['DepthOSCandPup',WHNWH{ww}],analysisfolder,groupnames{ff})
end

end

%%

for ff = 1:2
    for oo = 1:2
figure
for gi=1:length(groups{ff})
    gg = groups{ff}(gi);
    
for dd = 1:6


subplot(6,4,(dd-1)*4+gi)
        imagesc( SPECdepth.(genotypes{gg}).(LAYERS{dd}).(ONOFF{oo}).all.varbins,...
            log10(SPECdepth.(genotypes{gg}).freqs),...
            SPECdepth.(genotypes{gg}).(LAYERS{dd}).(ONOFF{oo}).all.mean)
        hold on; axis xy; box off
        plot([0 0],[0 max(SPECdepth.(genotypes{gg}).freqs)],'w')
        ColorbarWithAxis(speclim,'Power (med^-^1)')
        crameri vik
        LogScale('y',10)
        ylim([0 2.5])
        if dd == 6
        xlabel(['t - ',(ONOFF{oo})])
        end
        if gi == 1
            ylabel({LAYERS{dd},'Freq'})
        end
        if dd == 1
        title({genotypes{gg},ONOFF{oo}})
        end
end
 
end

NiceSave(['DepthSPECandWhisk',ONOFF{oo}],analysisfolder,groupnames{ff})

    end
end

%%
for ff = 1:2

figure
for gi=1:length(groups{ff})
    gg = groups{ff}(gi);
    
for dd = 1:6

subplot(6,4,(dd-1)*4+gi)
        imagesc( SPECdepth.(genotypes{gg}).(LAYERS{dd}).EMG.varbins,...
            log10(SPECdepth.(genotypes{gg}).freqs),...
            SPECdepth.(genotypes{gg}).(LAYERS{dd}).EMG.mean)
        hold on; axis xy; box off
        %plot(SPECdepth.(genotypes{gg}).(LAYERS{dd}).EMG.varbins,SPECdepth.(genotypes{gg}).(LAYERS{dd}).EMG.vardist*10,'w')
        %plot(log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],[0 max(SPECdepth.freqs)],'k--')
        ColorbarWithAxis(speclim,'Power (med^-^1)')
        crameri vik
        ylim([0 2.5])
        if gi == 1
            ylabel({LAYERS{dd},'Freq'})
        end
        LogScale('y',10)
                if dd == 6
        xlabel('EMG')
        end
        if dd == 1
        title(genotypes{gg})
        end
end

end
    
NiceSave('DepthSPECandEMG',analysisfolder,groupnames{ff})

end


%%

for ff = 1:2
    for oo = 1:2
figure
for gi=1:length(groups{ff})
    gg = groups{ff}(gi);
    
for dd = 1:6


subplot(6,4,(dd-1)*4+gi)
        imagesc( OSCdepth.(genotypes{gg}).(LAYERS{dd}).(ONOFF{oo}).all.varbins,...
            log10(OSCdepth.(genotypes{gg}).freqs),...
            OSCdepth.(genotypes{gg}).(LAYERS{dd}).(ONOFF{oo}).all.mean)
        hold on; axis xy; box off
        plot([0 0],[0 max(OSCdepth.(genotypes{gg}).freqs)],'w')
        %ColorbarWithAxis(speclim,'Power (med^-^1)')
        caxis([0 0.25])
        %crameri vik
        LogScale('y',10)
        %ylim([0 2.5])
        if dd == 6
        xlabel(['t - ',(ONOFF{oo})])
        end
        if gi == 1
            ylabel({LAYERS{dd},'Freq'})
        end
        if dd == 1
        title({genotypes{gg},ONOFF{oo}})
        end
end
 
end

NiceSave(['DepthOSCandWhisk',ONOFF{oo}],analysisfolder,groupnames{ff})

    end
end

%%
for ff = 1:2

figure
for gi=1:length(groups{ff})
    gg = groups{ff}(gi);
    
for dd = 1:6

subplot(6,4,(dd-1)*4+gi)
        imagesc( OSCdepth.(genotypes{gg}).(LAYERS{dd}).EMG.varbins,...
            log10(OSCdepth.(genotypes{gg}).freqs),...
            OSCdepth.(genotypes{gg}).(LAYERS{dd}).EMG.mean)
        hold on; axis xy; box off
        %plot(SPECdepth.(genotypes{gg}).(LAYERS{dd}).EMG.varbins,SPECdepth.(genotypes{gg}).(LAYERS{dd}).EMG.vardist*10,'w')
        %plot(log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],[0 max(SPECdepth.freqs)],'k--')
        %ColorbarWithAxis(speclim,'Power (med^-^1)')
        %crameri vik
        caxis([0 0.25])
        %ylim([0 2.5])
        if gi == 1
            ylabel({LAYERS{dd},'Freq'})
        end
        LogScale('y',10)
                if dd == 6
        xlabel('EMG')
        end
        if dd == 1
        title(genotypes{gg})
        end
end

end
    
NiceSave('DepthOSCandEMG',analysisfolder,groupnames{ff})

end
%%
for ff = 1:2

figure
for gi=1:length(groups{ff})
    gg = groups{ff}(gi);
    
for dd = 1:6

subplot(6,4,(dd-1)*4+gi)
        a = imagesc( SPECdepth.(genotypes{gg}).(LAYERS{dd}).PSS.varbins,...
            log10(SPECdepth.(genotypes{gg}).freqs),...
            SPECdepth.(genotypes{gg}).(LAYERS{dd}).PSS.mean)
        alpha(a,single(~isnan(SPECdepth.(genotypes{gg}).(LAYERS{dd}).PSS.mean)))
        hold on; axis xy; box off
        plot(SPECdepth.(genotypes{gg}).(LAYERS{dd}).PSS.varbins,SPECdepth.(genotypes{gg}).(LAYERS{dd}).PSS.vardist*10,'k','linewidth',2)
        %plot(SPECdepth.(genotypes{gg}).(LAYERS{dd}).EMG.varbins,SPECdepth.(genotypes{gg}).(LAYERS{dd}).EMG.vardist*10,'w')
        %plot(log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],[0 max(SPECdepth.freqs)],'k--')
        ColorbarWithAxis([0.8 1.2],'Power (med^-^1)')
        crameri vik
        ylim([0 2.5])
        if gi == 1
            ylabel({LAYERS{dd},'Freq'})
        end
        LogScale('y',10)
                if dd == 6
        xlabel('PSS')
        end
        if dd == 1
        title(genotypes{gg})
        end
end

end
    
NiceSave('DepthSPECandPSS',analysisfolder,groupnames{ff})

end


%% Figure: PSS distributions

genecolors = {[0 0 0],'r','b','g',[0.5 0 0],[0 0 0.5],[0 0.5 0]};

for ff = 1:2

figure
for gi=1:length(groups{ff})
    gg = groups{ff}(gi);
    
for dd = 1:6

subplot(6,4,(dd-1)*4+gi)
        hold on; axis xy; box off
        errorshade(SPECdepth.(genotypes{gg}).(LAYERS{dd}).PSS.varbins,SPECdepth.(genotypes{gg}).(LAYERS{dd}).PSS.vardist,...
            SPECdepth_std.(genotypes{gg}).(LAYERS{dd}).PSS.vardist,SPECdepth_std.(genotypes{gg}).(LAYERS{dd}).PSS.vardist,genecolors{gi},'scalar')
        plot(SPECdepth.(genotypes{7}).(LAYERS{dd}).PSS.varbins,SPECdepth.(genotypes{7}).(LAYERS{dd}).PSS.vardist,'k','linewidth',2)
        plot(SPECdepth.(genotypes{gg}).(LAYERS{dd}).PSS.varbins,SPECdepth.(genotypes{gg}).(LAYERS{dd}).PSS.vardist,'color',genecolors{gi},'linewidth',2)

        %plot(SPECdepth.(genotypes{gg}).(LAYERS{dd}).EMG.varbins,SPECdepth.(genotypes{gg}).(LAYERS{dd}).EMG.vardist*10,'w')
        %plot(log10(EMGwhisk.detectorparms.Whthreshold).*[1 1],[0 max(SPECdepth.freqs)],'k--')
        if gi == 1
            ylabel({LAYERS{dd},'Occupancy'})
        end
        axis tight
        %ylim([0 0.1])
        
        if dd == 6
        xlabel('PSS')
        end
        if dd == 1
        title(genotypes{gg})
        end
end

end
    
NiceSave('DepthPSSDists',analysisfolder,groupnames{ff})

end