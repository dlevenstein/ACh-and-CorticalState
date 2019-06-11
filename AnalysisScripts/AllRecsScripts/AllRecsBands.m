analysisfolder = '/Users/dlevenstein/Project Repos/ACh-and-CorticalState/AnalysisScripts/AnalysisFigs/FBandsAnalysis';
FBandsAll = GetMatResults(analysisfolder,'FBandsAnalysis');
%genotype = {PupilWhiskAll.genotype};
%[genotypes,~,genotypeidx] = unique(genotype);
%%
FBands = bz_CollapseStruct(FBandsAll);

%%
animalname = extractBefore(FBands.name, '_');
genotype = extractBetween(FBands.name,'_','.');
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
    BANDdepth.(genotypes{gg}) = bz_CollapseStruct(FBands.BANDdepth(genotypeidx==gg),3,'mean',true);
end

BANDdepth.AllWT = bz_CollapseStruct(FBands.BANDdepth(strcmp(WTKOtype,'WT')),3,'mean',true);


%%

WHNWH = {'Wh','NWh'};
HILO = {'lopup','hipup'};
LAYERS = {'L1','L23','L4','L5a','L5b6','L6'};
depthinfo.boundaries = [0 0.1 0.35 0.5 1];
BANDS = {'HiGamma','LoGamma','Theta','Delta','Slow'};
ONOFF = {'WhOn','WhOFF'};
%bandranges = [[100 312];[30 55];[6 12];[2 5];[0.5 2]];




%%
cosx = linspace(-pi,pi,100);
cospamp = [0.025 0.3];

bandranges = [[1.7 2.05];[1.85 2.1];[2.15 2.6];[2.15 2.75];[2.2 2.75]];

%gg = 7;

for ff = 1:2
    for bb = 1:length(BANDS)
figure

for gi=1:length(groups{ff})
    gg = groups{ff}(gi);
    
    

for ww = 1:2
    subplot(4,4,(ww-1)*4+gi)
        for pp = 1:2
        imagesc( BANDdepth.(genotypes{gg}).(BANDS{bb}).(HILO{pp}).(WHNWH{ww}).varbins+2*pi*(pp-1),...
            BANDdepth.(genotypes{gg}).(BANDS{bb}).depth,...
            BANDdepth.(genotypes{gg}).(BANDS{bb}).(HILO{pp}).(WHNWH{ww}).mean_interp)
        hold on; axis xy; box off
        plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp)-1,'k')
        end   
        ColorbarWithAxis(bandranges(bb,:),'Power (dB)')
        plot([-pi 3*pi],-depthinfo.boundaries'*[1 1],'w')
        colorbar
        xlim([-pi 3*pi])
        xlabel('Pupil Phase');
        if gi == 1
            ylabel({WHNWH{ww},'Depth'})
        end
        if ww == 1
            title({BANDS{bb},genotypes{gg}})
        end
        
        
    subplot(4,4,(ww-1)*4+8+gi)
        imagesc( BANDdepth.(genotypes{gg}).(BANDS{bb}).pup.(WHNWH{ww}).varbins,...
            BANDdepth.(genotypes{gg}).(BANDS{bb}).depth,...
            BANDdepth.(genotypes{gg}).(BANDS{bb}).pup.(WHNWH{ww}).mean_interp)
        hold on; axis xy; box off
        plot(BANDdepth.(genotypes{gg}).(BANDS{bb}).pup.(WHNWH{ww}).varbins([1 end]),-depthinfo.boundaries'*[1 1],'w')
        ColorbarWithAxis(bandranges(bb,:),'Power (dB)')
        colorbar
        xlabel('Pupil Size');ylabel('Depth')
        
        if gi == 1
            ylabel({WHNWH{ww},'Depth'})
        end

        
end
end
NiceSave(['FBandsPupil',(BANDS{bb})],analysisfolder,groupnames{ff},'figtype','pdf')

    end
end
%%

for ff = 1:2
    for bb = 1:length(BANDS)
figure

for gi=1:length(groups{ff})
    gg = groups{ff}(gi);

%for oo = 1:2
oo=1;
subplot(4,4,gi)
        imagesc( BANDdepth.(genotypes{gg}).(BANDS{bb}).(ONOFF{oo}).all.varbins,...
            BANDdepth.(genotypes{gg}).(BANDS{bb}).depth,...
            BANDdepth.(genotypes{gg}).(BANDS{bb}).(ONOFF{oo}).all.mean_interp)
        hold on; axis xy; box off
        plot([0 0],[-1 0],'w')
       ColorbarWithAxis(bandranges(bb,:),'Power (dB)')
        plot(BANDdepth.(genotypes{gg}).(BANDS{bb}).(ONOFF{oo}).all.varbins([1 end]),-depthinfo.boundaries'*[1 1],'w')

        xlabel(['t - ',(ONOFF{oo})]);
        %if oo == 1   
        if gi == 1
            ylabel('Depth')
        end

            title({BANDS{bb},genotypes{gg}})

        %end
        xlim([-0.5 0.5])

%end
 
subplot(4,4,gi+4)
        imagesc( BANDdepth.(genotypes{gg}).(BANDS{bb}).EMG.varbins,...
            BANDdepth.(genotypes{gg}).(BANDS{bb}).depth,...
            BANDdepth.(genotypes{gg}).(BANDS{bb}).EMG.mean_interp)
        hold on; axis xy; box off
        plot(BANDdepth.(genotypes{gg}).(BANDS{bb}).EMG.varbins([1 end]),-depthinfo.boundaries'*[1 1],'w')

        plot(BANDdepth.(genotypes{gg}).(BANDS{bb}).EMG.varbins,BANDdepth.(genotypes{gg}).(BANDS{bb}).EMG.vardist*10-1,'k')
       ColorbarWithAxis(bandranges(bb,:),'Power (dB)')
        xlabel('EMG');
        if gi == 1
            ylabel('Depth')
        end
        
        
end
NiceSave(['FBandsWhisk',(BANDS{bb})],analysisfolder,groupnames{ff},'figtype','pdf')

    end
end

