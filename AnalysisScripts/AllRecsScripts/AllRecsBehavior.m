analysisfolder = '/Users/dlevenstein/Project Repos/ACh-and-CorticalState/AnalysisScripts/AnalysisFigs/BehaviorAnalysis2';
BehaviorAll = GetMatResults(analysisfolder,'BehaviorAnalysis2');
%genotype = {PupilWhiskAll.genotype};
%[genotypes,~,genotypeidx] = unique(genotype);
%%
Behavior = bz_CollapseStruct(BehaviorAll);

%%
animalname = extractBefore(Behavior.name, '_');
genotype = extractBetween(Behavior.name,'_','.');
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
    EMGdist.(genotypes{gg}) = bz_CollapseStruct(Behavior.EMGdist(genotypeidx==gg),3,'mean',true);
    EMGdur.(genotypes{gg}) = bz_CollapseStruct(Behavior.EMGdur(genotypeidx==gg),3,'mean',true);
    EMGdur_std.(genotypes{gg}) = bz_CollapseStruct(Behavior.EMGdur(genotypeidx==gg),3,'std',true);

    pupilphaseEMG.(genotypes{gg}) = bz_CollapseStruct(Behavior.pupilphaseEMG(genotypeidx==gg),3,'mean',true);
    pupildpEMG.(genotypes{gg}) = bz_CollapseStruct(Behavior.pupildpEMG(genotypeidx==gg),3,'mean',true);
end

EMGdist.AllWT = bz_CollapseStruct(Behavior.EMGdist(strcmp(WTKOtype,'WT')),3,'mean',true);
EMGdur.AllWT = bz_CollapseStruct(Behavior.EMGdur(strcmp(WTKOtype,'WT')),3,'mean',true);
EMGdur_std.AllWT = bz_CollapseStruct(Behavior.EMGdur(strcmp(WTKOtype,'WT')),3,'std',true);
pupilphaseEMG.AllWT = bz_CollapseStruct(Behavior.pupilphaseEMG(strcmp(WTKOtype,'WT')),3,'mean',true);
pupildpEMG.AllWT = bz_CollapseStruct(Behavior.pupildpEMG(strcmp(WTKOtype,'WT')),3,'mean',true);

%%
HILO = {'lopup','hipup'};
pupilcycle.pupthresh = -0.8;
%%
cosx = linspace(-pi,pi,100);
cospamp = [0.08 0.6];

for ff = 1:2
figure
for gi=1:length(groups{ff})
    gg = groups{ff}(gi);


subplot(4,4,gi)
a = imagesc(pupilphaseEMG.(genotypes{gg}).Xbins,pupilphaseEMG.(genotypes{gg}).Ybins,pupilphaseEMG.(genotypes{gg}).meanZ');
hold on
alpha(a,double(~isnan(pupilphaseEMG.(genotypes{gg}).meanZ')))

%imagesc(pupilphaseEMG.Xbins+2*pi,pupilphaseEMG.Ybins,pupilphaseEMG.meanZ')
crameri lapaz
plot([-pi 3*pi],pupilcycle.pupthresh.*[1 1],'w--')
plot(cosx,(cos(cosx)+1).*cospamp(2)-2,'w')
axis xy
box off
%xlim([-pi 3*pi])
ylim([-2 0.1])
ColorbarWithAxis([-0.7 0.7],'Mean EMG')
LogScale('c',10)
LogScale('y',10)
        bz_piTickLabel('x',0.5)
title(genotypes{gg})
xlabel('Pupil Phase');ylabel('Pupil Amplitude')


subplot(4,4,4+gi)
        for pp = 1:2
        imagesc( EMGdist.(genotypes{gg}).(HILO{pp}).Xbins+2*pi*(pp-1),...
            EMGdist.(genotypes{gg}).(HILO{pp}).Ybins,...
            EMGdist.(genotypes{gg}).(HILO{pp}).pYX')
        hold on; axis xy; box off
        plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp)-1.75,'k')
        end   
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        %colorbar
        LogScale('y',10)
        xlim([-pi 3*pi])
                plot(pi*[1 1],get(gca,'ylim'),'k--')
        bz_piTickLabel('x')
        xlabel('Pupil Phase');ylabel('EMG')
        crameri bilbao
        



subplot(4,4,8+gi)
        for pp = 1:2
        imagesc( EMGdur.(genotypes{gg}).(HILO{pp}).Xbins+2*pi*(pp-1),...
            EMGdur.(genotypes{gg}).(HILO{pp}).Ybins,...
            EMGdur.(genotypes{gg}).(HILO{pp}).pYX')
        hold on; axis xy; box off
        %plot(EMGdur.(genotypes{gg}).(HILO{pp}).Xbins+2*pi*(pp-1),EMGdur.(genotypes{gg}).(HILO{pp}).pWhisk,'k')
        plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp)-1,'k')
        end   
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        %colorbar
        xlim([-pi 3*pi])
        LogScale('y',10)
        plot(pi*[1 1],get(gca,'ylim'),'k--')
        bz_piTickLabel('x')
        xlabel('Pupil Phase');ylabel('Duration (s)')
        crameri bilbao
        

subplot(4,4,12+gi)
        for pp = 1:2
        hold on; box off
        plot(EMGdur.(genotypes{gg}).(HILO{pp}).Xbins+2*pi*(pp-1),EMGdur.(genotypes{gg}).(HILO{pp}).pWhisk,'k')
        errorshade(EMGdur.(genotypes{gg}).(HILO{pp}).Xbins+2*pi*(pp-1),EMGdur.(genotypes{gg}).(HILO{pp}).pWhisk,...
            EMGdur_std.(genotypes{gg}).(HILO{pp}).pWhisk,EMGdur_std.(genotypes{gg}).(HILO{pp}).pWhisk,'k','scalar')
        %plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp),'k')
        end   
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        %colorbar
        xlim([-pi 3*pi])
        ylim([0 0.6])
                plot(pi*[1 1],get(gca,'ylim'),'k--')
        bz_piTickLabel('x')
        xlabel('Pupil Phase');ylabel('P[Whisk] (Hz)')
        crameri bilbao  
   
%subplot(12,4,4)

end
       NiceSave('EMGbyPupilPhase',analysisfolder,groupnames{ff},'figtype','pdf')

end


%%

for ff = 1:2
figure
for gi=1:length(groups{ff})
    gg = groups{ff}(gi);

subplot(4,4,gi)
    a = imagesc(pupildpEMG.(genotypes{gg}).Xbins,pupildpEMG.(genotypes{gg}).Ybins,pupildpEMG.(genotypes{gg}).meanZ');
    hold on
    alpha(a,double(~isnan(pupildpEMG.(genotypes{gg}).meanZ')))
    plot(pupildpEMG.(genotypes{gg}).Xbins([1 end]),[0 0],'k--')
    crameri lapaz
    ColorbarWithAxis([-0.7 0.7],'Mean EMG')
    LogScale('c',10)
    axis xy
    box off
    LogScale('c',10)
    LogScale('x',10)
    xlabel('Pupil Size');ylabel('dp/dt')
    

subplot(4,4,4+gi)

        imagesc( EMGdist.(genotypes{gg}).pup.Xbins,...
            EMGdist.(genotypes{gg}).pup.Ybins,...
            EMGdist.(genotypes{gg}).pup.pYX')
        hold on; axis xy; box off
        %colorbar
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        xlabel('Pupil Size');ylabel('EMG')
   crameri bilbao
   
subplot(4,4,8+gi)

        imagesc( EMGdur.(genotypes{gg}).pup.Xbins,...
            EMGdur.(genotypes{gg}).pup.Ybins,...
            EMGdur.(genotypes{gg}).pup.pYX')
        hold on; axis xy; box off
        %plot(EMGdur.(genotypes{gg}).pup.Xbins,EMGdur.(genotypes{gg}).pup.pWhisk,'k')
        %colorbar
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        xlabel('Pupil Size');ylabel('Dur')
        LogScale('y',10)
   crameri bilbao

subplot(4,4,12+gi)

        plot(EMGdur.(genotypes{gg}).pup.Xbins,EMGdur.(genotypes{gg}).pup.pWhisk,'k')
        errorshade(EMGdur.(genotypes{gg}).pup.Xbins,EMGdur.(genotypes{gg}).pup.pWhisk,...
            EMGdur_std.(genotypes{gg}).pup.pWhisk,EMGdur_std.(genotypes{gg}).pup.pWhisk,'k','scalar')
        %plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp),'k') 
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        %colorbar
        %xlim([-pi 3*pi])
        xlim(EMGdur.(genotypes{gg}).pup.Xbins([1 end]))
        ylim([0 0.6])
        xlabel('Pupil Size');ylabel('P[Whisk] (Hz)')
        %crameri bilbao  
   
   
end
NiceSave('EMGbyPupilSize',analysisfolder,groupnames{ff},'figtype','pdf')

end

