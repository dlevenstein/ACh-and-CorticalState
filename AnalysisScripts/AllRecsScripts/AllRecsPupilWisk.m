analysisfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/S1State/AnalysisScripts/AnalysisFigs/PupilWhiskAnalysis';
PupilWhiskAll = GetMatResults(analysisfolder,'_PupilWhiskAnalysis');
genotype = {PupilWhiskAll.genotype};
[genotypes,~,genotypeidx] = unique(genotype);

%%
pupilwhisk = CollapseStruct(PupilWhiskAll);

pupilhist = CollapseStruct(pupilwhisk.puphist,1);
acg_t = pupilwhisk.pupACG.tlag;
pupilacg = CollapseStruct(pupilwhisk.pupACG);

emghist = CollapseStruct(pupilwhisk.EMGhist,1);

%geneotype = CollapseStruct(pupilwhisk.genotype);

whdurhist = CollapseStruct(pupilwhisk.Whdurhist,1);
%%
colors = {'r','g','k'};
figure
% subplot(2,2,1)
%     imagesc(pupilhist.counts)
subplot(3,2,2)
hold on
for gg = 1:length(genotypes)
    plot(pupilhist.bins(1,:),pupilhist.counts(genotypeidx==gg,:)','linewidth',2,'color',colors{gg})
    xlabel('Pupil Diameter')
end
subplot(3,2,4)
hold on
for gg = 1:length(genotypes)
    plot(acg_t,pupilacg.ACG(:,genotypeidx==gg),'linewidth',2,'color',colors{gg})
    xlim([-60 60])
end
subplot(3,2,6)
hold on
for gg = 1:length(genotypes)
    plot(emghist.logbins(1,:),emghist.logcounts(genotypeidx==gg,:)','linewidth',2,'color',colors{gg})
    xlabel('EMG')
    LogScale('x',10)
end
%legend(genotypes)

subplot(3,2,1)
hold on
for gg = 1:length(genotypes)
    plot(whdurhist.bins(1,:),whdurhist.Whdurs(genotypeidx==gg,:),'linewidth',2,'color',colors{gg})
end
title('Whisk Durations')
LogScale('x',10)
xlabel('duration (s)')

subplot(3,2,3)
hold on
for gg = 1:length(genotypes)
    plot(whdurhist.bins(1,:),whdurhist.InterWhdurs(genotypeidx==gg,:),'linewidth',2,'color',colors{gg})
end
title('Inter-Whisk Durations')
xlabel('duration (s)')
LogScale('x',10)

NiceSave('PupEMGHet',analysisfolder,[])
%%
pupilemgxcorr = CollapseStruct(pupilwhisk.pwCCG,2,'justcat',true)

figure
plot(pupilemgxcorr.pupil.WhOn)