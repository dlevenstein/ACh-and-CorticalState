%analysisfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/S1State/AnalysisScripts/AnalysisFigs/PupilWhiskAnalysis';
analysisfolder = '/Users/dlevenstein/Project Repos/ACh-and-CorticalState/AnalysisScripts/AnalysisFigs/PupilWhiskAnalysis';

PupilWhiskAll = GetMatResults(analysisfolder,'PupilWhiskAnalysis');
PupilWhiskAll = bz_CollapseStruct(PupilWhiskAll);

%%
animalname = extractBefore(PupilWhiskAll.name, '_');
genotype = extractBetween(PupilWhiskAll.name,'_','_');
genotype = genotype(:,:,1);
[genotypes,~,genotypeidx] = unique(genotype);
%genotypes = char(genotypes);

genotypecolors = {'r','k'};

%%
pupilwhisk = CollapseStruct(PupilWhiskAll);


%% Pupil Stats

pupilPSD = CollapseStruct(pupilwhisk.pupPSD,2);
pupilPSD.freqs = pupilwhisk.pupPSD.freqs;

pupdthist = CollapseStruct(pupilwhisk.pupdthist,3);
pupdthist.bins = pupilwhisk.pupdthist.bins;

pupilhist = CollapseStruct(pupilwhisk.puphist,1);

pupilacg = CollapseStruct(pupilwhisk.pupACG);
pupilacg.tlag = pupilwhisk.pupACG.tlag;

for gg = 1:length(genotypes)
    pupilhist.(genotypes{gg}) = mean(pupilhist.counts(genotypeidx==gg,:),1);
    pupilhist.std.(genotypes{gg}) = std(pupilhist.counts(genotypeidx==gg,:),[],1);

    pupdthist.(genotypes{gg}) = mean(pupdthist.counts(:,:,genotypeidx==gg),3);
    
    pupilPSD.(genotypes{gg}) = mean(pupilPSD.psd(:,genotypeidx==gg),2);
    pupilPSD.std.(genotypes{gg}) = std(pupilPSD.psd(:,genotypeidx==gg),[],2);
    
    pupilacg.(genotypes{gg}) = mean(pupilacg.ACG(:,genotypeidx==gg),2);
    pupilacg.std.(genotypes{gg}) =std(pupilacg.ACG(:,genotypeidx==gg),[],2);
end

%% Whisk Stats

emghist = CollapseStruct(pupilwhisk.EMGhist,1);
emghist.logbins = pupilwhisk.EMGhist.logbins;
%geneotype = CollapseStruct(pupilwhisk.genotype);

whdurhist = CollapseStruct(pupilwhisk.Whdurhist,1);
whdurhist.bins = pupilwhisk.Whdurhist.bins;


for gg = 1:length(genotypes)
    emghist.(genotypes{gg}) = mean(emghist.logcounts(genotypeidx==gg,:),1);
    emghist.std.(genotypes{gg}) = std(emghist.logcounts(genotypeidx==gg,:),[],1);

    whdurhist.(genotypes{gg}).Wh = mean(whdurhist.Whdurs(genotypeidx==gg,:),1);
    whdurhist.(genotypes{gg}).IWh = mean(whdurhist.InterWhdurs(genotypeidx==gg,:),1);
    
    whdurhist.std.(genotypes{gg}).Wh = std(whdurhist.Whdurs(genotypeidx==gg,:),[],1);
    whdurhist.std.(genotypes{gg}).IWh = std(whdurhist.InterWhdurs(genotypeidx==gg,:),[],1);
end

%% Pupil-Whisk Joint Stats

pupilEMGdist = CollapseStruct(pupilwhisk.pupilEMGdist,3);
pupilEMGdist.bins = pupilwhisk.pupilEMGdist.bins;

pwCCG = CollapseStruct(pupilwhisk.pwCCG,2,'justcat',true);
pwCCG.t = ([1:length(pupilwhisk.pwCCG(1).EMG.WhOn)]-0.5.*length(pupilwhisk.pwCCG(1).EMG.WhOn))./50;

pupildynamicsEMG = CollapseStruct(pupilwhisk.pupildynamicsEMG,3,'justcat',true);
pupildynamicsEMG.bins = pupilwhisk.pupildynamicsEMG.bins;

pupilphaseEMG = CollapseStruct(pupilwhisk.pupilphaseEMG,3,'justcat',true);
pupilphaseEMG.bins = pupilwhisk.pupilphaseEMG.bins;

WPcoupling = CollapseStruct(pupilwhisk.WPcoupling);
WPcoupling.freqs = pupilwhisk.WPcoupling.freqs;

for gg = 1:length(genotypes)
    pupilEMGdist.(genotypes{gg}) = mean(pupilEMGdist.counts(:,:,genotypeidx==gg),3);

    pwCCG.(genotypes{gg}).EMG = mean(pwCCG.EMG.WhOn(:,genotypeidx==gg),2);
    pwCCG.(genotypes{gg}).pupil = mean(pwCCG.pupil.WhOn(:,genotypeidx==gg),2);
    
    pwCCG.std.(genotypes{gg}).EMG = std(pwCCG.EMG.WhOn(:,genotypeidx==gg),[],2);
    pwCCG.std.(genotypes{gg}).pupil = std(pwCCG.pupil.WhOn(:,genotypeidx==gg),[],2);
    
    pupildynamicsEMG.(genotypes{gg}).meanEMG = nanmean(pupildynamicsEMG.mean(:,:,genotypeidx==gg),3);
    pupildynamicsEMG.(genotypes{gg}).startrate = mean(pupildynamicsEMG.Whstartrate(:,:,genotypeidx==gg),3);
    pupildynamicsEMG.(genotypes{gg}).Whdur = nanmean(pupildynamicsEMG.meanWhdur(:,:,genotypeidx==gg),3);
    
    pupilphaseEMG.(genotypes{gg}).meanEMG = nanmean(pupilphaseEMG.mean(:,:,genotypeidx==gg),3);
    pupilphaseEMG.(genotypes{gg}).startrate = mean(pupilphaseEMG.Whstartrate(:,:,genotypeidx==gg),3);
    pupilphaseEMG.(genotypes{gg}).Whdur = nanmean(pupilphaseEMG.meanWhdur(:,:,genotypeidx==gg),3);

    WPcoupling.(genotypes{gg}).coupling = mean(WPcoupling.coupling(:,genotypeidx==gg),2);
    WPcoupling.std.(genotypes{gg}).coupling = std(WPcoupling.coupling(:,genotypeidx==gg),[],2);

end


%% PupiWhisk: Coupling
figure
subplot(2,2,1)
    hold on
        for gg = 1:length(genotypes)
            plot(log10(WPcoupling.freqs),WPcoupling.(genotypes{gg}).coupling,...
                'color',genotypecolors{gg},'linewidth',1)
            errorshade(log10(WPcoupling.freqs),WPcoupling.(genotypes{gg}).coupling,...
                WPcoupling.std.(genotypes{gg}).coupling,WPcoupling.std.(genotypes{gg}).coupling ,genotypecolors{gg},'scalar')
        end
      
       LogScale('x',10)
       % ylim([0 0.35])
        xlabel('Pupil freq (Hz)');ylabel('Pupil-EMG Coupling')

%% PupilWhisk: PupilSpace
figure
    for gg = 1:length(genotypes)   
    subplot(3,3,gg)
        imagesc((pupilphaseEMG.bins),pupilphaseEMG.bins,pupilphaseEMG.(genotypes{gg}).meanEMG)
        hold on
        imagesc((pupilphaseEMG.bins)+2*pi,pupilphaseEMG.bins,pupilphaseEMG.(genotypes{gg}).meanEMG)

        colorbar
        title(genotypes{gg})
        axis xy
        caxis([-1 0.5])
        xlim([-pi 3*pi]);ylim([-2 1])

        LogScale('c',10)
        xlabel('Pupil Phase');ylabel('Amp.')
    end

    for gg = 1:length(genotypes)   
    subplot(3,3,gg+3)
        imagesc((pupilphaseEMG.bins),pupilphaseEMG.bins,pupilphaseEMG.(genotypes{gg}).startrate)
        hold on
        imagesc((pupilphaseEMG.bins)+2*pi,pupilphaseEMG.bins,pupilphaseEMG.(genotypes{gg}).startrate)

        colorbar
        title(genotypes{gg})
        axis xy
        caxis([0 1])
        xlim([-pi 3*pi]);ylim([-2 1])

        %LogScale('c',10)
        xlabel('Pupil Phase');ylabel('Amp.')
    end
    
    for gg = 1:length(genotypes)   
    subplot(3,3,gg+6)
        imagesc((pupilphaseEMG.bins),pupilphaseEMG.bins,log10(pupilphaseEMG.(genotypes{gg}).Whdur))
        hold on
        imagesc((pupilphaseEMG.bins)+2*pi,pupilphaseEMG.bins,log10(pupilphaseEMG.(genotypes{gg}).Whdur))

        colorbar
        title(genotypes{gg})
        axis xy
        caxis([-1 0.5])
        xlim([-pi 3*pi]);ylim([-2 1])
        %LogScale('x',10)
        LogScale('c',10)
        xlabel('Pupil Phase');ylabel('Amp.')
    end
    
NiceSave('PWPhaseSpace',analysisfolder,[])
    
%%    
    
figure
    for gg = 1:length(genotypes)   
    subplot(3,3,gg)
        imagesc((pupildynamicsEMG.bins),pupildynamicsEMG.bins,pupildynamicsEMG.(genotypes{gg}).meanEMG')
        colorbar
        title(genotypes{gg})
        axis xy
        caxis([-1 0.5])
        LogScale('x',10)
        LogScale('c',10)
        xlabel('Pupil Area');ylabel('dpdt')
    end

    for gg = 1:length(genotypes)   
    subplot(3,3,gg+3)
        imagesc((pupildynamicsEMG.bins),pupildynamicsEMG.bins,pupildynamicsEMG.(genotypes{gg}).startrate')
        colorbar
        title(genotypes{gg})
        axis xy
        caxis([0 1])
        LogScale('x',10)
        %LogScale('c',10)
        xlabel('Pupil Area');ylabel('dpdt')
    end
    
    for gg = 1:length(genotypes)   
    subplot(3,3,gg+6)
        imagesc((pupildynamicsEMG.bins),pupildynamicsEMG.bins,log10(pupildynamicsEMG.(genotypes{gg}).Whdur)')
        colorbar
        title(genotypes{gg})
        axis xy
        caxis([-1 0.5])
        LogScale('x',10)
        LogScale('c',10)
        xlabel('Pupil Area');ylabel('dpdt')
    end
    
NiceSave('PWPupilSpace',analysisfolder,[])


%% Figure Pupil-Whisk
figure
    for gg = 1:length(genotypes)   
    subplot(3,3,gg+6)
        imagesc((pupilEMGdist.bins{1}),pupilEMGdist.bins{2},log10(pupilEMGdist.(genotypes{gg}))')
        colorbar
        title(genotypes{gg})
        axis xy
        caxis([-5 -2.75])
        LogScale('xy',10)
        %LogScale('c',10)
        xlabel('Pupil Area');ylabel('EMG')
    end

    
    subplot(4,3,1)
    hold on
        for gg = 1:length(genotypes)
            plot(pwCCG.t,pwCCG.(genotypes{gg}).EMG,...
                'color',genotypecolors{gg},'linewidth',1)
            errorshade(pwCCG.t,pwCCG.(genotypes{gg}).EMG,...
                pwCCG.std.(genotypes{gg}).EMG,pwCCG.std.(genotypes{gg}).EMG ,genotypecolors{gg},'scalar')
        end
       xlim([-100 100])
       % LogScale('x',10)
       % ylim([0 0.35])
        xlabel('EMG');ylabel('Occupancy')
        
    subplot(4,3,1)
    hold on
        for gg = 1:length(genotypes)
            plot(pwCCG.t,pwCCG.(genotypes{gg}).EMG,...
                'color',genotypecolors{gg},'linewidth',1)
            errorshade(pwCCG.t,pwCCG.(genotypes{gg}).EMG,...
                pwCCG.std.(genotypes{gg}).EMG,pwCCG.std.(genotypes{gg}).EMG ,genotypecolors{gg},'scalar')
        end
       xlim([-4 4])
       % LogScale('x',10)
       % ylim([0 0.35])
        xlabel('t (s)');ylabel('EMG')
        
    subplot(4,3,2)
    hold on
        for gg = 1:length(genotypes)
            plot(pwCCG.t,pwCCG.(genotypes{gg}).pupil,...
                'color',genotypecolors{gg},'linewidth',1)
            errorshade(pwCCG.t,pwCCG.(genotypes{gg}).pupil,...
                pwCCG.std.(genotypes{gg}).pupil,pwCCG.std.(genotypes{gg}).pupil ,genotypecolors{gg},'scalar')
        end
        axis tight
       xlim([-10 10])
       % LogScale('x',10)
       % ylim([0 0.35])
        xlabel('t (s)');ylabel('Pupil')
        

%% Figure : Whisk stats
figure
    subplot(4,3,1)
    hold on
        for gg = 1:length(genotypes)
            plot(emghist.logbins,emghist.(genotypes{gg}),...
                'color',genotypecolors{gg},'linewidth',1)
            errorshade(emghist.logbins,emghist.(genotypes{gg}),...
                emghist.std.(genotypes{gg}),emghist.std.(genotypes{gg}),genotypecolors{gg},'scalar')
        end
        xlim(emghist.logbins(1,[1 end]))
        LogScale('x',10)
        ylim([0 0.35])
        xlabel('EMG');ylabel('Occupancy')
        
    subplot(4,3,2)
    hold on
        for gg = 1:length(genotypes)
            plot(whdurhist.bins,whdurhist.(genotypes{gg}).Wh,...
                'color',genotypecolors{gg},'linewidth',1)
            errorshade(whdurhist.bins,whdurhist.(genotypes{gg}).Wh,...
                whdurhist.std.(genotypes{gg}).Wh,whdurhist.std.(genotypes{gg}).Wh,...
                genotypecolors{gg},'scalar')
        end
        xlim(whdurhist.bins(1,[1 end]))
        LogScale('x',10)
        axis tight
        xlabel('Whisk Duration');ylabel('Rate (s^-^1)')
        
    subplot(4,3,3)
    hold on
        for gg = 1:length(genotypes)
            plot(whdurhist.bins,whdurhist.(genotypes{gg}).IWh,...
                'color',genotypecolors{gg},'linewidth',1)
            errorshade(whdurhist.bins,whdurhist.(genotypes{gg}).IWh,...
                whdurhist.std.(genotypes{gg}).IWh,whdurhist.std.(genotypes{gg}).IWh,...
                genotypecolors{gg},'scalar')
        end
        xlim(whdurhist.bins(1,[1 end]))
        LogScale('x',10)
        axis tight
        ylim([0 0.07])
        xlabel('Inter-Whisk Duration');ylabel('Rate (s^-^1)')
        
NiceSave('WhiskStats',analysisfolder,[])

%% Figure: pupil stats
figure
    subplot(4,3,1)
    hold on
        for gg = 1:length(genotypes)
            plot(pupilhist.bins(1,:),pupilhist.(genotypes{gg}),...
                'color',genotypecolors{gg},'linewidth',1)
            errorshade(pupilhist.bins(1,:),pupilhist.(genotypes{gg}),...
                pupilhist.std.(genotypes{gg}),pupilhist.std.(genotypes{gg}),genotypecolors{gg},'scalar')
        end
        xlim(pupilhist.bins(1,[1 end]))
        ylim([0 0.35])
        xlabel('Pupil Area');ylabel('Occupancy')
    subplot(4,3,2)
    hold on
        for gg = 1:length(genotypes)
            plot(log10(pupilPSD.freqs),pupilPSD.(genotypes{gg}),...
                'color',genotypecolors{gg},'linewidth',1)
            errorshade(log10(pupilPSD.freqs),pupilPSD.(genotypes{gg}),...
                pupilPSD.std.(genotypes{gg}),pupilPSD.std.(genotypes{gg}),genotypecolors{gg},'scalar')
        end
        LogScale('x',10)
        axis tight
        xlabel('Frequency (Hz)');ylabel('Power (dB)')
    
    subplot(4,3,3)
        plot(pupilacg.tlag([1 end]),[0 0],'k')
    hold on
        for gg = 1:length(genotypes)
            plot(pupilacg.tlag,pupilacg.(genotypes{gg}),...
                'color',genotypecolors{gg},'linewidth',1)
            errorshade(pupilacg.tlag,pupilacg.(genotypes{gg}),...
                pupilacg.std.(genotypes{gg}),pupilacg.std.(genotypes{gg}),genotypecolors{gg},'scalar')
        end
        axis tight
        xlabel('t lag (s)');ylabel('XCorr')
        xlim([-150 150])
        
    for gg = 1:length(genotypes)   
    subplot(3,3,6+gg)
        imagesc((pupdthist.bins{1}),pupdthist.bins{2},log10(pupdthist.(genotypes{gg}))')
        colorbar
        title(genotypes{gg})
        axis xy
        caxis([-5 -2])
        LogScale('x',10)
        LogScale('c',10)
        xlabel('Pupil Area');ylabel('dpdt')
    end
    
NiceSave('PupilStats',analysisfolder,[])

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