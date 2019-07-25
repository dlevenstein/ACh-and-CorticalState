analysisfolder = '/Users/dlevenstein/Project Repos/ACh-and-CorticalState/AnalysisScripts/AnalysisFigs/PupilFilterOptimization';
PupilFilterAll = GetMatResults(analysisfolder,'PupilFilterOptimization');
%genotype = {PupilWhiskAll.genotype};
%[genotypes,~,genotypeidx] = unique(genotype);
%%
PupilFilter = bz_CollapseStruct(PupilFilterAll);

%%
animalname = extractBefore(PupilFilter.name, '_');
genotype = extractBetween(PupilFilter.name,'_','.');
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
    PWcoupling.(genotypes{gg}) = bz_CollapseStruct(PupilFilter.PWcoupling(genotypeidx==gg),4,'median',true);
    pupilphaseEMG.(genotypes{gg})= bz_CollapseStruct(PupilFilter.pupilphaseEMG(genotypeidx==gg),4,'median',true);
	pupilphaseDUR.(genotypes{gg})= bz_CollapseStruct(PupilFilter.pupilphaseDUR(genotypeidx==gg),4,'median',true);
	filtercoupling.(genotypes{gg})= bz_CollapseStruct(PupilFilter.filtercoupling(genotypeidx==gg),4,'median',true);

end

PWcoupling.AllWT = bz_CollapseStruct(PupilFilter.PWcoupling(strcmp(WTKOtype,'WT')),4,'median',true);
pupilphaseEMG.AllWT = bz_CollapseStruct(PupilFilter.pupilphaseEMG(strcmp(WTKOtype,'WT')),4,'median',true);
pupilphaseDUR.AllWT = bz_CollapseStruct(PupilFilter.pupilphaseDUR(strcmp(WTKOtype,'WT')),4,'median',true);
filtercoupling.AllWT = bz_CollapseStruct(PupilFilter.filtercoupling(strcmp(WTKOtype,'WT')),4,'median',true);

%%
%Stored in PWOptimze later
lowerbounds = logspace(-2,-1,20);
upperbounds = logspace(-1.5,0,20);
orders = 1:4;

trybounds = [0.02 0.33];
%trybounds = [0.05 0.3];
gg = 7; %WT
figure
for oo = orders
    subplot(2,4,oo)
imagesc(log10(lowerbounds),log10(upperbounds),PWcoupling.(genotypes{gg}).EMG(:,:,oo)')
hold on
axis xy
plot(log10(trybounds(1)),log10(trybounds(2)),'r+')
LogScale('xy',10)
clim([0.25 0.8])
xlabel('Lower Bound (Hz)');ylabel('Upper Bounr (Hz)')
colorbar
title(num2str(oo))


    subplot(2,4,oo+4)
imagesc(log10(lowerbounds),log10(upperbounds),PWcoupling.(genotypes{gg}).dur(:,:,oo)')
hold on
axis xy
plot(log10(trybounds(1)),log10(trybounds(2)),'r+')
LogScale('xy',10)
clim([0.1 0.5])
xlabel('Lower Bound (Hz)');ylabel('Upper Bounr (Hz)')
colorbar
title(num2str(oo))
end

NiceSave('FilterCompare',analysisfolder,'WT')


%%
figure
cosx = linspace(-pi,pi,100);
cospamp = [0.08 0.8];
pupilcycle.pupthresh = -0.8;
figure
pp =1;

subplot(3,2,1)
a = imagesc(pupilphaseEMG.(genotypes{gg}).Xbins,pupilphaseEMG.(genotypes{gg}).Ybins,pupilphaseEMG.(genotypes{gg}).meanZ');
hold on
alpha(a,double(~isnan(pupilphaseEMG.(genotypes{gg}).meanZ')))

%imagesc(pupilphaseEMG.Xbins+2*pi,pupilphaseEMG.Ybins,pupilphaseEMG.meanZ')
crameri lapaz
plot([-pi 3*pi],pupilcycle.pupthresh.*[1 1],'w--')
plot(cosx,(cos(cosx)+1).*cospamp(pp)-2,'k')
axis xy
box off
%xlim([-pi 3*pi])
ColorbarWithAxis([-0.7 0.7],'Mean EMG')
LogScale('c',10)
LogScale('y',10)
xlabel('Pupil Phase');ylabel('Pupil Amplitude')
 bz_piTickLabel('x')
 ylim([-2 0.1])

 
subplot(3,2,2)
a = imagesc(pupilphaseDUR.(genotypes{gg}).Xbins,pupilphaseDUR.(genotypes{gg}).Ybins,pupilphaseDUR.(genotypes{gg}).meanZ');
hold on
alpha(a,double(~isnan(pupilphaseDUR.(genotypes{gg}).meanZ')))

%imagesc(pupilphaseEMG.Xbins+2*pi,pupilphaseEMG.Ybins,pupilphaseEMG.meanZ')
%crameri lapaz
plot([-pi 3*pi],pupilcycle.pupthresh.*[1 1],'w--')
plot(cosx,(cos(cosx)+1).*cospamp(pp)-2,'k')
axis xy
box off
%xlim([-pi 3*pi])
ColorbarWithAxis([-1 0.25],'Mean Dur')
LogScale('c',10)
LogScale('y',10)
xlabel('Pupil Phase');ylabel('Pupil Amplitude')
 bz_piTickLabel('x')
 ylim([-2 0.1])
 
subplot(3,2,3)
plot(filtercoupling.(genotypes{gg}).amp1bins,filtercoupling.(genotypes{gg}).EMGskew,'k')
xlabel('PupilAmplitude');ylabel('Coupling')
 
NiceSave('CouplingtoFilter',analysisfolder,'WT')
