analysisfolder = '/Users/dlevenstein/Project Repos/ACh-and-CorticalState/AnalysisScripts/AnalysisFigs/UPDOWNandPupilAnalysis';
UPDOWNPupilAll = GetMatResults(analysisfolder,'UPDOWNandPupilAnalysis');
%genotype = {PupilWhiskAll.genotype};
%[genotypes,~,genotypeidx] = unique(genotype);

%%
UPDOWNpupil = CollapseStruct(UPDOWNPupilAll);

%%
pupcycleUPDOWN = CollapseStruct(UPDOWNpupil.pupcycleUPDOWN,3);
