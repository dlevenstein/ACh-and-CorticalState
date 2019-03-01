basePath = pwd;
baseName = bz_BasenameFromBasepath(basePath);
figfolder = fullfile(basePath,'PSSstructs');

%%
lowerbound = [1:20];
upperbound = [60:120];

%%

grandWhiskPSScorr.EMG = zeros(length(lowerbound),length(upperbound));
grandWhiskPSScorr.pup = zeros(length(lowerbound),length(upperbound));
grandWhiskPSScorr.dpdt = zeros(length(lowerbound),length(upperbound));
grandWhiskPSScorr.phasecoupling = zeros(length(lowerbound),length(upperbound));

for f = 1:length(lowerbound)
    for ff = 1:length(upperbound)
        load(fullfile(figfolder,[baseName,'.',num2str(f),'.',num2str(ff),'.WhiskPSScorr.lfp.mat']),'WhiskPSScorr');
        
        grandWhiskPSScorr.EMG(f,ff) = max(WhiskPSScorr.EMG);
        grandWhiskPSScorr.pup(f,ff) = max(WhiskPSScorr.pup);
        grandWhiskPSScorr.dpdt(f,ff) = max(WhiskPSScorr.dpdt);
        grandWhiskPSScorr.phasecoupling(f,ff) = max(WhiskPSScorr.phasecoupling);
        
    end
end

save(fullfile(basePath,[baseName,'win2.dt0.5.grandWhiskPSScorr.lfp.mat']),'grandWhiskPSScorr');