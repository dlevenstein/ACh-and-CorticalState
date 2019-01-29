basePath = pwd;

[baseFolder,baseName] = fileparts(basePath);
savefile = fullfile(basePath,[baseName,'.LayerID.lfp.mat']);
figfolder = fullfile(basePath,'AnalysisFigures');

savfile = fullfile(basePath,[baseName,'.sessionInfo.mat']);
load(savfile,'sessionInfo');
channels = sessionInfo.channels;
usechannels = sessionInfo.AnatGrps.Channels;

%%
LOSPEC = []; HISPEC = []; MUAdepth = [];
for ii = 1:length(channels)
    [lof,lospec,hif,hispec,MUA] = bz_MUAGammafromDat(basePath,'channels',channels(ii));
    ii
    LOSPEC = cat(2,LOSPEC,lospec);
    HISPEC = cat(2,HISPEC,hispec);
    MUAdepth = cat(2,MUAdepth,MUA.data);
end

%% Saving laminar pspectrum
LayerID.LoPSpec = LOSPEC;
LayerID.Lofreq = lof;
LayerID.HiPSpec = HISPEC;
LayerID.Hifreq = hif;
LayerID.MUA.data = MUAdepth;
LayerID.MUA.timestamps = MUA.timestamps;

save(savefile,'LayerID');
