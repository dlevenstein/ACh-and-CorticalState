basePath = pwd;

[baseFolder,baseName] = fileparts(basePath);
savefile = fullfile(basePath,[baseName,'.LayerID.lfp.mat']);
figfolder = fullfile(basePath,'AnalysisFigures');

sessionInfo = bz_getSessionInfo(basePath, 'noPrompts', true);
channels = sessionInfo.channels;
usechannels = sessionInfo.AnatGrps.Channels;

%% 
mergefile = fullfile(basePath,[baseName,'.MergePoints.events.mat']);
load(mergefile,'MergePoints');

spontidx = find(startsWith(MergePoints.foldernames,"Spont"));
sponttimes = [MergePoints.timestamps(spontidx(1),1) MergePoints.timestamps(spontidx(end),2)];
spontendsample = MergePoints.timestamps_samples(spontidx(end),2)+1;

%%
LOSPEC = []; HISPEC = []; MUAdepth = [];
tLOSPEC = []; tHISPEC = []; 
for ii = 1:length(channels)
    ii
    [lof,lospec,t_lo,hif,hispec,t_hi,MUA] = bz_MUAGammafromDat(basePath,'channels',channels(ii));
    
    tempsidx = find(t_lo < sponttimes(2));
    LOSPEC = cat(2,LOSPEC,mean(lospec(:,tempsidx),2));
    tempsidx = find(t_hi < sponttimes(2));
    HISPEC = cat(2,HISPEC,mean(hispec(:,tempsidx),2));
   
    MUAdepth = cat(2,MUAdepth,MUA.data);
    
    if spontendsample < MergePoints.timestamps_samples(end,end)+1
        tempsidx = find(t_lo > sponttimes(2));
        tLOSPEC = cat(2,tLOSPEC,mean(lospec(:,tempsidx),2));
        tempsidx = find(t_hi > sponttimes(2));
        tHISPEC = cat(2,tHISPEC,mean(hispec(:,tempsidx),2));
    end
end

%% Saving laminar pspectrum
LayerID.LoPSpec = LOSPEC;
LayerID.touchLoPSpec = tLOSPEC;
LayerID.Lofreq = lof;
LayerID.HiPSpec = HISPEC;
LayerID.touchHiPSpec = tHISPEC;
LayerID.Hifreq = hif;
LayerID.MUA.data = MUAdepth;
LayerID.MUA.timestamps = MUA.timestamps;
LayerID.channels = channels;

save(savefile,'LayerID');
