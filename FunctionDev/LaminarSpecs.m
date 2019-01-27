basePath = pwd;

[baseFolder,baseName] = fileparts(basePath);
savefile = fullfile(basePath,[baseName,'.LaminarPowerSpec.lfp.mat']);
figfolder = fullfile(basePath,'AnalysisFigures');

savfile = fullfile(basePath,[baseName,'.sessionInfo.mat']);
load(savfile,'sessionInfo');
badchannels = sessionInfo.badchannels;
usechannels = sessionInfo.AnatGrps.Channels;
usechannels(ismember(usechannels,badchannels))=[];

lfp = bz_GetLFP('all','noPrompts',true);

[~,chanidx] = ismember(usechannels,lfp.channels);

clear lfp

%%
lamwavespec = [];
lammua = [];
for ii = 1:length(usechannels)
    [f, pspec, MUA] = bz_MUAGammafromDat(basePath,'channels',chanidx(ii));
    ii
    lamwavespec = cat(2,lamwavespec,pspec);
    lammua = cat(2,lammua,[nanmean(MUA); nanmedian(MUA); max(MUA); min(MUA)]);
end

%% Saving laminar pspectrum
lampspectrum.data = lamwavespec;
lampspectrum.muapow = lammua;
lampspectrum.freq = f;

save(savefile,'lampspectrum')
