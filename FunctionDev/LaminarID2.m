basePath = pwd;

[baseFolder,baseName] = fileparts(basePath);
savefile = fullfile(basePath,[baseName,'.',num2str(i),'.LayerID.lfp.mat']);

[lof,lospec,t_lo,hif,hispec,t_hi,MUA] = MUAGammafromDat(basePath,'channels',i);

%% Saving laminar pspectrum
LayerID.lof = lof;
LayerID.lospec = lospec;
LayerID.t_lo = t_lo;
LayerID.hif = hif;
LayerID.hispec = hispec;
LayerID.t_hi = t_hi;
LayerID.MUA = MUA;

save(savefile,'LayerID');