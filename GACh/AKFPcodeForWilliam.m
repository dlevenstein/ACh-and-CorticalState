%% pre-load data structure
%need variables rawFs (sampling rate) and rawFP (photometry signal)

rawFs = 1250;
basePath = pwd;
baseName = bz_BasenameFromBasepath(basePath);
abfname = fullfile(basePath,[baseName,'.abf']);
[abffile,si,file_info] = abfload(abfname);

data = struct;
data.acq.Fs = rawFs;        %sampling rate for photometry signal
data.acq.FPnames = {'ACh'}; 
data.acq.nFPchan = 1; 
data.acq.FP = cell(1,1); 
    data.acq.FP{1} = GACh; %photometry signal
data.acq.time = [1:length(GACh)]/rawFs;

%% open parameters
params_FPonly

%% process photometry signal
data = processFP(data, params);
%processed photometry signal will be in data.final.FP{1}
%corresponding time vector will be data.final.time