function [rescaled] = rescaleCx( basePath, varargin )
%
%   RESCALES laminar .lfp data according to normalized cortical column, as in
%   Munoz et al. 2017 Science Sst IN paper. Assumes H3 Cambridge NeuroTech
%   probe.
%
%   INPUTS
%    'basePath'           - folder in which .lfp.mat file will be found (default
%                           is pwd)
%                           folder should follow buzcode standard:
%                           whateverPath/baseName
%                           and contain file baseName.lfp
%    'BADOUT'             - to indicate if bad channels had to be excluded 
%                           for analysis, as in laminar CSD analysis (default: true)
%   OUTPUT
%     rescaled       struct of lfp data. 
%    .depth          [1 x Nd] vector of absolute distance to match LFP data
%    .ndepth         [1 x Nd] vector of normalized depth (0-1) to match LFP data
%    .channels       [1 x Nd] vector of channel ID's
%    .lnames         cell array of layernames
%    .layer          layer assigned to each channel
%
% WMunoz - 02/28/2019
%
%% Parse the inputs!
% parse args
p = inputParser;
addParameter(p,'BADOUT',true,@islogical);

parse(p,varargin{:})
BADOUT = p.Results.BADOUT;

%% Scales according to Cambridge NeuroTech H3 probe
chandist = [0:20:1275];
%lnorm = [0 0.1 0.35 0.5 0.6 0.75 0.9 1]; %including 6a-6b division
%lnorm = [0 0.1 0.35 0.5 0.6 0.75 1]; %proper 5b and 6 layers
lnorm = [0 0.1 0.35 0.5 0.6 0.9 1]; %tentative given lack of 5b-6 boundary
lnames = {'L1','L23','L4','L5a','L5b6','L6'};

%% Loading data
%[baseFolder,baseName] = fileparts(basePath);
%data = load(fullfile(basePath,[baseName,'.',analysisName,'.lfp.mat']));

sessionInfo = bz_getSessionInfo(basePath,'noPrompts',true);
usechannels = sessionInfo.AnatGrps.Channels;
%lborders = [NaN 20 18 62 16 36 33 32]; %DEV
lborders = sessionInfo.layerborders;

%% Rescaling...
truelayer1 = chandist(usechannels == lborders(end))*lnorm(2);
chandist = chandist-(chandist(usechannels == lborders(2))-truelayer1); 
truecolumn = chandist(usechannels == lborders(end));
normdepth = NaN(1,length(usechannels));
layer = cell(size(normdepth));

for i = 2:length(lnorm)-1
    lb1 = find(usechannels == lborders(i));
    lb2 = find(usechannels == lborders(i+1));
    lfactor1 = lnorm(i);
    lfactor2 = lnorm(i+1);
    normdepth(lb1:lb2) = ((chandist(lb1:lb2) - chandist(lb1)).* (lfactor2 - lfactor1))./(chandist(lb2) - chandist(lb1)) + lfactor1;
end

%Dealing w/ pia:L1
if chandist(1) < 0
    lb1 = find(chandist >= 0);
else
    lb1 = 1;
end
lb2 = find(usechannels == lborders(2));
lfactor1 = lnorm(1);
lfactor2 = lnorm(2);
temp = (([0 chandist(lb1(1):lb2)] - 0).* (lfactor2 - lfactor1))./(chandist(lb2) - 0) + lfactor1;
normdepth(lb1(1):lb2) = temp(2:end);

% Tagging bad channels, if needed
if BADOUT
    badchannels = sessionInfo.badchannels;
    badidx = ismember(usechannels,badchannels);
    usechannels(badidx) = NaN;
    chandist(badidx) = NaN;
    normdepth(badidx) = NaN;
end


for ii = 1:length(lnames)
    inlayerchans = normdepth>=lnorm(ii) & normdepth<lnorm(ii+1);
    layer(inlayerchans) = lnames(ii);
end

%% Rescaled data to struct
rescaled.depth = chandist;
rescaled.ndepth = normdepth;
rescaled.channels = usechannels;
rescaled.layer = layer;
rescaled.lnames = lnames;
rescaled.boundaries = lnorm;