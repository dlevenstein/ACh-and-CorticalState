function [ eventMUA ] = eventMUA (mua, events, varargin)

% [ eventMUA ] = bz_eventMUA (lfp, events, varargin)
% Calculates event-triggered (i.e. SWRs) CSD map from a linear array of LFPs
%
% INPUT
%    lfp            a buzcode structure with fields lfp.data,
%                                                   lfp.timestamps
%                                                   lfp.samplingRate
%                   -lfp can also be a [t x 1] timeseries signal. in which
%                   case you need to input 'samplingRate'
%    events          events timestamps (in sec)
%
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%       channels    vector with channels to inlcude. (default 1:64)
%       twin        time window around events to calculate average. Default [0.1 0.1]
%       spat_sm     degree of spatial smoothing. Default = 0.
%       temp_sm     degree of temporal smoothing. Default = 0.
%       doDetrend   Default false.
%       plotCSD     true/false. Default true.
%       plotLFP     true/false. Default true.
%       saveName    Default eventCSD.
%    =========================================================================
%
% OUTPUT:
%    eventCSD            a buzcode structure with fields eventcsd.data,
%                                                        eventcsd.timestamps
%                                                        eventcsd.samplingRate
%                                                        eventcsd.channels
%                                                        eventcsd.params
%   
%
% Antonio FR, Levenstein D, Munoz W - 7/18
%
% TODO: Make so it can read from the binary lfp file in chunks instead of from a lfp.mat
%
%% Parse inputs

p = inputParser;
addParameter(p,'channels',[1:64],@isvector);
addParameter(p,'samplingRate',1250,@isnumeric);
addParameter(p,'twin',[0.1 0.1],@isnumeric);
addParameter(p,'spat_sm',0,@isnumeric);
addParameter(p,'temp_sm',0,@isnumeric);
addParameter(p,'doDetrend',false,@islogical);
addParameter(p,'plotCSD',false,@islogical);
addParameter(p,'plotLFP',false,@islogical);
addParameter(p,'cwin',[]);
addParameter(p,'saveMat',true,@islogical)
addParameter(p,'saveName','eventMUA')

parse(p,varargin{:});
channels = p.Results.channels;
samplingRate = p.Results.samplingRate;
spat_sm = p.Results.spat_sm;
temp_sm = p.Results.temp_sm;
doDetrend = p.Results.doDetrend;
plotCSD = p.Results.plotCSD;
plotLFP = p.Results.plotLFP;
cwin = p.Results.cwin;
saveMat = p.Results.saveMat;
saveName = p.Results.saveName;

%Defaults
if ~exist('basePath','var')
    basePath = pwd;
end
[baseFolder,baseName] = fileparts(basePath);

%% File Management
baseName = bz_BasenameFromBasepath(basePath);
figfolder = fullfile(basePath,'DetectionFigures');

savefile = fullfile(basePath,[baseName,'.',saveName,'.lfp.mat']);

%% Load the LFP
if isempty(mua)
    load(fullfile(basePath,[baseName,'.MUA.lfp.mat']));
    mua = MUA;
end

%lfp input
if isstruct(mua)
    data = mua.data;
    timestamps = mua.timestamps;
    samplingRate = mua.samplingRate;
elseif iscell(mua) %for multiple trials
    celllengths = cellfun(@length,mua);
    data = vertcat(mua{:});
elseif isnumeric(mua)
    data = mua;
    timestamps = [1:length(mua)]'./samplingRate;
end

twin = p.Results.twin*samplingRate;
events = round(events*samplingRate);

%% Conpute event-triggered LFP average

mua_temp = nan(twin(1)+twin(2)+1,length(channels),length(events));

[~,chanidx] = ismember(channels,mua.channels);
for e = 1:length(events)
    if events(e)-twin(1) > 0 && events(e)+twin(2) < size(data,1)
        temp = data(events(e)-twin(1):events(e)+twin(2),chanidx);
        mua_temp(:,:,e) = temp;
    else
    end
end

mua_avg = nanmean(mua_temp,3);

%% Compute CSD

% detrend
if doDetrend
    mua_avg = detrend(mua_avg')';
end

% temporal smoothing
if temp_sm > 0
    for ch = 1:size(mua_avg,2)
        mua_avg(:,ch) = smooth(mua_avg(:,ch),temp_sm,'sgolay');
    end
end

% spatial smoothing
if spat_sm > 0
    for t = 1:size(mua_avg,1)
        mua_avg(t,:) = smooth(mua_avg(t,:),spat_sm,'lowess');
    end
end

% calculate CSD
%CSD = diff(lfp_avg,2,2);

% generate output structure
eventMUA.data = mua_avg;
eventMUA.timestamps = -twin(1):twin(2);
eventMUA.samplingRate = samplingRate;
eventMUA.channels = channels;
eventMUA.params.spat_sm = spat_sm;
eventMUA.params.temp_sm = temp_sm;
eventMUA.params.detrend = doDetrend;

if saveMat
    save(savefile,'eventMUA')
end

%% Plot
taxis = (-(twin(1)/samplingRate):(1/samplingRate):(twin(2)/samplingRate))*1e3;
cmax = max(max(mua_avg));

if ~isempty(cwin)
    cmax = cwin;
end

if plotLFP
    
    figure;
    contourf(taxis,1:size(mua_avg,2),mua_avg',40,'LineColor','none');hold on;
    colormap jet; caxis([0 cmax]); 
    c = colorbar;
    c.Label.String = 'MUA power (0.5-5 KHz)';
    ylim([1 size(mua_avg,2)]);
    set(gca,'YDir','reverse');xlabel('time (ms)');ylabel('channel');title('MUA power');
    set(gca,'Ytick',[1:1:size(mua_avg,2)]);
    set(gca,'Yticklabels',channels);
    set(gca,'YGrid','on','Layer','top','GridColor',[0 0 0]);
    plot([0 0],[1 size(mua_avg,2)],'--k');hold on;
        
elseif plotCSD
    
    figure;
    subplot(1,2,1);
    contourf(taxis,1:size(CSD,2),CSD',40,'LineColor','none');hold on;
    colormap jet; caxis([-cmax cmax]);
    set(gca,'YDir','reverse');xlabel('time (ms)');ylabel('channel');title('CSD');
    plot([0 0],[1 size(CSD,2)],'--k');hold on;
    
end

NiceSave(saveName,figfolder,baseName)

end
