function [ eventCSD ] = eventCSD (lfp, events, varargin)

% [ eventCSD ] = bz_eventCSD (lfp, events, varargin)
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
addParameter(p,'plotCSD',true,@islogical);
addParameter(p,'plotLFP',true,@islogical);
addParameter(p,'cwin',[]);
addParameter(p,'saveMat',true,@islogical)
addParameter(p,'saveName','eventCSD')

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
if isempty(lfp)
    lfp = bz_GetLFP('all','noPrompts',true);
end

%lfp input
if isstruct(lfp)
    data = lfp.data;
    timestamps = lfp.timestamps;
    samplingRate = lfp.samplingRate;
elseif iscell(lfp) %for multiple trials
    celllengths = cellfun(@length,lfp);
    data = vertcat(lfp{:});
elseif isnumeric(lfp)
    data = lfp;
    timestamps = [1:length(lfp)]'./samplingRate;
end

twin = p.Results.twin*samplingRate;
events = round(events*samplingRate);

%% Conpute event-triggered LFP average

lfp_temp = nan(twin(1)+twin(2)+1,length(channels),length(events));
csd_temp = nan(twin(1)+twin(2)+1,length(channels)-2,length(events));

[~,chanidx] = ismember(channels,lfp.channels);
for e = 1:length(events)
    if events(e)-twin(1) > 0 && events(e)+twin(2) < size(data,1)
        temp = data(events(e)-twin(1):events(e)+twin(2),chanidx)*-1;
        lfp_temp(:,:,e) = temp; 
        csd_temp(:,:,e) = diff(temp,2,2);
    else
    end
end

lfp_avg = nanmean(lfp_temp,3);
CSD = nanmean(csd_temp,3);

%% Compute CSD

% detrend
if doDetrend
    lfp_avg = detrend(lfp_avg')';
    CSD = detrend(CSD')';
end

% temporal smoothing
if temp_sm > 0
    for ch = 1:size(lfp_avg,2)
        lfp_avg(:,ch) = smooth(lfp_avg(:,ch),temp_sm,'sgolay');
    end
    for ch = 1:size(CSD,2)
        CSD(:,ch) = smooth(CSD(:,ch),temp_sm,'sgolay');
    end
end

% spatial smoothing
if spat_sm > 0
    for t = 1:size(lfp_avg,1)
        lfp_avg(t,:) = smooth(lfp_avg(t,:),spat_sm,'lowess');
    end
    for t = 1:size(CSD,1)
        CSD(t,:) = smooth(CSD(t,:),spat_sm,'lowess');
    end
end

% calculate CSD
%CSD = diff(lfp_avg,2,2);

% generate output structure
eventCSD.CSDdata = CSD;
eventCSD.LFPdata = lfp_avg;
eventCSD.timestamps = -twin(1):twin(2);
eventCSD.samplingRate = samplingRate;
eventCSD.channels = channels;
eventCSD.params.spat_sm = spat_sm;
eventCSD.params.temp_sm = temp_sm;
eventCSD.params.detrend = doDetrend;

if saveMat
    save(savefile,'eventCSD')
end

%% Plot
taxis = (-(twin(1)/samplingRate):(1/samplingRate):(twin(2)/samplingRate))*1e3;
cmax = max(max(CSD));

if ~isempty(cwin)
    cmax = cwin;
end

if plotLFP
    
    figure;
    subplot(1,5,1:3);
    contourf(taxis,1:size(CSD,2),CSD',40,'LineColor','none');hold on;
    colormap jet; caxis([-cmax cmax]); 
    c = colorbar;
    c.Label.String = 'sink -> source';
    ylim([1 size(CSD,2)]);
    set(gca,'YDir','reverse');xlabel('time (ms)');ylabel('channel');title('CSD');
    set(gca,'Ytick',[1:1:size(CSD,2)]);
    set(gca,'Yticklabels',channels(2:end-1));
    set(gca,'YGrid','on','Layer','top','GridColor',[0 0 0]);
    plot([0 0],[1 size(CSD,2)],'--k');hold on;
    
    subplot(1,5,4:5);
    for ch=1:size(lfp_avg,2)
        offset = 300*(ch-1);
        sh_tmp = 2.5.*(lfp_avg(:,ch)) + offset;
        plot(taxis,sh_tmp,'k','LineWidth',1); hold on;
        clear sh_tmp
    end
    set(gca,'YDir','reverse','YTickLabel',[]);
    ylim([-1500 offset+1500]);xlim([taxis(1) taxis(end)]);
    xlabel('time (ms)'); title('average LFP');
    plot([0 0],ylim,'--r');hold on;
    
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
