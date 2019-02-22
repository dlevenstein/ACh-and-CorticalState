function [ spec ] = eventSpec (events, varargin)

% [ spec ] = eventSpec (events, varargin)
% Calculates event-triggered (i.e. SWRs) spectrogram map from a linear array of LFPs
%
% INPUT
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
%    eventSpec            a buzcode structure with fields eventcsd.data,
%                                                         eventcsd.timestamps
%                                                         eventcsd.samplingRate
%                                                         eventcsd.channels
%                                                         eventcsd.params
%   
%
% Munoz W - 2/19

%% Parse inputs
p = inputParser;
addParameter(p,'channels',[1:64],@isvector);
addParameter(p,'twin',[0.1 0.1],@isnumeric);
addParameter(p,'spat_sm',0,@isnumeric);
addParameter(p,'temp_sm',0,@isnumeric);
addParameter(p,'doDetrend',false,@islogical);
addParameter(p,'plotCSD',true,@islogical);
addParameter(p,'plotLFP',true,@islogical);
addParameter(p,'cwin',[]);
addParameter(p,'lfp',[]);
addParameter(p,'saveMat',true,@islogical)
addParameter(p,'saveName','eventCSD')

parse(p,varargin{:});
lfp = p.Results.lfp;
channels = p.Results.channels;
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

lotwin = round(p.Results.twin*(1/4));
loevents = round(events*(1/4));
hitwin = round(p.Results.twin*(1/0.02));
hievents = round(events*(1/0.02));

%% Conpute event-triggered spectrogram
LOspec = []; HIspec = [];

[~,chanidx] = ismember(channels,lfp.channels);
clear lfp

for i = chanidx
    
    %%
    load(fullfile(basePath,[baseName,'.',num2str(i),'.LayerID.lfp.mat']),'LayerID');
    
    lof = LayerID.lof;
    lospec = LayerID.lospec;
    t_lo = LayerID.t_lo;
    
    hif = LayerID.hif;
    hispec = LayerID.hispec;
    t_hi = LayerID.t_hi;
    
    %%
    lospec_temp = nan(lotwin(1)+lotwin(2)+1,size(lospec,1),length(loevents));
    hispec_temp = nan(hitwin(1)+hitwin(2)+1,size(hispec,1),length(hievents));
     
    %%
    
    for e = 1:length(loevents)
        if loevents(e)-lotwin(1) > 0 && loevents(e)+lotwin(2) < size(lospec,2)
            lospec_temp(:,:,e) = lospec(:,loevents(e)-lotwin(1):loevents(e)+lotwin(2))';
        else
        end
    end
    
    for e = 1:length(hievents)
        if hievents(e)-hitwin(1) > 0 && hievents(e)+hitwin(2) < size(hispec,2)
            hispec_temp(:,:,e) = hispec(:,hievents(e)-hitwin(1):hievents(e)+hitwin(2))';
        else
        end
    end
    
    %%
    LOspec = cat(LOspec,nanmean(lospec_temp,3),3); 
    HIspec = cat(HIspec,nanmean(hispec_temp,3),3); 
    
end

lospec = nanmean(LOspec,3);
hispec = nanmean(HIspec,3);

%% Smoothening
% detrend
% if doDetrend
%     lfp_avg = detrend(lfp_avg')';
%     CSD = detrend(CSD')';
% end

% temporal smoothing
% if temp_sm > 0
%     for ch = 1:size(lfp_avg,2)
%         lfp_avg(:,ch) = smooth(lfp_avg(:,ch),temp_sm,'sgolay');
%     end
%     for ch = 1:size(CSD,2)
%         CSD(:,ch) = smooth(CSD(:,ch),temp_sm,'sgolay');
%     end
% end

% spatial smoothing
% if spat_sm > 0
%     for t = 1:size(lfp_avg,1)
%         lfp_avg(t,:) = smooth(lfp_avg(t,:),spat_sm,'lowess');
%     end
%     for t = 1:size(CSD,1)
%         CSD(t,:) = smooth(CSD(t,:),spat_sm,'lowess');
%     end
% end

%% Generate output structure
% eventCSD.CSDdata = CSD;
% eventCSD.LFPdata = lfp_avg;
% eventCSD.timestamps = -twin(1):twin(2);
% eventCSD.samplingRate = samplingRate;
% eventCSD.channels = channels;
% eventCSD.params.spat_sm = spat_sm;
% eventCSD.params.temp_sm = temp_sm;
% eventCSD.params.detrend = doDetrend;
% 
% if saveMat
%     save(savefile,'eventCSD')
% end

%% Plot
% taxis = (-(twin(1)/samplingRate):(1/samplingRate):(twin(2)/samplingRate))*1e3;
% cmax = max(max(CSD));
% 
% if ~isempty(cwin)
%     cmax = cwin;
% end
% 
% if plotLFP
%     
%     figure;
%     subplot(1,5,1:3);
%     contourf(taxis,1:size(CSD,2),CSD',40,'LineColor','none');hold on;
%     colormap jet; caxis([-cmax cmax]); 
%     c = colorbar;
%     c.Label.String = 'sink -> source';
%     ylim([1 size(CSD,2)]);
%     set(gca,'YDir','reverse');xlabel('time (ms)');ylabel('channel');title('CSD');
%     set(gca,'Ytick',[1:1:size(CSD,2)]);
%     set(gca,'Yticklabels',channels(2:end-1));
%     set(gca,'YGrid','on','Layer','top','GridColor',[0 0 0]);
%     plot([0 0],[1 size(CSD,2)],'--k');hold on;
%     
%     subplot(1,5,4:5);
%     for ch=1:size(lfp_avg,2)
%         offset = 300*(ch-1);
%         sh_tmp = 2.5.*(lfp_avg(:,ch)) + offset;
%         plot(taxis,sh_tmp,'k','LineWidth',1); hold on;
%         clear sh_tmp
%     end
%     set(gca,'YDir','reverse','YTickLabel',[]);
%     ylim([-1500 offset+1500]);xlim([taxis(1) taxis(end)]);
%     xlabel('time (ms)'); title('average LFP');
%     plot([0 0],ylim,'--r');hold on;
%     
% elseif plotCSD
%     
%     figure;
%     subplot(1,2,1);
%     contourf(taxis,1:size(CSD,2),CSD',40,'LineColor','none');hold on;
%     colormap jet; caxis([-cmax cmax]);
%     set(gca,'YDir','reverse');xlabel('time (ms)');ylabel('channel');title('CSD');
%     plot([0 0],[1 size(CSD,2)],'--k');hold on;
%     
% end
% 
% NiceSave(saveName,figfolder,baseName)

end