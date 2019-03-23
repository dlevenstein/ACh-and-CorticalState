%% Laminar MUA power (0.15 - 2 kHz)
load(fullfile(basePath,[baseName,'.MUA.lfp.mat']));

cLayerMUA_all = nanmean(MUA.data,1);

% for Wh
events = interp1(MUA.timestamps,MUA.timestamps,...
    EMGwhisk.ints.Wh,'nearest');
allidx = [];
for e = 1:size(events,1)
    allidx = cat(1,allidx,wavespec.timestamps(find([wavespec.timestamps >= events(e,1)...
        & wavespec.timestamps <= events(e,2)])));
end
allidx = allidx.*MUA.samplingRate;
cLayerMUA_Wh = nanmean(MUA.data(round(allidx),:),1);
        
% for NWh
events = interp1(MUA.timestamps,MUA.timestamps,...
    EMGwhisk.ints.NWh,'nearest');
allidx = [];
for e = 1:size(events,1)
    allidx = cat(1,allidx,wavespec.timestamps(find([wavespec.timestamps >= events(e,1)...
        & wavespec.timestamps <= events(e,2)])));
end
allidx = allidx.*MUA.samplingRate;
cLayerMUA_NWh = nanmean(MUA.data(round(allidx),:),1);
        
% if ~isempty(Piezotouch)
%     % for Touch
%     events = interp1(MUA.timestamps,MUA.timestamps,...
%         Piezotouch.ints.Touch,'nearest');
%     allidx = [];
%     for e = 1:size(events,1)
%         allidx = cat(1,allidx,wavespec.timestamps(find([wavespec.timestamps >= events(e,1)...
%             & wavespec.timestamps <= events(e,2)])));
%     end
%     allidx = allidx.*MUA.samplingRate;
%     cLayerMUA_T = nanmean(MUA.data(round(allidx),:),1);
% else
% end

events = interp1(MUA.timestamps,MUA.timestamps,...
    lowpupildata.timestamps(highpup),'nearest').*MUA.samplingRate;
cLayerMUA_hiP = nanmean(MUA.data(round(events)),1);

events = interp1(MUA.timestamps,MUA.timestamps,...
    lowpupildata.timestamps(~highpup),'nearest').*MUA.samplingRate;
cLayerMUA_loP = nanmean(MUA.data(round(events)),1);

% event MUA
eventMUA_Wh = eventMUA(MUA,EMGwhisk.ints.Wh(:,1),...
    'channels',channels,'twin',[0.75 0.75],'spat_sm',0,'saveMat',false);

% Saving to struct
LayerSpectral.cLayerMUA_all = cLayerMUA_all;
LayerSpectral.cLayerMUA_Wh = cLayerMUA_Wh;
LayerSpectral.cLayerMUA_NWh = cLayerMUA_NWh;
LayerSpectral.cLayerMUA_hiP = cLayerMUA_hiP;
LayerSpectral.cLayerMUA_loP = cLayerMUA_loP;
LayerSpectral.eventMUA_Wh = eventMUA_Wh;

%% Laminar eventCSD/MUA
% sort thourgh high/low Wh/pupil/SW amplitude events
% adjust by phase of SlowWaves
% easyyy lines but specifying better intervals, trial selection
% establish and z score to baseline

% eventCSD = eventCSD (lfp, events, varargin);
% eventMUA = eventMUA (mua, events, varargin);
% eventSpec = eventSpec (spec, events, varargin);

%% Laminar CSD-CSD/MUA-MUA cross-corr
% thin about this...
