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

%% State LFP-LFP cross-corr
SpecXcorr_all = zeros(size(cLayerSpec_all,2),size(cLayerSpec_all,2));
SpecXcorr_all_p = zeros(size(cLayerSpec_all,2),size(cLayerSpec_all,2));
for x = 1:size(cLayerSpec_all,2)
    for y = 1:size(cLayerSpec_all,2)
        [SpecXcorr_all(x,y),SpecXcorr_all_p(x,y)] = corr(cLayerSpec_all(:,x),cLayerSpec_all(:,y),...
            'type','spearman','rows','complete');
    end
end

SpecXcorr_NWh = zeros(size(cLayerSpec_NWh,2),size(cLayerSpec_NWh,2));
SpecXcorr_NWh_p = zeros(size(cLayerSpec_NWh,2),size(cLayerSpec_NWh,2));
for x = 1:size(cLayerSpec_NWh,2)
    for y = 1:size(cLayerSpec_NWh,2)
        [SpecXcorr_NWh(x,y),SpecXcorr_NWh_p(x,y)] = corr(cLayerSpec_NWh(:,x),cLayerSpec_NWh(:,y),...
            'type','spearman','rows','complete');
    end
end

SpecXcorr_Wh = zeros(size(cLayerSpec_Wh,2),size(cLayerSpec_Wh,2));
SpecXcorr_Wh_p = zeros(size(cLayerSpec_Wh,2),size(cLayerSpec_Wh,2));
for x = 1:size(cLayerSpec_Wh,2)
    for y = 1:size(cLayerSpec_Wh,2)
        [SpecXcorr_Wh(x,y),SpecXcorr_Wh_p(x,y)] = corr(cLayerSpec_Wh(:,x),cLayerSpec_Wh(:,y),...
            'type','spearman','rows','complete');
    end
end

SpecXcorr_loP = zeros(size(cLayerSpec_loP,2),size(cLayerSpec_loP,2));
SpecXcorr_loP_p = zeros(size(cLayerSpec_loP,2),size(cLayerSpec_loP,2));
for x = 1:size(cLayerSpec_loP,2)
    for y = 1:size(cLayerSpec_loP,2)
        [SpecXcorr_loP(x,y),SpecXcorr_loP_p(x,y)] = corr(cLayerSpec_loP(:,x),cLayerSpec_loP(:,y),...
            'type','spearman','rows','complete');
    end
end

SpecXcorr_hiP = zeros(size(cLayerSpec_hiP,2),size(cLayerSpec_hiP,2));
SpecXcorr_hiP_p = zeros(size(cLayerSpec_hiP,2),size(cLayerSpec_hiP,2));
for x = 1:size(cLayerSpec_hiP,2)
    for y = 1:size(cLayerSpec_hiP,2)
        [SpecXcorr_hiP(x,y),SpecXcorr_hiP_p(x,y)] = corr(cLayerSpec_hiP(:,x),cLayerSpec_hiP(:,y),...
            'type','spearman','rows','complete');
    end
end

% Saving to struct
LayerSpectral.SpecXcorr_all = SpecXcorr_all;
LayerSpectral.SpecXcorr_all_p = SpecXcorr_all_p;
LayerSpectral.SpecXcorr_NWh = SpecXcorr_NWh;
LayerSpectral.SpecXcorr_NWh_p = SpecXcorr_NWh_p;
LayerSpectral.SpecXcorr_Wh = SpecXcorr_Wh;
LayerSpectral.SpecXcorr_Wh_p = SpecXcorr_Wh_p;
LayerSpectral.SpecXcorr_loP = SpecXcorr_loP;
LayerSpectral.SpecXcorr_loP_p = SpecXcorr_loP_p;
LayerSpectral.SpecXcorr_hiP = SpecXcorr_hiP;
LayerSpectral.SpecXcorr_hiP_p = SpecXcorr_hiP_p;


%% FIGURE:
cmin = min([min(min(SpecXcorr_all(usechannels+1,usechannels+1)))...
    min(min(SpecXcorr_NWh(usechannels+1,usechannels+1)))...
    min(min(SpecXcorr_Wh(usechannels+1,usechannels+1)))]);
cmax = max([max(max(SpecXcorr_all(usechannels+1,usechannels+1)))...
    max(max(SpecXcorr_NWh(usechannels+1,usechannels+1)))...
    max(max(SpecXcorr_Wh(usechannels+1,usechannels+1)))]);

figure;
subplot(1,3,1);
imagesc(normdepth,normdepth,SpecXcorr_all(usechannels+1,usechannels+1));
colormap(gca,'jet')
axis square
ColorbarWithAxis([cmin cmax],['Spearman corr'])
caxis([cmin cmax])
xlabel('normalized depth');ylabel('normalized depth');
title('LFPxcorr ALL');

subplot(1,3,2);
imagesc(normdepth,normdepth,SpecXcorr_NWh(usechannels+1,usechannels+1));
colormap(gca,'jet')
axis square
ColorbarWithAxis([cmin cmax],['Spearman corr'])
caxis([cmin cmax])
xlabel('normalized depth');ylabel('normalized depth');
title('LFPxcorr NWh');

subplot(1,3,3);
imagesc(normdepth,normdepth,SpecXcorr_Wh(usechannels+1,usechannels+1));
colormap(gca,'jet')
axis square
ColorbarWithAxis([cmin cmax],['Spearman corr'])
caxis([cmin cmax])
xlabel('normalized depth');ylabel('normalized depth');
title('LFPxcorr Wh');

NiceSave('LaminarLFPXcorr_Wh_NWh',figfolder,baseName)

%% FIGURE:
cmin = min([min(min(SpecXcorr_all(usechannels+1,usechannels+1)))...
    min(min(SpecXcorr_loP(usechannels+1,usechannels+1)))...
    min(min(SpecXcorr_hiP(usechannels+1,usechannels+1)))]);
cmax = max([max(max(SpecXcorr_all(usechannels+1,usechannels+1)))...
    max(max(SpecXcorr_loP(usechannels+1,usechannels+1)))...
    max(max(SpecXcorr_hiP(usechannels+1,usechannels+1)))]);

figure;
subplot(1,3,1);
imagesc(normdepth,normdepth,SpecXcorr_all(usechannels+1,usechannels+1));
colormap(gca,'jet')
axis square
ColorbarWithAxis([cmin cmax],['Spearman corr'])
caxis([cmin cmax])
xlabel('normalized depth');ylabel('normalized depth');
title('LFPxcorr ALL');

subplot(1,3,2);
imagesc(normdepth,normdepth,SpecXcorr_loP(usechannels+1,usechannels+1));
colormap(gca,'jet')
axis square
ColorbarWithAxis([cmin cmax],['Spearman corr'])
caxis([cmin cmax])
xlabel('normalized depth');ylabel('normalized depth');
title('LFPxcorr <median Pupil area');

subplot(1,3,3);
imagesc(normdepth,normdepth,SpecXcorr_hiP(usechannels+1,usechannels+1));
colormap(gca,'jet')
axis square
ColorbarWithAxis([cmin cmax],['Spearman corr'])
caxis([cmin cmax])
xlabel('normalized depth');ylabel('normalized depth');
title('LFPxcorr >median Pupil area');

NiceSave('LaminarLFPXcorr_lo_hiPup',figfolder,baseName)