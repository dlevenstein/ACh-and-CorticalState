function [ ] = BehaviorAChAnalysis(basePath,figfolder)

%Initiate Paths
reporoot = '/home/dlevenstein/ProjectRepos/ACh-and-CorticalState/';
%reporoot = '/Users/dlevenstein/Project Repos/ACh-and-CorticalState/';
%basePath = '/mnt/proraidDL/Database/WMData/AChPupil/EM1M3/171209_WT_EM1M3/';
%basePath = '/mnt/proraidDL/Database/WMData/AChPupil/180706_WT_EM1M3/';
%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/BehaviorAChAnalysis'];


%%
baseName = bz_BasenameFromBasepath(basePath);
%sessionInfo = bz_getSessionInfo(basePath,'noPrompts',true);


%figfolder = fullfile(basePath,'AnalysisFigures');
%savefile = fullfile(basePath,[baseName,'.BehaviorAnalysis.mat']);
%% Loading behavior...
% Pupil diameter
[ pupilcycle ] = ExtractPupilCycle( basePath);
GACh = bz_LoadBehavior(basePath,'GACh');

% %pupildilation = bz_LoadBehavior(basePath,'pupildiameter');
% 
% smoothwin = 0.5; %s
% pupilcycle.data = smooth(pupilcycle.data,smoothwin.*pupilcycle.samplingRate,'moving');
% nantimes = isnan(pupilcycle.data);
% pupilcycle.data = pupilcycle.data(~isnan(pupilcycle.data));
% 
% if length(pupilcycle.data) < 1
%     warning('Not enough pupil data >)');
%     return
% end
% 
% smoothwin = 2; %s
% pupilcycle.dpdt = diff(smooth(pupilcycle.data,smoothwin.*pupilcycle.samplingRate,'moving')).*pupilcycle.samplingRate;
% pupilcycle.dpdt = smooth(pupilcycle.dpdt,smoothwin.*pupilcycle.samplingRate,'moving');
% pupilcycle.dpdt = [pupilcycle.dpdt; nan]; %To align to timestamps
% pupilcycle.timestamps = pupilcycle.timestamps(~nantimes);
% 
% % Filtered Pupil
% lowfilter = [0.01 0.1]; %old. order 3
% lowfilter = [0.02 0.2]; %new: EMG coupled. order 1
% lowfilter = [0.033 0.33]; %new: dur coupled. order 2
% %highfilter = [0.3 0.8];
% 
% pupil4filter = pupilcycle;
% pupilcycle = bz_Filter(pupil4filter,'passband',lowfilter,'filter' ,'fir1','order',2);
% %highpupildata = bz_Filter(pupil4filter,'passband',highfilter,'filter' ,'fir1');
% pupilcycle.pupthresh = -0.8;
% pupilcycle.highpup = log10(pupilcycle.amp)>pupilcycle.pupthresh; 

% EMG
EMGwhisk = bz_LoadBehavior(basePath,'EMGwhisk');


%Specifying SPONT whisking/ACH
load(fullfile(basePath,[baseName,'.MergePoints.events.mat']),'MergePoints');
sidx = find(startsWith(MergePoints.foldernames,"Baseline"));
sponttimes = [MergePoints.timestamps(sidx(1),1) MergePoints.timestamps(sidx(end),2)];

spontidx = find(EMGwhisk.ints.Wh(:,2) < sponttimes(2));
EMGwhisk.ints.Wh = EMGwhisk.ints.Wh(spontidx,:);

spontidx = find(EMGwhisk.ints.NWh(:,2) < sponttimes(2));
EMGwhisk.ints.NWh = EMGwhisk.ints.NWh(spontidx,:);

spontidx = find(EMGwhisk.timestamps < sponttimes(2));
EMGwhisk.timestamps = EMGwhisk.timestamps(spontidx);
EMGwhisk.EMGenvelope = EMGwhisk.EMGenvelope(spontidx);
EMGwhisk.EMG = EMGwhisk.EMG(spontidx);
EMGwhisk.EMGsm = EMGwhisk.EMGsm(spontidx);

EMGwhisk.iswhisk = InIntervals(EMGwhisk.timestamps,EMGwhisk.ints.Wh);

spontidx = find(pupilcycle.timestamps < sponttimes(2));
pupilcycle.amp(spontidx);
pupilcycle.data(spontidx);
pupilcycle.dpdt(spontidx);
pupilcycle.phase(spontidx);
pupilcycle.states(spontidx);
pupilcycle.timestamps(spontidx);

spontidx = find(pupilcycle.ints.highpupstate(:,2) < sponttimes(2));
pupilcycle.ints.highpupstate = pupilcycle.ints.highpupstate(spontidx,:);

spontidx = find(pupilcycle.ints.lowpupstate(:,2) < sponttimes(2));
pupilcycle.ints.lowpupstate = pupilcycle.ints.lowpupstate(spontidx,:);

spontidx = find(GACh.timestamps < sponttimes(2));
GACh.raw = GACh.raw(spontidx);
GACh.timestamps = GACh.timestamps(spontidx);

GACh.raw = GACh.raw./nanmedian(GACh.raw);
%% Get rid of recording start/stop artifact

maxtimejump = 1; %s
pupilcycle.amp = NanPadJumps( pupilcycle.timestamps,pupilcycle.amp,maxtimejump );
pupilcycle.phase = NanPadJumps( pupilcycle.timestamps,pupilcycle.phase,maxtimejump );
pupilcycle.data = NanPadJumps( pupilcycle.timestamps,pupilcycle.data,maxtimejump );
pupilcycle.dpdt = NanPadJumps( pupilcycle.timestamps,pupilcycle.dpdt,maxtimejump );
EMGwhisk.EMGsm = NanPadJumps( EMGwhisk.timestamps,EMGwhisk.EMGsm,maxtimejump );

GACh.raw = NanPadJumps( GACh.timestamps,GACh.raw,maxtimejump );



%% Pupil phase/amp and duration at whisks

EMGwhisk.whisks.pupphase = interp1(pupilcycle.timestamps,pupilcycle.phase,EMGwhisk.ints.Wh(:,1),'nearest');
EMGwhisk.whisks.pupamp = interp1(pupilcycle.timestamps,pupilcycle.amp,EMGwhisk.ints.Wh(:,1),'nearest');
EMGwhisk.whisks.pup = interp1(pupilcycle.timestamps,pupilcycle.data,EMGwhisk.ints.Wh(:,1),'nearest');
EMGwhisk.whisks.dur = diff(EMGwhisk.ints.Wh,[],2);
EMGwhisk.whisks.hipup = InIntervals(EMGwhisk.ints.Wh(:,1),pupilcycle.ints.highpupstate);
EMGwhisk.whisks.lopup = InIntervals(EMGwhisk.ints.Wh(:,1),pupilcycle.ints.lowpupstate);

EMGwhisk.pupphase = interp1(pupilcycle.timestamps,pupilcycle.phase,EMGwhisk.timestamps,'nearest');
EMGwhisk.pupamp = interp1(pupilcycle.timestamps,pupilcycle.amp,EMGwhisk.timestamps,'nearest');
EMGwhisk.pup = interp1(pupilcycle.timestamps,pupilcycle.data,EMGwhisk.timestamps,'nearest');
EMGwhisk.dpdt = interp1(pupilcycle.timestamps(1:end-1),pupilcycle.dpdt,EMGwhisk.timestamps,'nearest');
EMGwhisk.hipup = InIntervals(EMGwhisk.timestamps,pupilcycle.ints.highpupstate);
EMGwhisk.lopup = InIntervals(EMGwhisk.timestamps,pupilcycle.ints.lowpupstate);

GACh.pupphase = interp1(pupilcycle.timestamps,pupilcycle.phase,GACh.timestamps,'nearest');
GACh.pupamp = interp1(pupilcycle.timestamps,pupilcycle.amp,GACh.timestamps,'nearest');
GACh.pup = interp1(pupilcycle.timestamps,pupilcycle.data,GACh.timestamps,'nearest');
GACh.dpdt = interp1(pupilcycle.timestamps(1:end-1),pupilcycle.dpdt,GACh.timestamps,'nearest');
GACh.hipup = InIntervals(GACh.timestamps,pupilcycle.ints.highpupstate);
GACh.lopup = InIntervals(GACh.timestamps,pupilcycle.ints.lowpupstate);
GACh.Wh = InIntervals(GACh.timestamps,EMGwhisk.ints.Wh);
GACh.NWh = InIntervals(GACh.timestamps,EMGwhisk.ints.NWh);
GACh.EMG = interp1(EMGwhisk.timestamps,EMGwhisk.EMGsm,GACh.timestamps,'nearest');

HILO = {'lopup','hipup'};
WHNWH = {'Wh','NWh'};
ONOFF = {'WhOn','WhOFF'};
LONGSHORT = {'long','short'};
%% GACh relative to whsking time
EMGwhisk.dur = diff(EMGwhisk.ints.Wh,[],2);
window = 5; %s
durthresh = 1; %
GACh.whtime.WhOn = nan(size(GACh.timestamps));
GACh.whtime.WhOFF = nan(size(GACh.timestamps));
GACh.long.WhOn = false(size(GACh.timestamps));
GACh.long.WhOFF = false(size(GACh.timestamps));
GACh.short.WhOn = false(size(GACh.timestamps));
GACh.short.WhOFF = false(size(GACh.timestamps));
for wh = 2:(size(EMGwhisk.ints.Wh,1)-1)
    bz_Counter(wh-1,(size(EMGwhisk.ints.Wh,1)-2),'Whisk');
    longwhisk = EMGwhisk.dur(wh)>durthresh;
    
    for oo = 1:2
        if oo == 1
            inwintimes = GACh.timestamps > EMGwhisk.ints.Wh(wh-1,2) & ...
                GACh.timestamps <= EMGwhisk.ints.Wh(wh,2);  
        elseif oo==2
            inwintimes = GACh.timestamps > EMGwhisk.ints.Wh(wh,1) & ...
                GACh.timestamps <= EMGwhisk.ints.Wh(wh+1,1); 
        end
        GACh.whtime.(ONOFF{oo})(inwintimes) = GACh.timestamps(inwintimes)-EMGwhisk.ints.Wh(wh,oo);
        GACh.long.(ONOFF{oo})(inwintimes) = longwhisk;
        GACh.short.(ONOFF{oo})(inwintimes) = ~longwhisk;
    end
    
end
%% Mean GACh in pupil space

[pupilphaseGACh.meanZ,pupilphaseGACh.N,pupilphaseGACh.Xbins,pupilphaseGACh.Ybins] = ...
    ConditionalHist3( GACh.pupphase,...
    log10(GACh.pupamp),GACh.raw,...
    'minXY',500,'Xbounds',[-pi pi],'Ybounds',[-2.25 0.5],...
    'numXbins',30,'numYbins',30);

[pupildpGACh.meanZ,pupildpGACh.N,pupildpGACh.Xbins,pupildpGACh.Ybins] = ...
    ConditionalHist3( log10(GACh.pup),...
    (GACh.dpdt),(GACh.raw),...
    'minXY',250,'Xbounds',[-0.5 0.5],'Ybounds',[-0.5 0.5],...
    'numXbins',40,'numYbins',40);




%% GACh and pupil 
GAChrange = [0.85 1.2];
for pp= 1:2
    GAChdist.(HILO{pp})  = ConditionalHist(...
        GACh.pupphase(GACh.(HILO{pp})),...
        (GACh.raw(GACh.(HILO{pp}))),...
        'Ybounds',GAChrange,'numYbins',80,'numXbins',20,'Xbounds',[-pi pi]);
    
    for ww = 1:2
        GAChdist.(HILO{pp}).(WHNWH{ww})  = ...
            ConditionalHist( GACh.pupphase(GACh.(WHNWH{ww})&GACh.(HILO{pp})),...
                (GACh.raw(GACh.(WHNWH{ww})&GACh.(HILO{pp}))),...
                'Ybounds',GAChrange,'numYbins',80,'numXbins',20,'Xbounds',[-pi pi]);
    end
end

GAChdist.pup  = ConditionalHist( log10(GACh.pup),...
        (GACh.raw),...
        'numXbins',20,'Xbounds',[-0.25 0.25],'Ybounds',GAChrange,'numYbins',80);
    
for ww = 1:2
    GAChdist.pup.(WHNWH{ww})  = ConditionalHist( log10(GACh.pup(GACh.(WHNWH{ww}))),...
            (GACh.raw(GACh.(WHNWH{ww}))),...
            'numXbins',20,'Xbounds',[-0.25 0.25],'Ybounds',GAChrange,'numYbins',80);
end

%% GACh and Whisking
GAChdist.EMG  = ConditionalHist( log10(GACh.EMG),...
        (GACh.raw),...
        'numXbins',20,'Xbounds',[-1.7 0.9],'Ybounds',GAChrange,'numYbins',80);
    
for oo = 1:2
    for ll = 1:2
        GAChdist.(ONOFF{oo}).(LONGSHORT{ll})  = ConditionalHist( GACh.whtime.(ONOFF{oo})(GACh.(LONGSHORT{ll}).(ONOFF{oo})),...
                (GACh.raw(GACh.(LONGSHORT{ll}).(ONOFF{oo}))),...
                'numXbins',100,'Xbounds',[-5 5],'Ybounds',GAChrange,'numYbins',80);
        
    end
        GAChdist.(ONOFF{oo}).all = ConditionalHist( GACh.whtime.(ONOFF{oo}),...
                (GACh.raw),...
                'numXbins',100,'Xbounds',[-5 5],'Ybounds',GAChrange,'numYbins',80);
            
end
%%
figure
for oo = 1:2
    subplot(4,2,oo)
        imagesc( GAChdist.(ONOFF{oo}).all.Xbins,...
            GAChdist.(ONOFF{oo}).all.Ybins,...
            GAChdist.(ONOFF{oo}).all.pYX')
        hold on; axis xy; box off
        colorbar
        plot([0 0],ylim(gca),'k')
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        xlabel(['t - ',(ONOFF{oo})]);ylabel('GACh')
   crameri bilbao
   
   for ll = 1:2
        subplot(4,2,oo+ll*2)
            imagesc( GAChdist.(ONOFF{oo}).(LONGSHORT{ll}).Xbins,...
                GAChdist.(ONOFF{oo}).(LONGSHORT{ll}).Ybins,...
                GAChdist.(ONOFF{oo}).(LONGSHORT{ll}).pYX')
            hold on; axis xy; box off
            colorbar
            plot([0 0],ylim(gca),'k')
            %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
            xlabel(['t - ',(ONOFF{oo})]);ylabel('GACh')
       crameri bilbao
   end
end

    subplot(4,2,7)
        imagesc( GAChdist.EMG.Xbins,...
            GAChdist.EMG.Ybins,...
            GAChdist.EMG.pYX')
        hold on; axis xy; box off
        colorbar
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        xlabel('EMG');ylabel('GACh')
   crameri bilbao
   NiceSave('GAChWhisk',figfolder,baseName)
%%
cosx = linspace(-pi,pi,100);
cospamp = [0.08 0.8];
figure


subplot(3,2,1)
a = imagesc(pupilphaseGACh.Xbins,pupilphaseGACh.Ybins,pupilphaseGACh.meanZ');
hold on
alpha(a,double(~isnan(pupilphaseGACh.meanZ')))

%imagesc(pupilphaseEMG.Xbins+2*pi,pupilphaseEMG.Ybins,pupilphaseEMG.meanZ')
crameri lapaz
plot([-pi 3*pi],pupilcycle.detectionparms.pupthresh.*[1 1],'w--')
plot(cosx,(cos(cosx)+1).*cospamp(pp)-2,'k')
axis xy
box off
%xlim([-pi 3*pi])
ColorbarWithAxis(GAChrange,'Mean GACh')
%LogScale('c',10)
LogScale('y',10)
xlabel('Pupil Phase');ylabel('Pupil Amplitude')

subplot(3,2,2)
a = imagesc(pupildpGACh.Xbins,pupildpGACh.Ybins,pupildpGACh.meanZ');
hold on
alpha(a,double(~isnan(pupildpGACh.meanZ')))
plot(pupildpGACh.Xbins([1 end]),[0 0],'k--')
crameri lapaz
ColorbarWithAxis(GAChrange,'Mean GACh')
%LogScale('c',10)
axis xy
box off
LogScale('c',10)
LogScale('x',10)



subplot(3,2,3)
        for pp = 1:2
        imagesc( GAChdist.(HILO{pp}).Xbins+2*pi*(pp-1),...
            GAChdist.(HILO{pp}).Ybins,...
            GAChdist.(HILO{pp}).pYX')
        hold on; axis xy; box off
        plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp)-2,'k')
        end   
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        colorbar
        xlim([-pi 3*pi])
        xlabel('Pupil Phase');ylabel('GACh')
        crameri bilbao
        
subplot(3,2,4)

        imagesc( GAChdist.pup.Xbins,...
            GAChdist.pup.Ybins,...
            GAChdist.pup.pYX')
        hold on; axis xy; box off
        colorbar
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        xlabel('Pupil Size');ylabel('GACh')
   crameri bilbao
        NiceSave('GAChPupil',figfolder,baseName)

   %% Figure: whisk/nonwhisk
   figure
for ww = 1:2
 subplot(3,2,ww)

        imagesc( GAChdist.pup.(WHNWH{ww}).Xbins,...
            GAChdist.pup.(WHNWH{ww}).Ybins,...
            GAChdist.pup.(WHNWH{ww}).pYX')
        hold on; axis xy; box off
        colorbar
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        xlabel('Pupil Size');ylabel('GACh')
        title((WHNWH{ww}))
   crameri bilbao
   
subplot(3,2,2+ww)
        for pp = 1:2
        imagesc( GAChdist.(HILO{pp}).(WHNWH{ww}).Xbins+2*pi*(pp-1),...
            GAChdist.(HILO{pp}).(WHNWH{ww}).Ybins,...
            GAChdist.(HILO{pp}).(WHNWH{ww}).pYX')
        hold on; axis xy; box off
        plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp)-2,'k')
        end   
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        colorbar
        xlim([-pi 3*pi])
        xlabel('Pupil Phase');ylabel('GACh')
        crameri bilbao
        
end
   
        NiceSave('GAChPupil_WhNWh',figfolder,baseName)

%% Example Figure
windows(1,:) = [100 250];
windows(2,:) = bz_RandomWindowInIntervals(pupilcycle.timestamps([1 end]),150);
windows(3,:) = bz_RandomWindowInIntervals(pupilcycle.timestamps([1 end]),150);

figure;
for ww = 1:3
subplot(3,1,ww);
%plot(pupildilation.timestamps,pupildilation.data,'k','linewidth',2); hold on;
scatter(pupilcycle.timestamps,pupilcycle.data,3,pupilcycle.phase,'filled')
hold on
plot(GACh.timestamps,bz_NormToRange(GACh.raw,0.6,GAChrange))
scatter(pupilcycle.timestamps,pupilcycle.data,3,pupilcycle.phase,'filled')
%scatter(pupilcycle.timestamps,pupilcycle.amp,3,pupilcycle.phase,'filled')
plot(EMGwhisk.timestamps,EMGwhisk.EMG./max(EMGwhisk.EMG),...
    'color',[0.5 0.5 0.5],'linewidth',0.5);
plot(EMGwhisk.timestamps,5*EMGwhisk.EMGsm./(max(EMGwhisk.EMG)),...
    'color','k','linewidth',0.5);
plot(pupilcycle.timestamps(pupilcycle.states==2),3*ones(sum(pupilcycle.states==2),1),'r.')
plot(pupilcycle.timestamps(pupilcycle.states==1),3*ones(sum(pupilcycle.states==1),1),'k.')
plot(EMGwhisk.timestamps(EMGwhisk.iswhisk),2.8*ones(sum(EMGwhisk.iswhisk),1),'b.')
%plot(highpupildata.timestamps,highpupildata.data+nanmean(pupildilation.data),'r')
colormap(gca,hsv)
ColorbarWithAxis([min(pupilcycle.phase) max(pupilcycle.phase)],['pupil phase'])
%h1 = plot(get(gca,'xlim'),[0 0],'k-');
xlim(windows(ww,:)); ylim([-0.05 3]);
bz_ScaleBar('s')
ylabel('Pupil');
%set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%legend({'diameter','phase','dPdt'},'location','northeast');

end

NiceSave('BehaviorExamples',figfolder,baseName)


%%
%% Measure coupling with a range of lower/uppwer bounds, at a few different ncyc resolutions

lowerbounds = logspace(-2,-1,16);
upperbounds = logspace(-1.5,0,16);
orders = 1:3;

%For each combination of bands

%Filter the pupil

amp2 = interp1(GACh.timestamps,GACh.raw,pupilcycle.timestamps,'nearest');
PAcoupling.lowerbounds = lowerbounds;
PAcoupling.upperbounds = upperbounds;
PAcoupling.orders = orders;
coupling = nan(length(lowerbounds),length(upperbounds),length(orders));
PAcorr = nan(length(lowerbounds),length(upperbounds),length(orders));
for oo = orders
for ll = 1:length(lowerbounds)
    bz_Counter(ll,length(lowerbounds),'lower bound')
    parfor uu = 1:length(upperbounds)
        
        if upperbounds(uu)<=lowerbounds(ll)
            %PAcoupling.coupling(ll,uu,oo)=nan;
            %PAcoupling.corr(ll,uu,oo) = nan;
            continue
        end
        [ pupilcycle_filttest ] = ExtractPupilCycle( basePath,'redetect',true,'saveMat',false,...
            'filterbounds',[lowerbounds(ll) upperbounds(uu)],'filterorder',oo);


    %Normalize power to the median (mean?)
        amp1 = pupilcycle_filttest.amp./nanmean(pupilcycle_filttest.amp);
        
        %Calculate the coupling with EMG
        coupling(ll,uu,oo) = abs(nanmean(amp1.*amp2.*exp(1i.*pupilcycle_filttest.phase)));
        PAcorr(ll,uu,oo) = corr(pupilcycle_filttest.data,amp2,'type','spearman','rows','pairwise');
        

    end
end
end
PAcoupling.coupling = coupling;
PAcoupling.corr = PAcorr;
%%
figure
for oo = orders
    subplot(2,4,oo)
imagesc(log10(lowerbounds),log10(upperbounds),PAcoupling.coupling(:,:,oo)')
hold on
axis xy
%plot(log10(trybounds(1)),log10(trybounds(2)),'r+')
LogScale('xy',10)
%clim([0.25 0.7])
xlabel('Lower Bound (Hz)');ylabel('Upper Bound (Hz)')
colorbar
title(num2str(oo))


    subplot(2,4,oo+4)
imagesc(log10(lowerbounds),log10(upperbounds),PAcoupling.corr(:,:,oo)')
hold on
axis xy
%plot(log10(trybounds(1)),log10(trybounds(2)),'r+')
LogScale('xy',10)
%clim([0.1 0.4])
xlabel('Lower Bound (Hz)');ylabel('Upper Bound (Hz)')
colorbar
title(num2str(oo))
end

NiceSave('FilterCompare',figfolder,baseName)