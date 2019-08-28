function [  ] = MakeFigures(basePath,figfolder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% Load Header
%Initiate Paths
%reporoot = '/home/dlevenstein/ProjectRepos/ACh-and-CorticalState/';
%reporoot = '/Users/dlevenstein/Project Repos/ACh-and-CorticalState/';
%basePath = '/mnt/proraidDL/Database/WMData/AChPupil/171209_WT_EM1M3/';
%basePath = '/mnt/proraidDL/Database/WMData/AChPupil/180706_WT_EM1M3/';
%basePath = pwd;
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/MakeFigures'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath,'noPrompts',true);

%%
%Restricting SPONT UP/DOWNs
load(fullfile(basePath,[baseName,'.MergePoints.events.mat']),'MergePoints');
sidx = find(startsWith(MergePoints.foldernames,"Spont"));
sponttimes = [MergePoints.timestamps(sidx(1),1) MergePoints.timestamps(sidx(end),2)];

%% Loading behavior...
% Pupil diameter
pupildilation = bz_LoadBehavior(basePath,'pupildiameter');
pupildilation.data_raw = pupildilation.data;
smoothwin = 0.5; %s
pupildilation.data = smooth(pupildilation.data,smoothwin.*pupildilation.samplingRate,'moving');
nantimes = isnan(pupildilation.data);
pupildilation.data = pupildilation.data(~isnan(pupildilation.data));
pupildilation.data_raw = pupildilation.data_raw(~isnan(pupildilation.data));

if length(pupildilation.data) < 1
    warning('Not enough pupil data >)');
    return
end

smoothwin = 2; %s
pupildilation.dpdt = diff(smooth(pupildilation.data,smoothwin.*pupildilation.samplingRate,'moving')).*pupildilation.samplingRate;
pupildilation.dpdt = smooth(pupildilation.dpdt,smoothwin.*pupildilation.samplingRate,'moving');
pupildilation.timestamps = pupildilation.timestamps(~nantimes);

% Filtered Pupil
lowfilter = [0.01 0.1]; %old. order 3
lowfilter = [0.02 0.2]; %new: EMG coupled. order 1
%highfilter = [0.3 0.8];

pupil4filter = pupildilation;
pupilcycle = bz_Filter(pupil4filter,'passband',lowfilter,'filter' ,'fir1','order',1);
%highpupildata = bz_Filter(pupil4filter,'passband',highfilter,'filter' ,'fir1');
pupilcycle.pupthresh = -0.8;
pupilcycle.highpup = log10(pupilcycle.amp)>pupilcycle.pupthresh; 

% EMG
EMGwhisk = bz_LoadBehavior(basePath,'EMGwhisk');

spontidx = find(EMGwhisk.ints.Wh(:,2) < sponttimes(2));
EMGwhisk.ints.Wh = EMGwhisk.ints.Wh(spontidx,:);

spontidx = find(EMGwhisk.ints.NWh(:,2) < sponttimes(2));
EMGwhisk.ints.NWh = EMGwhisk.ints.NWh(spontidx,:);

spontidx = find(EMGwhisk.timestamps < sponttimes(2));
EMGwhisk.timestamps = EMGwhisk.timestamps(spontidx);
EMGwhisk.EMGenvelope = EMGwhisk.EMGenvelope(spontidx);
EMGwhisk.EMG = EMGwhisk.EMG(spontidx);
EMGwhisk.EMGsm = EMGwhisk.EMGsm(spontidx);
%% Get the depth info

depthinfo = rescaleCx(basePath);

%Get example channels that are closest to midpoints between bondaries (for
%LFP)
exdepths = depthinfo.boundaries(1:end-1) + 0.5*diff(depthinfo.boundaries);
for ee = 1:length(exdepths)
    distfromdepth = abs(exdepths(ee)-depthinfo.ndepth);
    exChanIDX(ee) = find(distfromdepth == min(distfromdepth),1);
end

exChan = depthinfo.channels(exChanIDX);

%Get depth info for all cortical channels (for PSS)
inCTX = find(~isnan(depthinfo.ndepth));
CTXchans = depthinfo.channels(inCTX);
CTXdepth = -depthinfo.ndepth(inCTX);



%% Load only the cortex channels in the spontaneous time window
downsamplefactor = 5;
lfp = bz_GetLFP(exChan,...
    'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor,...
    'intervals',sponttimes);
lfp.chanlayers = depthinfo.layer(exChanIDX);
lfp.chandepths = -depthinfo.ndepth(exChanIDX);
%% Calculate Spectrogram on all channels
clear spec
dt = 0.1;
spec.winsize = 1;
spec.frange = [0.5 256]; %Frequency lower than can be assessed for window because IRASA... but maybe this is bad for IRASA too
spec.nfreqs = 100;

noverlap = spec.winsize-dt;
spec.freqs = logspace(log10(spec.frange(1)),log10(spec.frange(2)),spec.nfreqs);
winsize_sf = round(spec.winsize .*lfp.samplingRate);
noverlap_sf = round(noverlap.*lfp.samplingRate);

spec.channels = lfp.channels;
for cc =1:length(spec.channels)
    cc
    
    load(fullfile(basePath,'WaveSpec_Downsampled',[baseName,'.',num2str(spec.channels(cc)),'.WaveSpec.lfp.mat']));
    inspont = InIntervals(wavespec.timestamps,sponttimes);
    spec.data(:,:,cc) = log10(abs(wavespec.data(inspont,wavespec.freqs>=1)));
    spec.data(:,:,cc) = NormToInt(spec.data(:,:,cc),'median');
    spec.timestamps = wavespec.timestamps(inspont,:);
    spec.freqs = wavespec.freqs(wavespec.freqs>=1); 

    clear wavespec
    
end
spec.chanlayers = lfp.chanlayers;




%% Load PSS and aligning behavior
%if strcmp(baseName,'171209_WT_EM1M3')
%    load(fullfile(basePath,[baseName,'.PowerSpectrumSlope.mat']));
%else
    load(fullfile(basePath,[baseName,'.PowerSpectrumSlope.lfp.mat']));
%end

%%

clear PSS
[~,~,PSS.CTXchans] = intersect(CTXchans,specslope.channels,'stable');
PSS.data = specslope.data(:,PSS.CTXchans);
PSS.timestamps = specslope.timestamps;
PSS.samplingRate = specslope.samplingRate;
PSS.winsize = specslope.detectionparms.winsize;
PSS.depth = CTXdepth;
PSS.chan = specslope.channels(PSS.CTXchans);
PSS.osci = specslope.resid(:,:,PSS.CTXchans);
PSS.freqs = specslope.freqs;
PSS.chanlayers = depthinfo.layer(inCTX);

LAYERS = depthinfo.lnames;

for ll = 1:length(LAYERS)
    PSS.Lchans.(LAYERS{ll}) = strcmp(LAYERS{ll},PSS.chanlayers);
end

% Interpolate PSS to normalized depth
PSS.interpdepth = linspace(-1,0,100);
PSS.depthinterp = interp1(PSS.depth',PSS.data',PSS.interpdepth')';
%% Get rid of recording start/stop artifact and restrict to spontaneous behavior



inspont = InIntervals(PSS.timestamps,sponttimes);
PSS.timestamps = PSS.timestamps(inspont);
PSS.data = PSS.data(inspont,:);

maxtimejump = 1; %s
pupilcycle.amp = NanPadJumps( pupilcycle.timestamps,pupilcycle.amp,maxtimejump );
pupilcycle.phase = NanPadJumps( pupilcycle.timestamps,pupilcycle.phase,maxtimejump );
pupildilation.data = NanPadJumps( pupildilation.timestamps,pupildilation.data,maxtimejump );
EMGwhisk.EMGsm = NanPadJumps( EMGwhisk.timestamps,EMGwhisk.EMGsm,maxtimejump );

PSS.Wh = InIntervals(PSS.timestamps,EMGwhisk.ints.Wh);
EMGwhisk.ints.ExpWh = bsxfun(@plus,EMGwhisk.ints.Wh,[-0.5 0.5]*PSS.winsize);
PSS.NWh = ~InIntervals(PSS.timestamps,EMGwhisk.ints.ExpWh);
PSS.hipup = interp1(pupilcycle.timestamps,single(pupilcycle.highpup),PSS.timestamps,'nearest')==1;
PSS.lopup = interp1(pupilcycle.timestamps,single(~pupilcycle.highpup),PSS.timestamps,'nearest')==1;


%% Pupil phase/amp and duration at whisks

EMGwhisk.whisks.pupphase = interp1(pupilcycle.timestamps,pupilcycle.phase,EMGwhisk.ints.Wh(:,1),'nearest');
EMGwhisk.whisks.pupamp = interp1(pupilcycle.timestamps,pupilcycle.amp,EMGwhisk.ints.Wh(:,1),'nearest');
EMGwhisk.whisks.pup = interp1(pupildilation.timestamps,pupildilation.data,EMGwhisk.ints.Wh(:,1),'nearest');
EMGwhisk.whisks.dur = diff(EMGwhisk.ints.Wh,[],2);
EMGwhisk.whisks.hipup = log10(EMGwhisk.whisks.pupamp)>pupilcycle.pupthresh;
EMGwhisk.whisks.lopup = ~EMGwhisk.whisks.hipup;

%%
LAYERS = {'L1','L23','L4','L5a','L5b6','L6'};
depthinfo.boundaries = [0 0.1 0.35 0.5 1];
%% Figure: long time scale (~100s)
winsize = [10:10:70];
winsize = reshape([winsize' winsize']',[],1);
for ww =1:length(winsize)
    exwins(ww,:) = bz_RandomWindowInIntervals(pupildilation.timestamps([1 end]),winsize(ww));
end

whiskwin = [-4 4];
%Long Whisk on dilation
dilationwhisk = EMGwhisk.whisks.hipup & EMGwhisk.whisks.pupphase<-0.25 & ...
    EMGwhisk.whisks.pupphase>-1.5 & EMGwhisk.whisks.dur>2 & EMGwhisk.whisks.pup>1.5;
dilationwhisk = EMGwhisk.ints.Wh(dilationwhisk,1);
exwins(ww+1,:) = randsample(dilationwhisk,1)+ whiskwin;
exwins(ww+2,:) = randsample(dilationwhisk,1)+ whiskwin;
%Whisk on constriction
constrictionwhisk = EMGwhisk.whisks.hipup & EMGwhisk.whisks.pupphase>0.25;
constrictionwhisk = EMGwhisk.ints.Wh(constrictionwhisk,1);
exwins(ww+3,:) = randsample(constrictionwhisk,1)+ whiskwin;
exwins(ww+4,:) = randsample(constrictionwhisk,1)+ whiskwin;

%Whisk on tonic
tonicwhisk = EMGwhisk.whisks.lopup & EMGwhisk.whisks.dur<1;
tonicwhisk = EMGwhisk.ints.Wh(tonicwhisk,1);
exwins(ww+5,:) = randsample(tonicwhisk,1)+ whiskwin;
exwins(ww+6,:) = randsample(tonicwhisk,1)+ whiskwin;
%%
for ww = 1:length(exwins)
    timewin = exwins(ww,:);

lfpinwin = InIntervals(lfp.timestamps,timewin);
%pssinwin = InIntervals(PSS.timestamps,timewin);
exchan_spec = strcmp(spec.chanlayers,'L5b6');
exchan_PSS = find(spec.channels(exchan_spec)==PSS.chan);

speclim = [0.65 1.35]; %Med norm
PSSrange = [-2.4 -1.2];
PSSrange = [-1.5 -0.25]; %new pss
%PSSrange = [-2.6 -0.8];

if mod(ww,2)==1
figure
end
subplot(8,2,[1 3]+mod(ww,2))
plot(pupildilation.timestamps,pupildilation.data,'r','Linewidth',2)
hold on
plot(EMGwhisk.timestamps,EMGwhisk.EMG./50+0.25,'color',0.5.*[1 1 1])
plot(EMGwhisk.timestamps,EMGwhisk.EMGsm./9+0.25,'k')
ylim([0.2 2.4])
if ww>14
   plot(mean(timewin).*[1 1], ylim(gca),'k')
end
ylabel('Pupil/EMG')
%colorbar
box off
xlim(timewin)
    bz_ScaleBar('s')

subplot(8,2,[13 15]+mod(ww,2))
imagesc(spec.timestamps,log10(spec.freqs),spec.data(:,:,exchan_spec)')
hold on
ylim([-1 2.5])
plot(PSS.timestamps,bz_NormToRange(PSS.data(:,exchan_PSS),[0 2.5]),'k','LineWidth',2)
plot(lfp.timestamps(lfpinwin),bz_NormToRange(single(lfp.data(lfpinwin,exchan_spec)),[-1 -0.1]),'k','LineWidth',0.1)
axis xy
ylabel('L5 LFP - f (Hz)')
if ww>14
   plot(mean(timewin).*[1 1], ylim(gca),'k')
end
LogScale('y',10)
        %ColorbarWithAxis(speclim,'Power (med^-^1)')
        clim(speclim)
       box off
       crameri vik
xlim(timewin)
bz_ScaleBar('s')



depthscalefact=0.5e5;
if ww>14
    lfpscalefact = 1.5;
else
    lfpscalefact = 1;
end
subplot(4,2,[3 5]+mod(ww,2))
imagesc(PSS.timestamps,PSS.interpdepth*depthscalefact,PSS.depthinterp')
axis xy
hold on
bz_MultiLFPPlot(lfp,'timewin',timewin,'LFPlabels',lfp.chanlayers,...
    'LFPmidpoints',lfp.chandepths*depthscalefact,'lfpcolor','k',...
    'lfpwidth',0.1,'scaleLFP',lfpscalefact);
%ColorbarWithAxis(PSSrange,'PSS')
hold on
plot(timewin,-depthinfo.boundaries'*[1 1].*depthscalefact,'w','linewidth',1)
clim(PSSrange)
ylim([-1.05 -0.025]*depthscalefact)
    %crameri bamako
    if ww>14
   plot(mean(timewin).*[1 1], ylim(gca),'k')
    end
xlim(timewin)
    bz_ScaleBar('s')
    
if mod(ww,2)==0
NiceSave(['ExampleWin',num2str(ww)],figfolder,baseName)
end
end

%% PSS example figure
if strcmp(baseName,'171209_WT_EM1M3')
    PSS.spec = specslope.specgram(:,:,PSS.CTXchans);
    %PSS.frac = PSpecSlope.Shortwin.FRAC(:,:,PSS.CTXchans);
%exwins = [600 775];
%exwins = [lfp.timestamps(1) lfp.timestamps(end)];
%exwins = [600 666];
%exwins = [100 250];
exwins = [1175 1300;1242 1270;1260 1266;1244 1250];
for ww = 1:length(exwins)
%for ww = 2
    timewin = exwins(ww,:);

lfpinwin = InIntervals(lfp.timestamps,timewin);
%pssinwin = InIntervals(PSS.timestamps,timewin);
exchan_spec = strcmp(spec.chanlayers,'L5b6');
exchan_PSS = find(spec.channels(exchan_spec)==PSS.chan);

speclim = [0.65 1.35]; %Med norm
PSSrange = [-2.4 -1.2];
PSSrange = [-1.5 -0.25]; %new pss
%PSSrange = [-2.6 -0.8];

figure
subplot(8,1,1:2)
%plot(pupildilation.timestamps,pupildilation.data_raw,'k--','Linewidth',2)
hold on
plot(pupildilation.timestamps,pupildilation.data,'r','Linewidth',2)
hold on
plot(EMGwhisk.timestamps,EMGwhisk.EMG./50+0.4,'color',0.5.*[1 1 1])
plot(EMGwhisk.timestamps,EMGwhisk.EMGsm./9+0.4,'k')
ylim([0.35 2.5])
if ww>7
   plot(mean(timewin).*[1 1], ylim(gca),'k')
end

xlim(timewin)
%colorbar
box off
bz_ScaleBar('s')

subplot(8,1,3:4)
imagesc(spec.timestamps,log10(spec.freqs),spec.data(:,:,exchan_spec)')
hold on
ylim([-1 2.5])
%plot(PSS.timestamps,bz_NormToRange(PSS.data(:,exchan_PSS),[0 2.5]),'k','LineWidth',2)
plot(lfp.timestamps(lfpinwin),bz_NormToRange(single(lfp.data(lfpinwin,exchan_spec)),[-1 -0.1]),'k','LineWidth',0.1)
axis xy
if ww>7
   plot(mean(timewin).*[1 1], ylim(gca),'k')
end
LogScale('y',10)
        %ColorbarWithAxis(speclim,'Power (med^-^1)')
        clim(speclim);
       box off
       crameri vik

yyaxis right
plot(PSS.timestamps,PSS.data(:,exchan_PSS),'k','LineWidth',2)
ylim([-3.5 -0.85])
xlim(timewin)
    bz_ScaleBar('s')

depthscalefact=0.5e5;
if ww>2
    lfpscalefact = 1.5;
else
    lfpscalefact = 1;
end


subplot(2,1,2)
imagesc(PSS.timestamps,PSS.interpdepth*depthscalefact,PSS.depthinterp')
axis xy
hold on
bz_MultiLFPPlot(lfp,'timewin',timewin,'LFPlabels',lfp.chanlayers,...
    'LFPmidpoints',lfp.chandepths*depthscalefact,'lfpcolor','k',...
    'lfpwidth',0.1,'scaleLFP',lfpscalefact);
%ColorbarWithAxis(PSSrange,'PSS')
clim(PSSrange)
ylim([-1.05 -0.025]*depthscalefact)

    %crameri bamako
if ww>7
   plot(mean(timewin).*[1 1], ylim(gca),'k')
end    
xlim(timewin)
    bz_ScaleBar('s')   

NiceSave(['IllustrateExample',num2str(ww)],figfolder,baseName)

if ww == 2
figure
    subplot(8,1,7:8)
        imagesc(spec.timestamps,log10(spec.freqs),spec.data(:,:,exchan_spec)')
        hold on
        ylim([-1 2.5])
        %plot(PSS.timestamps,bz_NormToRange(PSS.data(:,exchan_PSS),[0 2.5]),'k','LineWidth',2)
        plot(lfp.timestamps(lfpinwin),bz_NormToRange(single(lfp.data(lfpinwin,exchan_spec)),[-1 -0.1]),...
            'k','LineWidth',0.1)
        axis xy
        if ww>7
           plot(mean(timewin).*[1 1], ylim(gca),'k')
        end
        LogScale('y',10)
                %ColorbarWithAxis(speclim,'Power (med^-^1)')
                clim(speclim);
               box off
               crameri vik
        xlim(timewin)
        bz_ScaleBar('s')
        yyaxis right
        plot(PSS.timestamps,PSS.data(:,exchan_PSS),'k','LineWidth',2)
        ylim([-3.5 -0.85])
    
        
        times = (0.5+timewin(1)):2:timewin(2);
        usefreqs = PSS.freqs>=1 & PSS.freqs<=100;
        fracfreqs = PSS.freqs>=2 & PSS.freqs<=100;
       
    subplot(4,1,1)
    for tt = 1:length(times)
        hold on
        [~,timepoint] = min(abs(PSS.timestamps-(times(tt))));
        
        %linefit.validfreq = PSS.freqs; linefit.frac = PSS.frac(:,timepoint,exchan_PSS)';
        %[linefit] = WaveIRASA_plawfit( linefit, [2.5 100] );
        
        plot(log10(PSS.freqs(usefreqs))+times(tt)-0.5,(PSS.spec(timepoint,usefreqs,exchan_PSS)),'k','linewidth',1)
        plot(log10(PSS.freqs(usefreqs))+times(tt)-0.5,specslope.data(timepoint,exchan_PSS)*log10(PSS.freqs(usefreqs))+specslope.intercept(timepoint,exchan_PSS),'r','linewidth',2)
    
    end
    %clim(PSSrange)
        axis tight
        box off
        xlabel('f (Hz)');ylabel('Power (dB)')
        %LogScale('x',10)
        xlim(timewin)
        %xlim(timewin)
            NiceSave(['IllustrateExamplePSS',num2str(ww)],figfolder,baseName)   

end  
end
end

end
