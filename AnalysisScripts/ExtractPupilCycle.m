function [ pupilcycle ] = ExtractPupilCycle( basePath,varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%basePath = '/mnt/proraidDL/Database/WMData/AChPupil/EM1M3/171209_WT_EM1M3';

figfolder = fullfile(basePath,'DetectionFigures');
baseName = bz_BasenameFromBasepath(basePath);
savefilename = fullfile(basePath,[baseName,'.pupilcycle.behavior.mat']);


%%
p = inputParser;
addParameter(p,'redetect',false,@islogical);
addParameter(p,'saveMat',true,@islogical);

parse(p,varargin{:})

redetect = p.Results.redetect;
saveMat = p.Results.saveMat;
%%
if exist(savefilename,'file') & ~redetect
   display('Loading Pupil Cycle')
   load(savefilename)
   return
end
%%
pupildilation = bz_LoadBehavior(basePath,'pupildiameter');

smoothwin_pup = 0.5; %s
pupildilation.data = smooth(pupildilation.data,smoothwin_pup.*pupildilation.samplingRate,'moving');
nantimes = isnan(pupildilation.data);
pupildilation.data = pupildilation.data(~isnan(pupildilation.data));

if length(pupildilation.data) < 1
    warning('Not enough pupil data >)');
    return
end

smoothwin_dpdt = 2; %s
pupildilation.dpdt = diff(smooth(pupildilation.data,smoothwin_dpdt.*pupildilation.samplingRate,'moving')).*pupildilation.samplingRate;
pupildilation.dpdt = smooth(pupildilation.dpdt,smoothwin_dpdt.*pupildilation.samplingRate,'moving');
pupildilation.timestamps = pupildilation.timestamps(~nantimes);

% Filtered Pupil
lowfilter = [0.02 0.2]; %new: EMG coupled. order 1


pupil4filter = pupildilation;
filteredpupil = bz_Filter(pupil4filter,'passband',lowfilter,'filter' ,'fir1','order',1);
%highpupildata = bz_Filter(pupil4filter,'passband',highfilter,'filter' ,'fir1');
filteredpupil.pupthresh = -0.8;
filteredpupil.highpup = log10(filteredpupil.amp)>filteredpupil.pupthresh; 
filteredpupil.lowpup = log10(filteredpupil.amp)<=filteredpupil.pupthresh; 


%%
pupilcycle.data = pupildilation.data;
pupilcycle.timestamps = filteredpupil.timestamps;
pupilcycle.dpdt = pupildilation.dpdt;
pupilcycle.phase = filteredpupil.phase;
pupilcycle.amp = filteredpupil.amp;
pupilcycle.states = zeros(size(filteredpupil.timestamps));
pupilcycle.statenames = {'lowpup','highpup'};
pupilcycle.states(filteredpupil.highpup) = 2;
pupilcycle.states(filteredpupil.lowpup) = 1;


pupilcycle.ints = bz_IDXtoINT(pupilcycle);
%%
pupilcycle.dur.highpup = diff(pupilcycle.ints.highpupstate,[],2);
pupilcycle.dur.lowpup = diff(pupilcycle.ints.lowpupstate,[],2);

%%
durhist.bins = linspace(0,2,25);
durhist.pre.low = hist(log10(pupilcycle.dur.lowpup),durhist.bins);
durhist.pre.high = hist(log10(pupilcycle.dur.highpup),durhist.bins);

%% Merge short windows

highpupdurthresh = 5; %s

shortHIints = pupilcycle.dur.highpup<=highpupdurthresh;

switchIDX = InIntervals(pupilcycle.timestamps,pupilcycle.ints.highpupstate(shortHIints,:));
pupilcycle.states(switchIDX) = 1;   
pupilcycle.ints = bz_IDXtoINT(pupilcycle);

pupilcycle.dur.highpup = diff(pupilcycle.ints.highpupstate,[],2);
pupilcycle.dur.lowpup = diff(pupilcycle.ints.lowpupstate,[],2);


lowpupdurthresh = 5; %s

shortLOints = pupilcycle.dur.lowpup<=lowpupdurthresh;

switchIDX = InIntervals(pupilcycle.timestamps,pupilcycle.ints.lowpupstate(shortLOints,:));
pupilcycle.states(switchIDX) = 2;   
pupilcycle.ints = bz_IDXtoINT(pupilcycle);

pupilcycle.dur.highpup = diff(pupilcycle.ints.highpupstate,[],2);
pupilcycle.dur.lowpup = diff(pupilcycle.ints.lowpupstate,[],2);


%% duration histogram after

durhist.post.low = hist(log10(pupilcycle.dur.lowpup),durhist.bins);
durhist.post.high = hist(log10(pupilcycle.dur.highpup),durhist.bins);
%% Buzcode output

pupilcycle.detectionparms.smoothwin_pupil = smoothwin_pup;
pupilcycle.detectionparms.smoothwin_dpdt = smoothwin_dpdt;
pupilcycle.detectionparms.filterparms = filteredpupil.filterparms;
pupilcycle.detectionparms.lowpupdurthresh = lowpupdurthresh;
pupilcycle.detectionparms.highpupdurthresh = highpupdurthresh;
pupilcycle.detectionparms.pupthresh = filteredpupil.pupthresh;

if saveMat
    display(['Saving ',baseName,'.pupilcycle.behavior.mat for next time'])
    save(savefilename,'pupilcycle')
end
%% Example Figure
if saveMat
%windows(1,:) = [100 400];
windows(2,:) = bz_RandomWindowInIntervals(pupildilation.timestamps([1 end]),300);
windows(3,:) = bz_RandomWindowInIntervals(pupildilation.timestamps([1 end]),300);

figure;
subplot(3,3,1)
hist(log10(pupilcycle.amp))
hold on
%axis tight
box off 
LogScale('x',10)
xlabel('Pupil Amplitude')
plot(pupilcycle.detectionparms.pupthresh.*[1 1],get(gca,'ylim'),'r')
subplot(3,3,2)
bar(durhist.bins,durhist.post.low,'facecolor',[0.5 0.5 0.5])
hold on
bar(durhist.bins,durhist.pre.low,'facealpha',0,'edgecolor','k')
xlabel('Low Pupil Durations (s)')
axis tight
box off
LogScale('x',10)
subplot(3,3,3)
bar(durhist.bins,durhist.post.high,'facecolor',[1 0.5 0.5])
hold on
bar(durhist.bins,durhist.pre.high,'facealpha',0,'edgecolor','r')
xlabel('High Pupil Durations (s)')
axis tight
box off
LogScale('x',10)

for ww = 2:3
subplot(3,1,ww);
scatter(pupilcycle.timestamps,pupilcycle.data,3,pupilcycle.phase,'filled')
hold on
plot(pupilcycle.timestamps,pupilcycle.amp.*2)
% plot(EMGwhisk.timestamps,EMGwhisk.EMG./max(EMGwhisk.EMG),...
%     'color',[0.5 0.5 0.5],'linewidth',0.5);
% plot(EMGwhisk.timestamps,5*EMGwhisk.EMGsm./(max(EMGwhisk.EMG)),...
%     'color','k','linewidth',0.5);
%plot(pupilcycle.timestamps(pupilcycle.highpup),3*ones(sum(pupilcycle.highpup),1),'r.')
%plot(pupilcycle.timestamps(~pupilcycle.highpup),3*ones(sum(~pupilcycle.highpup),1),'k.')

plot(pupilcycle.ints.highpupstate',3*ones(size(pupilcycle.ints.highpupstate')),'r')
plot(pupilcycle.ints.lowpupstate',3*ones(size(pupilcycle.ints.lowpupstate')),'k')

%plot(EMGwhisk.timestamps(EMGwhisk.iswhisk),2.8*ones(sum(EMGwhisk.iswhisk),1),'b.')
%plot(highpupildata.timestamps,highpupildata.data+nanmean(pupildilation.data),'r')
colormap(gca,hsv)
%ColorbarWithAxis([min(filteredpupil.phase) max(filteredpupil.phase)],['pupil phase'])
%h1 = plot(get(gca,'xlim'),[0 0],'k-');
xlim(windows(ww,:)); ylim([-0.05 3]);
bz_ScaleBar('s')
ylabel('Pupil');
%set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%legend({'diameter','phase','dPdt'},'location','northeast');

end
NiceSave('PupilCycle',figfolder,baseName)
end


end
