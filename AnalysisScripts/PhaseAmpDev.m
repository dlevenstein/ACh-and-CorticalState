basePath = '/mnt/proraidDL/Database/WMProbeData/180213_WT_M1M3_LFP_Layers_Pupil_EMG_Pole/180213_WT_M1M3_LFP_Layers_Pupil_EMG_180213_113045';

baseName = bz_BasenameFromBasepath(basePath);


%% PUPIL %%
pupildilation = bz_LoadBehavior(basePath,'pupildiameter');

smoothwin =2;%s
pupildilation.dpdt = diff(smooth(pupildilation.data,smoothwin.*pupildilation.samplingRate,'moving')).*pupildilation.samplingRate;
pupildilation.dpdt = smooth(pupildilation.dpdt,smoothwin.*pupildilation.samplingRate,'moving');

nantimes = isnan(pupildilation.data);
pupildilation.interpdata = interp1(pupildilation.timestamps(~nantimes),...
    pupildilation.data(~nantimes),pupildilation.timestamps);


lowfilter = [0.01 0.1];
highfilter = [0.3 0.8];

pupil4filter = pupildilation;
pupil4filter.data = pupildilation.interpdata(~isnan(pupildilation.interpdata));
%pupil4filter.t
pupil4filter.timestamps = pupil4filter.timestamps(~isnan(pupildilation.interpdata));
lowpupildata = bz_Filter(pupil4filter,'passband',lowfilter,'filter' ,'fir1','order',3);

%% WHISK %%

EMGwhisk = bz_LoadStates( basePath,'EMGwhisk');


%%
EMGwhisk.pupiltime = interp1(EMGwhisk.t,EMGwhisk.EMGsm,lowpupildata.timestamps);


%%
figure
subplot(2,1,1)
    plot(pupil4filter.timestamps,pupil4filter.data,'k')
    hold on
    plot(lowpupildata.timestamps,lowpupildata.data,'b')
subplot(2,1,2)
    plot(lowpupildata.timestamps,log10(EMGwhisk.pupiltime),'k')

%% Phase Amp coupling

PhaseAmpCouplingByAmp( lowpupildata.phase,log10(lowpupildata.amp),...
    (EMGwhisk.pupiltime),'showFig',true,'AmpBounds',[-1.75 -0.5],...
    'shufflesig',true,'AmpZNorm',false,'numAmpbins',15);
