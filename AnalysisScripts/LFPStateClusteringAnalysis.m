%LFPStateAnalysis

basePath= '/mnt/proraidDL/Database/WMProbeData/170421_Layers_LFP_Pupil_EMG_Emx1M1M3/170421_Layers_LFP_Pupil_EMG_170421_180427';
usechans = 'all';

LFP = bz_GetLFP(usechans,'basepath',basePath);

%% Parms

specparms.type = 'wavelet';
specparms.frange = [1 128];
specparms.nfreqs = 50;
specparms.ncyc = 4;
specparms.space = 'log';
specparms.specnorm = 'log';


%specparms.type = 'FFT';
%specparms.winsize = 3; %s
%specparms.noverlap = 2;
%%
clear spec
numchannels = length(LFP.channels);
for cc = 1:numchannels
%for cc = 1:2
    cc
    %spectrogram
    %parameters: dt = 1;
    %window = 4;
    
    switch specparms.type
        case 'wavelet'
            %Calcualte the Wavelet Transform
            [freqs,~,spec(:,:,cc)] = WaveSpec(single(LFP.data(:,cc)),...
                specparms.frange,specparms.nfreqs,specparms.ncyc,...
                1/LFP.samplingRate,specparms.space);
            spectimestamps = LFP.timestamps; %Wavelet timestamp are same as LFP        

        case 'FFT'
            %Calculate the frequences to use
            switch specparms.space
                case 'log'
                    freqs = logspace(log10(specparms.frange(1)),...
                        log10(specparms.frange(2)),specparms.nfreqs);
                case 'lin'
                    freqs = linspace(specparms.frange(1),...
                        specparms.frange(2),specparms.nfreqs);  
            end

            %Calculate the FFT spectrogram parameters - covert from s to sf
            winsize = specparms.winsize*LFP.samplingRate;
            noverlap = specparms.noverlap*LFP.samplingRate; %Might want to calaulte overlap based on pupil...?
            %Calculate the FFT spectrogram
            [spec(:,:,cc),~,spectimestamps] = spectrogram(single(LFP.data(:,cc)),...
                winsize,noverlap,freqs,LFP.samplingRate);
            spectimestamps = spectimestamps'+LFP.timestamps(1); %account for any time offset

    end
end
spec = permute(spec,[2,1,3]);
%Normalize the spectrogram
switch specparms.specnorm 
    case 'log'
        spec = log10(abs(spec));
    case 'mean'
        spec = (abs(spec));
        spec = bsxfun(@(X,Y) X./Y,spec,mean(spec,1));
    case 'logmean'
        spec = (abs(spec));
        spec = bsxfun(@(X,Y) X./Y,spec,mean(spec,1));
        spec = log10(spec);
end
    


%%
%put all the spectrograms into a big time vector
spec = reshape(spec,[],specparms.nfreqs.*numchannels);
%%
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(spec);
%%
figure
plot(cumsum(EXPLAINED),'o-')
%%
figure
imagesc(spec)
%%
%use tsne....
%color by pupil d, d/dt
mappedX = tsne(spec);

%%
figure
plot(mappedX(:,1),mappedX(:,2),'.')