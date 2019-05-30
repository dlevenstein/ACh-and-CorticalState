function [frac,osci,validfreq] = WaveIRASA(spec,varargin)
%
% INPUT
% data       power spectrum array in buzcode format (complex-valued)
%
%(options)
%   'logamp'    true if spectrogram is already log-amplitude transformed
%   'freqs'     enter freqs if not using structure input
%
%
% R. Hardstone & W. Munoz - 2019
%% Parse the inputs
% Pending: maxRescaleFactor
%
p = inputParser;
addParameter(p,'logamp',false)
addParameter(p,'freqs',[])
parse(p,varargin{:})
logamp = p.Results.logamp;
freqs = p.Results.freqs;

%%
%fs = spec.samplingRate;
if isstruct(spec)
    numFreqs = spec.nfreqs;
    freqs = spec.freqs;
else
    temp = spec; clear spec;
    spec.data=temp;
    numFreqs = length(freqs);
end

%%
maxRescaleFactor = 2.9; %as per Muthukumaraswamy and Liley, NeuroImage 2018

numberRescalesfreq = maxRescaleFactor*freqs(1);
numberRescalesidx = find(freqs >= numberRescalesfreq);
numberRescales = numberRescalesidx(1)-1;

maxRescaleFactor = freqs(numberRescales+1)/freqs(1);
disp(['MaxRescaleFactor is ' num2str(maxRescaleFactor,2) '. Default IRASA is 2.9']);

%% Resampling
smoothingSamples = 1;
validFreqInds = numberRescales + 1:numFreqs - numberRescales - 1;
% ampData = gpuArray(zeros(size(spec.data)));
% resampledData = gpuArray(zeros(size(spec.data)));
ampData = zeros(size(spec.data));
resampledData = zeros(size(spec.data));

%Extract amplitude from the spectrogram
if logamp
    ampData = spec.data;
else
    for i_freq = 1:numFreqs
        i_freq;
        ampData(:,i_freq) = abs(spec.data(:,i_freq));
        %smoothedData(:,i_freq) = ampData(:,i_freq); %smooth(ampData(:,i_freq),smoothingSamples);
    end
end

%Median smooth the spectrum
for i_freq = validFreqInds
    i_freq;
    inds = [i_freq-numberRescales:i_freq-1 i_freq+1:i_freq+numberRescales];
    resampledData(:,i_freq) = nanmedian(ampData(:,inds),2);
end

% ampData = gather(ampData);
% resampledData = gather(resampledData);

%osci = ampData(:,validFreqInds)-resampledData(:,validFreqInds);
switch logamp
    case false
        osci = log10(ampData(:,validFreqInds))-log10(resampledData(:,validFreqInds));
    case true
        osci = (ampData(:,validFreqInds))-(resampledData(:,validFreqInds));
end
frac = resampledData(:,validFreqInds);
validfreq = freqs(validFreqInds);

%%
% figure
% for i_time = 1:800:size(datlfp.data,1)
%     plot(log10(wavespec.freqs),log10(ampData(i_time,:)));
%     hold on
%     plot(log10(wavespec.freqs),log10(resampledData(i_time,:)));
%     title(i_time/fs);
%     ylim([0 5])
%     drawnow;
%     cla
% end

%%
% figure
% h(1) = subplot(2,1,1);
% plot(log10(wavespec.freqs),log10(mean(ampData)))
% hold on
% plot(log10(wavespec.freqs),log10(mean(resampledData)))
% h(2) = subplot(2,1,2);
% plot(log10(wavespec.freqs(validFreqInds)),log10(mean(ampData(:,validFreqInds)))-log10(mean(resampledData(:,validFreqInds))))
% linkaxes(h,'x')

end
