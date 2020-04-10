function [bandpower] = GetBandFromOsci(spec,frange,depthrange)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%
freqs = spec.freqs;
normdepth = spec.depth;
    infreq = freqs>=frange(1) & freqs<=frange(2);
    indepth = normdepth>=depthrange(1) & normdepth<=depthrange(2);

if strcmp(frange,'PSS')
    allbandpower = spec.PSS(:,indepth);
    bandpower = median(allbandpower,2);
else
    allbandpower = spec.osci(:,infreq,indepth);
    bandpower = median(allbandpower,[2 3]);
end 
%%


end

