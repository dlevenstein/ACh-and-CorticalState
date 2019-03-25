%% fitting power-law function to scale-free power-spectrum
%   amri_sig_plawfit()
%
% Usage
%   spec = amri_sig_plawfit(spec, frange)
%
% Inputs
%   spec   - spec.freq: frequency points
%          - spec.frac: the scale-free power spectrum.
%   frange - given frequency range for power-law fitting
%
% Outputs
%   spec  - spectrum
%           .Freq = a vector of frequency points with given frequency range
%           .Plaw = power-law spectrum
%           .Beta = the power-law exponent within the given frequency range
%           .Cons = the power intersect of power-law line in log-log scale.
%           .rsq
%
% Version
%   0.01
%
% Reference
%   -Wen H. and Liu Z. (2015) Separating Fractal and Oscillatory Components in
%    the Power Spectrum of Neurophysiological Signals
%
%% History
% 0.01 - HGWEN - 12/20/2013 - original file
%

%%
function [spec] = WaveIRASA_plawfit( spec, frange )

% define frequency range
ff = spec.validfreq >= frange(1) & spec.validfreq <= frange(2);
freq = spec.validfreq(ff);
frac = spec.frac(:,ff)';

% convert to log-log scale
logfreq = log10(freq);
y1 = log10(frac);

% resample frac in equal space
x2 = linspace(min(logfreq),max(logfreq),length(logfreq)); x2 = x2(:);
y2 = interp1(logfreq,y1,x2);

% fitting power-law function
Nt = size(y2,2);
Nc = size(y2,3);
beta = zeros(Nt,Nc);
rsq = zeros(Nt,Nc);
cons = zeros(Nt,Nc);
plaw = zeros(size(frac));

for j = 1 : Nc
    for i = 1 : Nt
        % ordinary least square
        p = polyfit(x2,y2(:,i,j),1);
        y = squeeze(y2(:,i,j));
        
        %Calculate the residuals
        yfit =  p(1) * x2 + p(2);
        yresid = y - yfit;
        
        %Calculate the rsquared value
        SSresid = sum(yresid.^2);
        SStotal = (length(y)-1) * var(y);
        rsq(i,j) = 1 - SSresid/SStotal;
        beta(i,j) = -p(1);
        cons(i,j) = p(2);
        powlaw = 10.^(polyval(p,logfreq));
        plaw(:,i,j) = powlaw(:);
        
    end
end

% outputs
spec.Beta = beta.*-1;
spec.Cons = cons;
spec.Plaw = plaw;
spec.fitFreq = freq;
spec.rsq = rsq;
end