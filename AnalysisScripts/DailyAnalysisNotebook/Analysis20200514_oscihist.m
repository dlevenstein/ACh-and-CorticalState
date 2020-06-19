%%
oscihist.bins = linspace(-1,1.5,100);
%oscihist.median = 
for ff = 1:length(specslope.freqs)
    oscihist.hist(:,ff) = hist(specslope.resid(:,ff),oscihist.bins);
    oscihist.median(ff) = mean(specslope.resid(:,ff),1);
    oscihist.std(ff) = std(specslope.resid(:,ff),[],1);
end

%%

figure
imagesc(oscihist.bins,log10(specslope.freqs),oscihist.hist')
alpha(single(oscihist.hist'>5))
hold on
plot(oscihist.median(:),log10(specslope.freqs),'r.')
plot(oscihist.median(:)+oscihist.std(:),log10(specslope.freqs),'r--')
axis xy
xlabel('Osci');ylabel('f (Hz)')
LogScale('y',10)
title('HippoCampus')