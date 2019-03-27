function [ eventSpec ] = eventSpec (spec, events, varargin)
%
%
%
%%
twin = round([0.75 0.75].*spec.samplingRate);

% for WaveSpec
spec_temp = NaN(twin(1)+twin(2)+1,size(spec.dataz,2),length(events));
for e = 1:length(events)
    if events(e)-twin(1) > 0 && events(e)+twin(2) < size(spec.data,1)
        spec_temp(:,:,e) = spec.dataz(events(e)-twin(1):events(e)+twin(2),:);
    else
    end
end
eventSpec.spec = (spec_temp-nanmean(spec_temp(twin(1)-80:twin(1)-24,:,:),1))...
    ./nanstd(spec_temp(twin(1)-80:twin(1)-24,:,:),0,1);
eventSpec.spec = nanmean(eventSpec.spec,3);

% for Fractal
spec_temp = NaN(twin(1)+twin(2)+1,size(spec.frac,2),length(events));
for e = 1:length(events)
    if events(e)-twin(1) > 0 && events(e)+twin(2) < size(spec.frac,1)
        spec_temp(:,:,e) = log10(spec.frac(events(e)-twin(1):events(e)+twin(2),:));
    else
    end
end
eventSpec.frac = (spec_temp-nanmean(spec_temp(twin(1)-80:twin(1)-24,:,:),1))...
    ./nanstd(spec_temp(twin(1)-80:twin(1)-24,:,:),0,1);
eventSpec.frac = nanmean(eventSpec.frac,3);

% for PSS
% spec_temp = NaN(twin(1)+twin(2)+1,length(events));
% for e = 1:length(events)
%     if events(e)-twin(1) > 0 && events(e)+twin(2) < size(spec.frac,1)
%         spec_temp(:,e) = spec.Beta(events(e)-twin(1):events(e)+twin(2));
%     else
%     end
% end
% eventSpec.pss = nanmean(spec_temp,2);

% for Osci
spec_temp = NaN(twin(1)+twin(2)+1,size(spec.osci,2),length(events));
for e = 1:length(events)
    if events(e)-twin(1) > 0 && events(e)+twin(2) < size(spec.osci,1)
        spec_temp(:,:,e) = spec.osci(events(e)-twin(1):events(e)+twin(2),:);
    else
    end
end
eventSpec.osci = nanmean(spec_temp,3);

end