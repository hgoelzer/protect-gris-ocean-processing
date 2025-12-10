% hindcast for heiko
clear; close all;

% smoothing - choose backwards = 10 and forwards = 9
% for 20 year centred mean
% backwards
sb = 10;
% forwards
sf = 9;

% zero year
yr0 = 2014;

%% calculate past retreat according to parameterisation

% load glaciers
load('/Users/donaldslater/Downloads/glaciers.mat'); % downloaded from github
for ii=1:length(glaciers),
    iceflux(ii) = glaciers(ii).iceflux.final;
end
sectors = [glaciers.sectornum];
% sector indices
for l=1:7,
    ids(l).inds = find(sectors==l);
end

% load K samples
% this comes from create_K_samples.m
load('/Users/donaldslater/Downloads/Ksamples.mat'); % downloaded from github
N = length(k)/length(glaciers);
sectorname = {'SE','SW','CE','CW','NE','NW','NO'};

t = glaciers(1).melt.t;
for ii=1:length(glaciers),
    MA(ii,:) = movmean(glaciers(ii).melt.m,[sb,sf],'omitnan');
    MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
end

MA_N = repmat(MA,N,1);
sampleretreat = k.*MA_N;

% take flux-weighted means
for l=1:7,
    for j=1:N,
        fluxweightedmean(l,j,:) = nansum(iceflux(ids(l).inds)'.*sampleretreat(ids(l).inds+(j-1)*length(glaciers),:))/nansum(iceflux(ids(l).inds));
    end
end

% identify the ensemble members which have quartile retreats at 2100
% and extract retreat timeseries for these members
for l=1:7,
    fwm0 = fluxweightedmean(l,:,end);
    fwm0_sorted = sort(fwm0);
    id25(l) = min(find(fwm0==fwm0_sorted(floor(0.25*length(fwm0_sorted)))));
    id50(l) = min(find(fwm0==fwm0_sorted(floor(0.50*length(fwm0_sorted)))));
    id75(l) = min(find(fwm0==fwm0_sorted(floor(0.75*length(fwm0_sorted)))));
    fw50(l,:) = fluxweightedmean(l,id50(l),:);
    if fw50(l,end)<=0,
        fw25(l,:) = fluxweightedmean(l,id25(l),:);
        fw75(l,:) = fluxweightedmean(l,id75(l),:);
    elseif fw50(l,end)>0,
        fw75(l,:) = fluxweightedmean(l,id25(l),:);
        fw25(l,:) = fluxweightedmean(l,id75(l),:);
    end
end

% assign to output arrays
Rhigh = fw25';
Rmed = fw50';
Rlow = fw75';

% load regional definition
load('/Users/donaldslater/Downloads/ice_ocean_sectors.mat'); % downloaded from github

retreat.regions = regions;
for l=1:7,
    retreat.regions(l).name = sectorname{l};
end

retreat.time = t;
retreat.low = Rlow;
retreat.med = Rmed;
retreat.high = Rhigh;

save hindcast_retreat_290322.mat retreat