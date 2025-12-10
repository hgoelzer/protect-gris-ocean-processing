% script to do final ISMIP6 retreat projections
clear; close all;

% smoothing - choose backwards = 10 and forwards = 9
% for 20 year centred mean
% backwards
sb = 10;
% forwards
sf = 9;

% zero year
yr0 = 2014;

plot_flg = 1;

%% calculate retreat

% load glaciers
load ../glaciers/glaciers_SDBN1.mat
for ii=1:length(glaciers),
    iceflux(ii) = glaciers(ii).iceflux.final;
end
sectors = [glaciers.sectornum];
% sector indices
for l=1:7,
    ids(l).inds = find(sectors==l);
end
% make the 8th sector the largest glacier by ice flux in each region
ids(8).inds = [];
for l=1:7,
    ids(8).inds = [ids(8).inds,find(iceflux==max(iceflux(find(sectors==l))))];
end

% load K samples
% this comes from create_K_samples.m
load Ksamples.mat
N = length(k)/length(glaciers);

% loop over RCPs and model scenarios
% mm = 1 for CESM2-WACCM ssp585
% mm = 2 for CESM2-WACCM ssp126

sectorname = {'SE','SW','CE','CW','NE','NW','NO'};
for mm = [1, 2]
%for mm = [1]

    % make array of projected sample retreat
    if mm == 1,
        t = glaciers(1).CESM2WACCM.ssp585.tmelt;
        for ii=1:length(glaciers),
            MA(ii,:) = movmean(glaciers(ii).CESM2WACCM.ssp585.melt,[sb,sf]);
            MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
        end
    elseif mm == 2,
        t = glaciers(1).CESM2WACCM.ssp126.tmelt;
        for ii=1:length(glaciers),
            MA(ii,:) = movmean(glaciers(ii).CESM2WACCM.ssp126.melt,[sb,sf]);
            MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
        end
    end
    
    MA_N = repmat(MA,N,1);
    sampleretreat = k.*MA_N;

    % take flux-weighted means
    for l=1:8,
        for j=1:N,
            fluxweightedmean(l,j,:) = nansum(iceflux(ids(l).inds)'.*sampleretreat(ids(l).inds+(j-1)*length(glaciers),:))/nansum(iceflux(ids(l).inds));
        end
    end

    % identify the ensemble members which have quartile retreats at 2100
    % and extract retreat timeseries for these members
    for l=1:8,
        fwm0 = fluxweightedmean(l,:,end);
        fwm0_sorted = sort(fwm0);
        id05(l) = min(find(fwm0==fwm0_sorted(floor(0.05*length(fwm0_sorted)))));
        id25(l) = min(find(fwm0==fwm0_sorted(floor(0.25*length(fwm0_sorted)))));
        id50(l) = min(find(fwm0==fwm0_sorted(floor(0.50*length(fwm0_sorted)))));
        id75(l) = min(find(fwm0==fwm0_sorted(floor(0.75*length(fwm0_sorted)))));
        id95(l) = min(find(fwm0==fwm0_sorted(floor(0.95*length(fwm0_sorted)))));
        fw50(l,:) = fluxweightedmean(l,id50(l),:);
        if fw50(l,end)<=0,
            fw05(l,:) = fluxweightedmean(l,id05(l),:);
            fw25(l,:) = fluxweightedmean(l,id25(l),:);
            fw75(l,:) = fluxweightedmean(l,id75(l),:);
            fw95(l,:) = fluxweightedmean(l,id95(l),:);
        elseif fw50(l,end)>0,
            fw95(l,:) = fluxweightedmean(l,id05(l),:);
            fw75(l,:) = fluxweightedmean(l,id25(l),:);
            fw25(l,:) = fluxweightedmean(l,id75(l),:);
            fw05(l,:) = fluxweightedmean(l,id95(l),:);
        end
    end

    % assign to output arrays
    Rp05(mm,:,:) = fw05';
    Rhigh(mm,:,:) = fw25';
    Rmed(mm,:,:) = fw50';
    Rlow(mm,:,:) = fw75';
    Rp95(mm,:,:) = fw95';

end

% permute into more logical order
Rp05 = permute(Rp05,[1,3,2]);
Rhigh = permute(Rhigh,[1,3,2]);
Rmed = permute(Rmed,[1,3,2]);
Rlow = permute(Rlow,[1,3,2]);
Rp95 = permute(Rp95,[1,3,2]);

% assign to output arrays
retreat.time = t;
% CESM2 ssp585
retreat.CESM2WACCM.ssp585.p95 = squeeze(Rp95(1,:,:));
retreat.CESM2WACCM.ssp585.low = squeeze(Rlow(1,:,:));
retreat.CESM2WACCM.ssp585.med = squeeze(Rmed(1,:,:));
retreat.CESM2WACCM.ssp585.high = squeeze(Rhigh(1,:,:));
retreat.CESM2WACCM.ssp585.p05 = squeeze(Rp05(1,:,:));
% CESM2 ssp126
retreat.CESM2WACCM.ssp126.p95 = squeeze(Rp95(2,:,:));
retreat.CESM2WACCM.ssp126.low = squeeze(Rlow(2,:,:));
retreat.CESM2WACCM.ssp126.med = squeeze(Rmed(2,:,:));
retreat.CESM2WACCM.ssp126.high = squeeze(Rhigh(2,:,:));
retreat.CESM2WACCM.ssp126.p05 = squeeze(Rp05(2,:,:));

%% add sector TF for models with ice shelves
for l=1:7,    
    retreat.CESM2WACCM.ssp585.TF(l,:) = interp1(glaciers(ids(l).inds(1)).CESM2WACCM.ssp585.tTF,glaciers(ids(l).inds(1)).CESM2WACCM.ssp585.TF,retreat.time);
    retreat.CESM2WACCM.ssp126.TF(l,:) = interp1(glaciers(ids(l).inds(1)).CESM2WACCM.ssp126.tTF,glaciers(ids(l).inds(1)).CESM2WACCM.ssp126.TF,retreat.time);
end

%% save outputs

load ../final_region_def/ice_ocean_sectors.mat

retreat.regions = regions;
for l=1:7,
    retreat.regions(l).name = sectorname{l};
end

save projected_retreat_SDBN1.mat retreat

%% plot results
if (plot_flg)
    plot_projections_SDBN1
end
