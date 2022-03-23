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
load ../glaciers/glaciers_MARv312.mat
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
% mm = 1 for ACCESS RCP8.5
% mm = 2 for CESM2
% mm = 3 for CNRM-CM6-1 ssp585
% mm = 4 for CNRM-ESM2-1 ssp585
% mm = 5 for MPI-ESM1-2-H ssp585
% mm = 6 for MPI-ESM1-2-H ssp245
% mm = 7 for MPI-ESM1-2-H ssp126
% mm = 8 for UKESM1-0-LL ssp585

modelscenario = {'ACCESS_RCP8.5','CESM2_ssp585','CNRM-CM6-1_ssp585','MPIESM12HR-ssp585','MPIESM12HR-ssp245','MPIESM12HR-ssp126','UKESM1-0-LL-ssp585'};
sectorname = {'SE','SW','CE','CW','NE','NW','NO'};
for mm = [1, 2, 3, 5, 6, 7]

    % make array of projected sample retreat
    if mm == 1,
        t = glaciers(1).ACCESS.RCP85.tmelt;
        for ii=1:length(glaciers),
            MA(ii,:) = movmean(glaciers(ii).ACCESS.RCP85.melt,[sb,sf]);
            MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
        end
    elseif mm == 2,
        t = glaciers(1).CESM2.ssp585.tmelt;
        for ii=1:length(glaciers),
            MA(ii,:) = movmean(glaciers(ii).CESM2.ssp585.melt,[sb,sf]);
            MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
        end
    elseif mm == 3,
        t = glaciers(1).CNRMCM6.ssp585.tmelt;
        for ii=1:length(glaciers),
            MA(ii,:) = movmean(glaciers(ii).CNRMCM6.ssp585.melt,[sb,sf]);
            MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
        end
    elseif mm == 4,
        t = glaciers(1).CNRMESM2.ssp585.tmelt;
        for ii=1:length(glaciers),
            MA(ii,:) = movmean(glaciers(ii).CNRMESM2.ssp585.melt,[sb,sf]);
            MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
        end
    elseif mm == 5,
        t = glaciers(1).MPIESM12HR.ssp585.tmelt;
        for ii=1:length(glaciers),
            MA(ii,:) = movmean(glaciers(ii).MPIESM12HR.ssp585.melt,[sb,sf]);
            MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
        end
    elseif mm == 6,
        t = glaciers(1).MPIESM12HR.ssp245.tmelt;
        for ii=1:length(glaciers),
            MA(ii,:) = movmean(glaciers(ii).MPIESM12HR.ssp245.melt,[sb,sf]);
            MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
        end
    elseif mm == 7,
        t = glaciers(1).MPIESM12HR.ssp126.tmelt;
        for ii=1:length(glaciers),
            MA(ii,:) = movmean(glaciers(ii).MPIESM12HR.ssp126.melt,[sb,sf]);
            MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
        end
    elseif mm == 8,
        t = glaciers(1).UKESM1.ssp585.tmelt;
        for ii=1:length(glaciers),
            MA(ii,:) = movmean(glaciers(ii).UKESM1.ssp585.melt,[sb,sf]);
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
    Rhigh(mm,:,:) = fw25';
    Rmed(mm,:,:) = fw50';
    Rlow(mm,:,:) = fw75';

end

% permute into more logical order
Rhigh = permute(Rhigh,[1,3,2]);
Rmed = permute(Rmed,[1,3,2]);
Rlow = permute(Rlow,[1,3,2]);

% assign to output arrays
retreat.time = t;
% ACCESS RCP8.5
retreat.ACCESS.RCP85.low = squeeze(Rlow(1,:,:));
retreat.ACCESS.RCP85.med = squeeze(Rmed(1,:,:));
retreat.ACCESS.RCP85.high = squeeze(Rhigh(1,:,:));
% CESM2 ssp585
retreat.CESM2.ssp585.low = squeeze(Rlow(2,:,:));
retreat.CESM2.ssp585.med = squeeze(Rmed(2,:,:));
retreat.CESM2.ssp585.high = squeeze(Rhigh(2,:,:));
% CNRM-CM6-1 ssp585
retreat.CNRMCM6.ssp585.low = squeeze(Rlow(3,:,:));
retreat.CNRMCM6.ssp585.med = squeeze(Rmed(3,:,:));
retreat.CNRMCM6.ssp585.high = squeeze(Rhigh(3,:,:));
%% CNRM-ESM2-1 ssp585
%retreat.CNRMESM2.ssp585.low = squeeze(Rlow(4,:,:));
%retreat.CNRMESM2.ssp585.med = squeeze(Rmed(4,:,:));
%retreat.CNRMESM2.ssp585.high = squeeze(Rhigh(4,:,:));
% MPI-ESM1-2-HR ssp585
retreat.MPIESM12HR.ssp585.low = squeeze(Rlow(5,:,:));
retreat.MPIESM12HR.ssp585.med = squeeze(Rmed(5,:,:));
retreat.MPIESM12HR.ssp585.high = squeeze(Rhigh(5,:,:));
% MPI-ESM1-2-HR ssp245
retreat.MPIESM12HR.ssp245.low = squeeze(Rlow(6,:,:));
retreat.MPIESM12HR.ssp245.med = squeeze(Rmed(6,:,:));
retreat.MPIESM12HR.ssp245.high = squeeze(Rhigh(6,:,:));
% MPI-ESM1-2-HR ssp126
retreat.MPIESM12HR.ssp126.low = squeeze(Rlow(7,:,:));
retreat.MPIESM12HR.ssp126.med = squeeze(Rmed(7,:,:));
retreat.MPIESM12HR.ssp126.high = squeeze(Rhigh(7,:,:));
%% UKESM1-0-LL ssp585
%retreat.UKESM1.ssp585.low = squeeze(Rlow(8,:,:));
%retreat.UKESM1.ssp585.med = squeeze(Rmed(8,:,:));
%retreat.UKESM1.ssp585.high = squeeze(Rhigh(8,:,:));

%% add sector TF for models with ice shelves
for l=1:7,    
    retreat.ACCESS.RCP85.TF(l,:) = interp1(glaciers(ids(l).inds(1)).ACCESS.RCP85.tTF,glaciers(ids(l).inds(1)).ACCESS.RCP85.TF,retreat.time);
    retreat.CESM2.ssp585.TF(l,:) = interp1(glaciers(ids(l).inds(1)).CESM2.ssp585.tTF,glaciers(ids(l).inds(1)).CESM2.ssp585.TF,retreat.time);
    retreat.CNRMCM6.ssp585.TF(l,:) = interp1(glaciers(ids(l).inds(1)).CNRMCM6.ssp585.tTF,glaciers(ids(l).inds(1)).CNRMCM6.ssp585.TF,retreat.time);
    %retreat.CNRMESM2.ssp585.TF(l,:) = interp1(glaciers(ids(l).inds(1)).CNRMESM2.ssp585.tTF,glaciers(ids(l).inds(1)).CNRMESM2.ssp585.TF,retreat.time);
    retreat.MPIESM12HR.ssp585.TF(l,:) = interp1(glaciers(ids(l).inds(1)).MPIESM12HR.ssp585.tTF,glaciers(ids(l).inds(1)).MPIESM12HR.ssp585.TF,retreat.time);
    retreat.MPIESM12HR.ssp245.TF(l,:) = interp1(glaciers(ids(l).inds(1)).MPIESM12HR.ssp245.tTF,glaciers(ids(l).inds(1)).MPIESM12HR.ssp245.TF,retreat.time);
    retreat.MPIESM12HR.ssp126.TF(l,:) = interp1(glaciers(ids(l).inds(1)).MPIESM12HR.ssp126.tTF,glaciers(ids(l).inds(1)).MPIESM12HR.ssp126.TF,retreat.time);
    %retreat.UKESM1.ssp585.TF(l,:) = interp1(glaciers(ids(l).inds(1)).UKESM1.ssp585.tTF,glaciers(ids(l).inds(1)).UKESM1.ssp585.TF,retreat.time);
end

%% save outputs

load ../final_region_def/ice_ocean_sectors.mat

retreat.regions = regions;
for l=1:7,
    retreat.regions(l).name = sectorname{l};
end

save projected_retreat_MARv312.mat retreat

%% plot results
if (plot_flg)
    plot_projections_MARv312
end

