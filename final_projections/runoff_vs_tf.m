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

plot_flg = 0;

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
% mm =  1 for ACCESS RCP8.5
% mm =  2 for CESM2
% mm =  3 for CNRM-CM6-1 ssp585
% mm =  4 for CNRM-ESM2-1 ssp585
% mm =  5 for MPI-ESM1-2-H ssp585
% mm =  6 for MPI-ESM1-2-H ssp245
% mm =  7 for MPI-ESM1-2-H ssp126
% mm =  8 for UKESM1-0-LL ssp585
% mm =  9 for NorESM2 - ssp585
% mm = 10 for NorESM2 - ssp245
% mm = 11 for CESM2-CMIP6 - ssp585
% mm = 12 for CESM2-CMIP6 - ssp245
% mm = 13 for CESM2-CMIP6 - ssp126
% mm = 14 for UKESM1-0-LL-CMIP6 - ssp585
% mm = 15 for UKESM1-0-LL-CMIP6 - ssp245
% mm = 16 for IPSL-CM6A-LR - ssp585

sectorname = {'SE','SW','CE','CW','NE','NW','NO'};
%for mm = [1, 2, 3, 4, 5, 6, 7, 8]
for mm = [1:16]

    % make array of projected sample retreat
    if mm == 1,
        t = glaciers(1).ACCESS.RCP85.tmelt;
        tq = glaciers(1).ACCESS.RCP85.tJJA;
        for ii=1:length(glaciers),
            %MA(ii,:) = movmean(glaciers(ii).ACCESS.RCP85.melt,[sb,sf]);
            %MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
            QA(ii,:) = movmean(glaciers(ii).ACCESS.RCP85.QJJA,[sb,sf]);
            QA(ii,:) = QA(ii,:) - QA(ii,find(floor(tq)==yr0));
        end
    elseif mm == 2,
        t = glaciers(1).CESM2Leo.ssp585.tmelt;
        tq = glaciers(1).CESM2Leo.ssp585.tJJA;
        for ii=1:length(glaciers),
            %MA(ii,:) = movmean(glaciers(ii).CESM2Leo.ssp585.melt,[sb,sf]);
            %MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
            QA(ii,:) = movmean(glaciers(ii).CESM2Leo.ssp585.QJJA,[sb,sf]);
            QA(ii,:) = QA(ii,:) - QA(ii,find(floor(tq)==yr0));
        end
    elseif mm == 3,
        t = glaciers(1).CNRMCM6.ssp585.tmelt;
        tq = glaciers(1).CNRMCM6.ssp585.tJJA;
        for ii=1:length(glaciers),
            %MA(ii,:) = movmean(glaciers(ii).CNRMCM6.ssp585.melt,[sb,sf]);
            %MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
            QA(ii,:) = movmean(glaciers(ii).CNRMCM6.ssp585.QJJA,[sb,sf]);
            QA(ii,:) = QA(ii,:) - QA(ii,find(floor(tq)==yr0));
        end
    elseif mm == 4,
        t = glaciers(1).CNRMESM2.ssp585.tmelt;
        tq = glaciers(1).CNRMESM2.ssp585.tJJA;
        for ii=1:length(glaciers),
            %MA(ii,:) = movmean(glaciers(ii).CNRMESM2.ssp585.melt,[sb,sf]);
            %MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
            QA(ii,:) = movmean(glaciers(ii).CNRMESM2.ssp585.QJJA,[sb,sf]);
            QA(ii,:) = QA(ii,:) - QA(ii,find(floor(tq)==yr0));
        end
    elseif mm == 5,
        t = glaciers(1).MPIESM12HR.ssp585.tmelt;
        tq = glaciers(1).MPIESM12HR.ssp585.tJJA;
        for ii=1:length(glaciers),
            %MA(ii,:) = movmean(glaciers(ii).MPIESM12HR.ssp585.melt,[sb,sf]);
            %MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
            QA(ii,:) = movmean(glaciers(ii).MPIESM12HR.ssp585.QJJA,[sb,sf]);
            QA(ii,:) = QA(ii,:) - QA(ii,find(floor(tq)==yr0));
        end
    elseif mm == 6,
        t = glaciers(1).MPIESM12HR.ssp245.tmelt;
        tq = glaciers(1).MPIESM12HR.ssp245.tJJA;
        for ii=1:length(glaciers),
            %MA(ii,:) = movmean(glaciers(ii).MPIESM12HR.ssp245.melt,[sb,sf]);
            %MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
            QA(ii,:) = movmean(glaciers(ii).MPIESM12HR.ssp245.QJJA,[sb,sf]);
            QA(ii,:) = QA(ii,:) - QA(ii,find(floor(tq)==yr0));
        end
    elseif mm == 7,
        t = glaciers(1).MPIESM12HR.ssp126.tmelt;
        tq = glaciers(1).MPIESM12HR.ssp126.tJJA;
        for ii=1:length(glaciers),
            %MA(ii,:) = movmean(glaciers(ii).MPIESM12HR.ssp126.melt,[sb,sf]);
            %MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
            QA(ii,:) = movmean(glaciers(ii).MPIESM12HR.ssp126.QJJA,[sb,sf]);
            QA(ii,:) = QA(ii,:) - QA(ii,find(floor(tq)==yr0));
        end
    elseif mm == 8,
        % needs padding for missing values 1950-1959, repeat 1960-1969
        %t = glaciers(1).UKESM1Robin.ssp585.tmelt;
        t = [glaciers(1).UKESM1Robin.ssp585.tmelt(1:10)-10, glaciers(1).UKESM1Robin.ssp585.tmelt];
        tq = [glaciers(1).UKESM1Robin.ssp585.tJJA(1:10)-10, glaciers(1).UKESM1Robin.ssp585.tJJA];
        for ii=1:length(glaciers),
            %ukm = [glaciers(ii).UKESM1Robin.ssp585.melt(1:10), glaciers(ii).UKESM1Robin.ssp585.melt];
            %MA(ii,:) = movmean(ukm,[sb,sf]);
            %MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
            ukm = [glaciers(ii).UKESM1Robin.ssp585.QJJA(1:10), glaciers(ii).UKESM1Robin.ssp585.QJJA];
            QA(ii,:) = movmean(ukm,[sb,sf]);
            QA(ii,:) = QA(ii,:) - QA(ii,find(floor(tq)==yr0));
        end
    elseif mm == 9,
        t = glaciers(1).NorESM2.ssp585.tmelt;
        tq = glaciers(1).NorESM2.ssp585.tJJA;
        for ii=1:length(glaciers),
            %MA(ii,:) = movmean(glaciers(ii).NorESM2.ssp585.melt,[sb,sf]);
            %MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
            QA(ii,:) = movmean(glaciers(ii).NorESM2.ssp585.QJJA,[sb,sf]);
            QA(ii,:) = QA(ii,:) - QA(ii,find(floor(tq)==yr0));
        end
    elseif mm == 10,
        t = glaciers(1).NorESM2.ssp245.tmelt;
        tq = glaciers(1).NorESM2.ssp245.tJJA;
        for ii=1:length(glaciers),
            %MA(ii,:) = movmean(glaciers(ii).NorESM2.ssp245.melt,[sb,sf]);
            %MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
            QA(ii,:) = movmean(glaciers(ii).NorESM2.ssp245.QJJA,[sb,sf]);
            QA(ii,:) = QA(ii,:) - QA(ii,find(floor(tq)==yr0));
        end
    elseif mm == 11,
        t = glaciers(1).CESM2CMIP6.ssp585.tmelt;
        tq = glaciers(1).CESM2CMIP6.ssp585.tJJA;
        for ii=1:length(glaciers),
            %MA(ii,:) = movmean(glaciers(ii).CESM2CMIP6.ssp585.melt,[sb,sf]);
            %MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
            QA(ii,:) = movmean(glaciers(ii).CESM2CMIP6.ssp585.QJJA,[sb,sf]);
            QA(ii,:) = QA(ii,:) - QA(ii,find(floor(tq)==yr0));
        end
    elseif mm == 12,
        t = glaciers(1).CESM2CMIP6.ssp245.tmelt;
        tq = glaciers(1).CESM2CMIP6.ssp245.tJJA;
        for ii=1:length(glaciers),
            %MA(ii,:) = movmean(glaciers(ii).CESM2CMIP6.ssp245.melt,[sb,sf]);
            %MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
            QA(ii,:) = movmean(glaciers(ii).CESM2CMIP6.ssp245.QJJA,[sb,sf]);
            QA(ii,:) = QA(ii,:) - QA(ii,find(floor(tq)==yr0));
        end
    elseif mm == 13,
        t = glaciers(1).CESM2CMIP6.ssp126.tmelt;
        tq = glaciers(1).CESM2CMIP6.ssp126.tJJA;
        for ii=1:length(glaciers),
            %MA(ii,:) = movmean(glaciers(ii).CESM2CMIP6.ssp126.melt,[sb,sf]);
            %MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
            QA(ii,:) = movmean(glaciers(ii).CESM2CMIP6.ssp126.QJJA,[sb,sf]);
            QA(ii,:) = QA(ii,:) - QA(ii,find(floor(tq)==yr0));
        end
    elseif mm == 14,
        t = glaciers(1).UKESM1CMIP6.ssp585.tmelt;
        tq = glaciers(1).UKESM1CMIP6.ssp585.tJJA;
        for ii=1:length(glaciers),
            %MA(ii,:) = movmean(glaciers(ii).UKESM1CMIP6.ssp585.melt,[sb,sf]);
            %MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
            QA(ii,:) = movmean(glaciers(ii).UKESM1CMIP6.ssp585.QJJA,[sb,sf]);
            QA(ii,:) = QA(ii,:) - QA(ii,find(floor(tq)==yr0));
        end
    elseif mm == 15,
        t = glaciers(1).UKESM1CMIP6.ssp245.tmelt;
        tq = glaciers(1).UKESM1CMIP6.ssp245.tJJA;
        for ii=1:length(glaciers),
            %MA(ii,:) = movmean(glaciers(ii).UKESM1CMIP6.ssp245.melt,[sb,sf]);
            %MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
            QA(ii,:) = movmean(glaciers(ii).UKESM1CMIP6.ssp245.QJJA,[sb,sf]);
            QA(ii,:) = QA(ii,:) - QA(ii,find(floor(tq)==yr0));
        end
    elseif mm == 16,
        t = glaciers(1).IPSLCM6ALR.ssp585.tmelt;
        tq = glaciers(1).IPSLCM6ALR.ssp585.tJJA;
        for ii=1:length(glaciers),
            %MA(ii,:) = movmean(glaciers(ii).IPSLCM6ALR.ssp585.melt,[sb,sf]);
            %MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
            QA(ii,:) = movmean(glaciers(ii).IPSLCM6ALR.ssp585.QJJA,[sb,sf]);
            QA(ii,:) = QA(ii,:) - QA(ii,find(floor(tq)==yr0));
        end
    end

    %MA_N = repmat(MA,N,1);
    %sampleretreat = k.*MA_N;

    QA_N = repmat(QA,N,1);
    sampleQ = k.*QA_N;

    % take flux-weighted means
    for l=1:8,
        for j=1:N,
            %fluxweightedmean(l,j,:) = nansum(iceflux(ids(l).inds)'.*sampleretreat(ids(l).inds+(j-1)*length(glaciers),:))/nansum(iceflux(ids(l).inds));

            fwq(l,j,:) = nansum(iceflux(ids(l).inds)'.*sampleQ(ids(l).inds+(j-1)*length(glaciers),:))/nansum(iceflux(ids(l).inds));
        end
    end

    % identify the ensemble members which have quartile retreats at 2100
    % and extract retreat timeseries for these members
    for l=1:8,
        %fwm0 = fluxweightedmean(l,:,end);
        %fwm0_sorted = sort(fwm0);
        %id50(l) = min(find(fwm0==fwm0_sorted(floor(0.50*length(fwm0_sorted)))));
        %fw50(l,:) = fluxweightedmean(l,id50(l),:);
    end

    for l=1:8,
        fwq0 = fwq(l,:,end);
        fwq0_sorted = sort(fwq0);
        idq50(l) = min(find(fwq0==fwq0_sorted(floor(0.50*length(fwq0_sorted)))));
        fwq50(l,:) = fwq(l,idq50(l),:);
    end

    % assign to output arrays
    %Rmed(mm,:,:) = fw50';

    % assign to output arrays
    Qmed(mm,:,:) = fwq50';

end

% permute into more logical order
%Rmed = permute(Rmed,[1,3,2]);

Qmed = permute(Qmed,[1,3,2]);

%% assign to output arrays
%retreat.time = t;
%% ACCESS RCP8.5
%retreat.ACCESS.RCP85.med = squeeze(Rmed(1,:,:));
%% CESM2 ssp585
%retreat.CESM2Leo.ssp585.med = squeeze(Rmed(2,:,:));
%% CNRM-CM6-1 ssp585
%retreat.CNRMCM6.ssp585.med = squeeze(Rmed(3,:,:));
%% CNRM-ESM2-1 ssp585
%retreat.CNRMESM2.ssp585.med = squeeze(Rmed(4,:,:));
%% MPI-ESM1-2-HR ssp585
%retreat.MPIESM12HR.ssp585.med = squeeze(Rmed(5,:,:));
%% MPI-ESM1-2-HR ssp245
%retreat.MPIESM12HR.ssp245.med = squeeze(Rmed(6,:,:));
%% MPI-ESM1-2-HR ssp126
%retreat.MPIESM12HR.ssp126.med = squeeze(Rmed(7,:,:));
%% UKESM1-0-LL ssp585
%retreat.UKESM1Robin.ssp585.med = squeeze(Rmed(8,:,:));
%% NorESM2 ssp585
%retreat.NorESM2.ssp585.med = squeeze(Rmed(9,:,:));
%% NorESM2 ssp245
%retreat.NorESM2.ssp245.med = squeeze(Rmed(10,:,:));
%% CESM2-CMIP6 - ssp585
%retreat.CESM2CMIP6.ssp585.med = squeeze(Rmed(11,:,:));
%% CESM2-CMIP6 - ssp245
%retreat.CESM2CMIP6.ssp245.med = squeeze(Rmed(12,:,:));
%% CESM2-CMIP6 - ssp126
%retreat.CESM2CMIP6.ssp126.med = squeeze(Rmed(13,:,:));
%% UKESM1-0-LL-CMIP6 - ssp585
%retreat.UKESM1CMIP6.ssp585.med = squeeze(Rmed(14,:,:));
%% UKESM1-0-LL-CMIP6 - ssp245
%retreat.UKESM1CMIP6.ssp245.med = squeeze(Rmed(15,:,:));
%% IPSL-CM6A-LR - ssp585
%retreat.IPSLCM6ALR.ssp585.med = squeeze(Rmed(16,:,:));

% assign to output arrays
retreat.time = t;
% ACCESS RCP8.5
retreat.ACCESS.RCP85.Qmed = squeeze(Qmed(1,:,:));
% CESM2 ssp585
retreat.CESM2Leo.ssp585.Qmed = squeeze(Qmed(2,:,:));
% CNRM-CM6-1 ssp585
retreat.CNRMCM6.ssp585.Qmed = squeeze(Qmed(3,:,:));
% CNRM-ESM2-1 ssp585
retreat.CNRMESM2.ssp585.Qmed = squeeze(Qmed(4,:,:));
% MPI-ESM1-2-HR ssp585
retreat.MPIESM12HR.ssp585.Qmed = squeeze(Qmed(5,:,:));
% MPI-ESM1-2-HR ssp245
retreat.MPIESM12HR.ssp245.Qmed = squeeze(Qmed(6,:,:));
% MPI-ESM1-2-HR ssp126
retreat.MPIESM12HR.ssp126.Qmed = squeeze(Qmed(7,:,:));
% UKESM1-0-LL ssp585
retreat.UKESM1Robin.ssp585.Qmed = squeeze(Qmed(8,:,:));
% NorESM2 ssp585
retreat.NorESM2.ssp585.Qmed = squeeze(Qmed(9,:,:));
% NorESM2 ssp245
retreat.NorESM2.ssp245.Qmed = squeeze(Qmed(10,:,:));
% CESM2-CMIP6 - ssp585
retreat.CESM2CMIP6.ssp585.Qmed = squeeze(Qmed(11,:,:));
% CESM2-CMIP6 - ssp245
retreat.CESM2CMIP6.ssp245.Qmed = squeeze(Qmed(12,:,:));
% CESM2-CMIP6 - ssp126
retreat.CESM2CMIP6.ssp126.Qmed = squeeze(Qmed(13,:,:));
% UKESM1-0-LL-CMIP6 - ssp585
retreat.UKESM1CMIP6.ssp585.Qmed = squeeze(Qmed(14,:,:));
% UKESM1-0-LL-CMIP6 - ssp245
retreat.UKESM1CMIP6.ssp245.Qmed = squeeze(Qmed(15,:,:));
% IPSL-CM6A-LR - ssp585
retreat.IPSLCM6ALR.ssp585.Qmed = squeeze(Qmed(16,:,:));

%% add sector TF for models with ice shelves
for l=1:7,    
    retreat.ACCESS.RCP85.TF(l,:) = interp1(glaciers(ids(l).inds(1)).ACCESS.RCP85.tTF,glaciers(ids(l).inds(1)).ACCESS.RCP85.TF,retreat.time);
    retreat.CESM2Leo.ssp585.TF(l,:) = interp1(glaciers(ids(l).inds(1)).CESM2Leo.ssp585.tTF,glaciers(ids(l).inds(1)).CESM2Leo.ssp585.TF,retreat.time);
    retreat.CNRMCM6.ssp585.TF(l,:) = interp1(glaciers(ids(l).inds(1)).CNRMCM6.ssp585.tTF,glaciers(ids(l).inds(1)).CNRMCM6.ssp585.TF,retreat.time);
    retreat.CNRMESM2.ssp585.TF(l,:) = interp1(glaciers(ids(l).inds(1)).CNRMESM2.ssp585.tTF,glaciers(ids(l).inds(1)).CNRMESM2.ssp585.TF,retreat.time);

    retreat.MPIESM12HR.ssp585.TF(l,:) = interp1(glaciers(ids(l).inds(1)).MPIESM12HR.ssp585.tTF,glaciers(ids(l).inds(1)).MPIESM12HR.ssp585.TF,retreat.time);
    retreat.MPIESM12HR.ssp245.TF(l,:) = interp1(glaciers(ids(l).inds(1)).MPIESM12HR.ssp245.tTF,glaciers(ids(l).inds(1)).MPIESM12HR.ssp245.TF,retreat.time);
    retreat.MPIESM12HR.ssp126.TF(l,:) = interp1(glaciers(ids(l).inds(1)).MPIESM12HR.ssp126.tTF,glaciers(ids(l).inds(1)).MPIESM12HR.ssp126.TF,retreat.time);

    retreat.UKESM1Robin.ssp585.TF(l,:) = interp1(glaciers(ids(l).inds(1)).UKESM1Robin.ssp585.tTF,glaciers(ids(l).inds(1)).UKESM1Robin.ssp585.TF,retreat.time);

    retreat.NorESM2.ssp585.TF(l,:) = interp1(glaciers(ids(l).inds(1)).NorESM2.ssp585.tTF,glaciers(ids(l).inds(1)).NorESM2.ssp585.TF,retreat.time);
    retreat.NorESM2.ssp245.TF(l,:) = interp1(glaciers(ids(l).inds(1)).NorESM2.ssp245.tTF,glaciers(ids(l).inds(1)).NorESM2.ssp245.TF,retreat.time);

    retreat.CESM2CMIP6.ssp585.TF(l,:) = interp1(glaciers(ids(l).inds(1)).CESM2CMIP6.ssp585.tTF,glaciers(ids(l).inds(1)).CESM2CMIP6.ssp585.TF,retreat.time);
    retreat.CESM2CMIP6.ssp245.TF(l,:) = interp1(glaciers(ids(l).inds(1)).CESM2CMIP6.ssp245.tTF,glaciers(ids(l).inds(1)).CESM2CMIP6.ssp245.TF,retreat.time);
    retreat.CESM2CMIP6.ssp126.TF(l,:) = interp1(glaciers(ids(l).inds(1)).CESM2CMIP6.ssp126.tTF,glaciers(ids(l).inds(1)).CESM2CMIP6.ssp126.TF,retreat.time);
    
    retreat.UKESM1CMIP6.ssp585.TF(l,:) = interp1(glaciers(ids(l).inds(1)).UKESM1CMIP6.ssp585.tTF,glaciers(ids(l).inds(1)).UKESM1CMIP6.ssp585.TF,retreat.time);
    retreat.UKESM1CMIP6.ssp245.TF(l,:) = interp1(glaciers(ids(l).inds(1)).UKESM1CMIP6.ssp245.tTF,glaciers(ids(l).inds(1)).UKESM1CMIP6.ssp245.TF,retreat.time);

    retreat.IPSLCM6ALR.ssp585.TF(l,:) = interp1(glaciers(ids(l).inds(1)).IPSLCM6ALR.ssp585.tTF,glaciers(ids(l).inds(1)).IPSLCM6ALR.ssp585.TF,retreat.time);
end

% save outputs

load ../final_region_def/ice_ocean_sectors.mat

retreat.regions = regions;
for l=1:7,
    retreat.regions(l).name = sectorname{l};
end

save projected_retreat_MARv312.mat retreat

%% plot results
if (plot_flg)
    plot_rtf_MARv312
end

