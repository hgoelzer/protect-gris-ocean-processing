clear; close all;
% script to add future runoff, ocean TF and melt to glaciers structure
load glaciers_past.mat

% modelid   1 ACCESS1.3 - rcp85
% modelid   2 CESM2-Leo - ssp585
% modelid   3 CNRM-CM6 - ssp585
% modelid   4 CNRM-ESM2 - ssp585 
% modelid   5 MPIESM12HR - ssp585
% modelid   6 MPIESM12HR - ssp245
% modelid   7 MPIESM12HR - ssp126
% modelid   8 UKESM1-0-LL-Robin - ssp585
% modelid   9 NorESM2 - ssp585
% modelid  10 NorESM2 - ssp245
% modelid  11 CESM2-CMIP6 - ssp585
% modelid  12 CESM2-CMIP6 - ssp245
% modelid  13 CESM2-CMIP6 - ssp126
% modelid  14 UKESM1-0-LL-CMIP6 - ssp245
% modelid  15 UKESM1-0-LL-CMIP6 - ssp585
% modelid  16 IPSL-CM6A-LR - ssp585

%modelid = [2];
%modelid = [8];
%modelid = [1, 2, 3, 4, 5, 6, 7, 8];
modelid = [1:16];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(modelid == 1)
disp('Adding ACCESS1.3 - rcp85')

%% associate future thermal forcing from ACCESS RCP8.5 (note there is no ACCESS RCP2.6)

% load ACCESS RCP8.5 ocean
load ../process_CMIP/ACCESS_RCP85.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
ACCESS_baseline_inds = find(ismember(year,baseline));
TF0_ACCESS = nanmean(squeeze(TF_basins(:,ACCESS_baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% ACCESS bias
bias = TF0_ACCESS - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0_ACCESS(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).ACCESS.RCP85.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).ACCESS.RCP85.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).ACCESS.RCP85.bias_TF = bias(glaciers(ii).sectornum);
end

%% associate future runoff from MAR/ACCESS

% load ACCESS-1-3 RCP8.5 run
load ../runoff/MARv3.12-ACCESS1.3-rcp85-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat

runoff_mod = runoff;

% sort time variable
for ii=1:length(runoff_mod),
    runoff_mod(ii).time = [1950:2100];
end

% assign to glaciers
for ii=1:length(glaciers),
    % rcp8.5
    for jj=1:length(runoff_mod),
        if ismember(glaciers(ii).rignotid,runoff_mod(jj).rignotGlacierID),
            glaciers(ii).ACCESS.RCP85.tJJA = runoff_mod(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).ACCESS.RCP85.tJJA),
                yr = floor(glaciers(ii).ACCESS.RCP85.tJJA(kk));
                glaciers(ii).ACCESS.RCP85.QJJA(kk) = 3*runoff_mod(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).ACCESS.RCP85.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).ACCESS.RCP85.bias_QJJA = mean(glaciers(ii).ACCESS.RCP85.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).ACCESS.RCP85.QJJA = glaciers(ii).ACCESS.RCP85.QJJA - glaciers(ii).ACCESS.RCP85.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).ACCESS.RCP85.QJJA(find(glaciers(ii).ACCESS.RCP85.QJJA<0)) = 0;
end

%% calculate future ACCESS melt

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % RCP8.5
    glaciers(ii).ACCESS.RCP85.tmelt = glaciers(ii).ACCESS.RCP85.tJJA;
    glaciers(ii).ACCESS.RCP85.melt = (glaciers(ii).ACCESS.RCP85.QJJA.^0.4).*...
        interp1(glaciers(ii).ACCESS.RCP85.tTF,glaciers(ii).ACCESS.RCP85.TF,glaciers(ii).ACCESS.RCP85.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).ACCESS.RCP85.melt(find(glaciers(ii).ACCESS.RCP85.melt<0)) = 0;
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(modelid == 2)
disp('Adding CESM2-Leo - ssp585')
%% associate future thermal forcing from CESM2-Leo ssp585

% load CESM2-Leo ssp585 ocean
load ../process_CMIP/CESM2-Leo_ssp585.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
cesm2_baseline_inds = find(ismember(year,baseline));
TF0_cesm2 = nanmean(squeeze(TF_basins(:,cesm2_baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% cesm2 bias
bias = TF0_cesm2 - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0_cesm2(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).CESM2Leo.ssp585.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).CESM2Leo.ssp585.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).CESM2Leo.ssp585.bias_TF = bias(glaciers(ii).sectornum);
end

%% associate future runoff from MAR/CESM2-Leo ssp585

% load CESM2-Leo/MAR ssp585
load ../runoff/MARv3.12-CESM2-Leo-ssp585-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_mod = runoff;

% sort time variable
for ii=1:length(runoff_mod),
    runoff_mod(ii).time = [1950:2100];
end

% assign to glaciers
for ii=1:length(glaciers),
    % ssp585
    for jj=1:length(runoff_mod),
        if ismember(glaciers(ii).rignotid,runoff_mod(jj).rignotGlacierID),
            glaciers(ii).CESM2Leo.ssp585.tJJA = runoff_mod(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).CESM2Leo.ssp585.tJJA),
                yr = floor(glaciers(ii).CESM2Leo.ssp585.tJJA(kk));
                glaciers(ii).CESM2Leo.ssp585.QJJA(kk) = 3*runoff_mod(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).CESM2Leo.ssp585.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).CESM2Leo.ssp585.bias_QJJA = mean(glaciers(ii).CESM2Leo.ssp585.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).CESM2Leo.ssp585.QJJA = glaciers(ii).CESM2Leo.ssp585.QJJA - glaciers(ii).CESM2Leo.ssp585.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).CESM2Leo.ssp585.QJJA(find(glaciers(ii).CESM2Leo.ssp585.QJJA<0)) = 0;
end

%% calculate future CESM2-Leo melt (ssp585 only)

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % ssp585
    glaciers(ii).CESM2Leo.ssp585.tmelt = glaciers(ii).CESM2Leo.ssp585.tJJA;
    glaciers(ii).CESM2Leo.ssp585.melt = (glaciers(ii).CESM2Leo.ssp585.QJJA.^0.4).*...
        interp1(glaciers(ii).CESM2Leo.ssp585.tTF,glaciers(ii).CESM2Leo.ssp585.TF,glaciers(ii).CESM2Leo.ssp585.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).CESM2Leo.ssp585.melt(find(glaciers(ii).CESM2Leo.ssp585.melt<0)) = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(modelid == 3)
disp('Adding CNRM-CM6-1 ssp585')

%% associate future thermal forcing from CNRM-CM6-1 ssp585

% load CNRM-CM6-1 ssp585 ocean
load ../process_CMIP/CNRM-CM6-1_ssp585.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
cnrmcm6_baseline_inds = find(ismember(year,baseline));
TF0_cnrmcm6 = nanmean(squeeze(TF_basins(:,cnrmcm6_baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% cnrmcm6 bias
bias = TF0_cnrmcm6 - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0_cnrmcm6(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).CNRMCM6.ssp585.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).CNRMCM6.ssp585.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).CNRMCM6.ssp585.bias_TF = bias(glaciers(ii).sectornum);
end

%% associate future runoff from MAR/CNRM-CM6-1 ssp585
% load CNRM-CM6-1/MAR ssp585
load ../runoff/MARv3.12-CNRM-CM6-ssp585-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_mod = runoff;

% sort time variable
for ii=1:length(runoff_mod),
    runoff_mod(ii).time = [1950:2100];
end

% assign to glaciers
for ii=1:length(glaciers),
    % ssp585
    for jj=1:length(runoff_mod),
        if ismember(glaciers(ii).rignotid,runoff_mod(jj).rignotGlacierID),
            glaciers(ii).CNRMCM6.ssp585.tJJA = runoff_mod(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).CNRMCM6.ssp585.tJJA),
                yr = floor(glaciers(ii).CNRMCM6.ssp585.tJJA(kk));
                glaciers(ii).CNRMCM6.ssp585.QJJA(kk) = 3*runoff_mod(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).CNRMCM6.ssp585.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).CNRMCM6.ssp585.bias_QJJA = mean(glaciers(ii).CNRMCM6.ssp585.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).CNRMCM6.ssp585.QJJA = glaciers(ii).CNRMCM6.ssp585.QJJA - glaciers(ii).CNRMCM6.ssp585.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).CNRMCM6.ssp585.QJJA(find(glaciers(ii).CNRMCM6.ssp585.QJJA<0)) = 0;
end

%% calculate future CNRM-CM6-1 melt

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % ssp585
    glaciers(ii).CNRMCM6.ssp585.tmelt = glaciers(ii).CNRMCM6.ssp585.tJJA;
    glaciers(ii).CNRMCM6.ssp585.melt = (glaciers(ii).CNRMCM6.ssp585.QJJA.^0.4).*...
        interp1(glaciers(ii).CNRMCM6.ssp585.tTF,glaciers(ii).CNRMCM6.ssp585.TF,glaciers(ii).CNRMCM6.ssp585.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).CNRMCM6.ssp585.melt(find(glaciers(ii).CNRMCM6.ssp585.melt<0)) = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(modelid == 4)
disp('Adding CNRM-ESM2-1 ssp585')

%% associate future thermal forcing from CNRM-ESM2-1 ssp585

% load CNRM-ESM2-1 ssp585 ocean
load ../process_CMIP/CNRM-ESM2-1_ssp585.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
cnrmesm2_baseline_inds = find(ismember(year,baseline));
TF0_cnrmesm2 = nanmean(squeeze(TF_basins(:,cnrmesm2_baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% cnrmesm2 bias
bias = TF0_cnrmesm2 - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0_cnrmesm2(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).CNRMESM2.ssp585.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).CNRMESM2.ssp585.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).CNRMESM2.ssp585.bias_TF = bias(glaciers(ii).sectornum);
end
%% associate future runoff from MAR/CNRM-ESM2-1 (ssp 585 only)

% load CNRM-ESM2-1/MAR ssp585
load ../runoff/MARv3.12-CNRM-ESM2-ssp585-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_mod = runoff;

% sort time variable
for ii=1:length(runoff_mod),
    runoff_mod(ii).time = [1950:2100];
end

% assign to glaciers
for ii=1:length(glaciers),
    % ssp585
    for jj=1:length(runoff_mod),
        if ismember(glaciers(ii).rignotid,runoff_mod(jj).rignotGlacierID),
            glaciers(ii).CNRMESM2.ssp585.tJJA = runoff_mod(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).CNRMESM2.ssp585.tJJA),
                yr = floor(glaciers(ii).CNRMESM2.ssp585.tJJA(kk));
                glaciers(ii).CNRMESM2.ssp585.QJJA(kk) = 3*runoff_mod(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).CNRMESM2.ssp585.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).CNRMESM2.ssp585.bias_QJJA = mean(glaciers(ii).CNRMESM2.ssp585.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).CNRMESM2.ssp585.QJJA = glaciers(ii).CNRMESM2.ssp585.QJJA - glaciers(ii).CNRMESM2.ssp585.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).CNRMESM2.ssp585.QJJA(find(glaciers(ii).CNRMESM2.ssp585.QJJA<0)) = 0;
end

%% calculate future CNRM-ESM2-1 melt (ssp585 only)

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % ssp585
    glaciers(ii).CNRMESM2.ssp585.tmelt = glaciers(ii).CNRMESM2.ssp585.tJJA;
    glaciers(ii).CNRMESM2.ssp585.melt = (glaciers(ii).CNRMESM2.ssp585.QJJA.^0.4).*...
        interp1(glaciers(ii).CNRMESM2.ssp585.tTF,glaciers(ii).CNRMESM2.ssp585.TF,glaciers(ii).CNRMESM2.ssp585.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).CNRMESM2.ssp585.melt(find(glaciers(ii).CNRMESM2.ssp585.melt<0)) = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(modelid == 5)
disp('Adding MPIESM12HR ssp585')

%% associate future thermal forcing from MPIESM12HR ssp585

% load ocean data
load ../process_CMIP/MPIESM12HR_ssp585.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
cnrmcm6_baseline_inds = find(ismember(year,baseline));
TF0_model = nanmean(squeeze(TF_basins(:,cnrmcm6_baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% bias
bias = TF0_model - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0_model(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).MPIESM12HR.ssp585.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).MPIESM12HR.ssp585.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).MPIESM12HR.ssp585.bias_TF = bias(glaciers(ii).sectornum);
end

%% associate future runoff from MARv3.12/MPIESM12HR ssp585
% load runoff
load ../runoff/MARv3.12-MPI-ESM1-2-HR-ssp585-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_mod = runoff;

% sort time variable
for ii=1:length(runoff_mod),
    runoff_mod(ii).time = [1950:2100];
end

% assign to glaciers
for ii=1:length(glaciers),
    % ssp585
    for jj=1:length(runoff_mod),
        if ismember(glaciers(ii).rignotid,runoff_mod(jj).rignotGlacierID),
            glaciers(ii).MPIESM12HR.ssp585.tJJA = runoff_mod(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).MPIESM12HR.ssp585.tJJA),
                yr = floor(glaciers(ii).MPIESM12HR.ssp585.tJJA(kk));
                glaciers(ii).MPIESM12HR.ssp585.QJJA(kk) = 3*runoff_mod(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).MPIESM12HR.ssp585.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).MPIESM12HR.ssp585.bias_QJJA = mean(glaciers(ii).MPIESM12HR.ssp585.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).MPIESM12HR.ssp585.QJJA = glaciers(ii).MPIESM12HR.ssp585.QJJA - glaciers(ii).MPIESM12HR.ssp585.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).MPIESM12HR.ssp585.QJJA(find(glaciers(ii).MPIESM12HR.ssp585.QJJA<0)) = 0;
end

%% calculate future MPIESM12HR melt ssp585
% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    glaciers(ii).MPIESM12HR.ssp585.tmelt = glaciers(ii).MPIESM12HR.ssp585.tJJA;
    glaciers(ii).MPIESM12HR.ssp585.melt = (glaciers(ii).MPIESM12HR.ssp585.QJJA.^0.4).*...
        interp1(glaciers(ii).MPIESM12HR.ssp585.tTF,glaciers(ii).MPIESM12HR.ssp585.TF,glaciers(ii).MPIESM12HR.ssp585.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).MPIESM12HR.ssp585.melt(find(glaciers(ii).MPIESM12HR.ssp585.melt<0)) = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(modelid == 6)
disp('Adding MPIESM12HR ssp245')

%% associate future thermal forcing from MPIESM12HR ssp245

% load ocean data
load ../process_CMIP/MPIESM12HR_ssp245.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
cnrmcm6_baseline_inds = find(ismember(year,baseline));
TF0_model = nanmean(squeeze(TF_basins(:,cnrmcm6_baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% bias
bias = TF0_model - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0_model(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).MPIESM12HR.ssp245.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).MPIESM12HR.ssp245.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).MPIESM12HR.ssp245.bias_TF = bias(glaciers(ii).sectornum);
end

%% associate future runoff from MARv3.12/MPIESM12HR ssp245
% load runoff
load ../runoff/MARv3.12-MPI-ESM1-2-HR-ssp245-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_mod = runoff;

% sort time variable
for ii=1:length(runoff_mod),
    runoff_mod(ii).time = [1950:2100];
end

% assign to glaciers
for ii=1:length(glaciers),
    % ssp585
    for jj=1:length(runoff_mod),
        if ismember(glaciers(ii).rignotid,runoff_mod(jj).rignotGlacierID),
            glaciers(ii).MPIESM12HR.ssp245.tJJA = runoff_mod(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).MPIESM12HR.ssp245.tJJA),
                yr = floor(glaciers(ii).MPIESM12HR.ssp245.tJJA(kk));
                glaciers(ii).MPIESM12HR.ssp245.QJJA(kk) = 3*runoff_mod(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).MPIESM12HR.ssp245.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).MPIESM12HR.ssp245.bias_QJJA = mean(glaciers(ii).MPIESM12HR.ssp245.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).MPIESM12HR.ssp245.QJJA = glaciers(ii).MPIESM12HR.ssp245.QJJA - glaciers(ii).MPIESM12HR.ssp245.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).MPIESM12HR.ssp245.QJJA(find(glaciers(ii).MPIESM12HR.ssp245.QJJA<0)) = 0;
end

%% calculate future MPIESM12HR melt ssp245
% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    glaciers(ii).MPIESM12HR.ssp245.tmelt = glaciers(ii).MPIESM12HR.ssp245.tJJA;
    glaciers(ii).MPIESM12HR.ssp245.melt = (glaciers(ii).MPIESM12HR.ssp245.QJJA.^0.4).*...
        interp1(glaciers(ii).MPIESM12HR.ssp245.tTF,glaciers(ii).MPIESM12HR.ssp245.TF,glaciers(ii).MPIESM12HR.ssp245.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).MPIESM12HR.ssp245.melt(find(glaciers(ii).MPIESM12HR.ssp245.melt<0)) = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(modelid == 7)
disp('Adding MPIESM12HR ssp126')

%% associate future thermal forcing from MPIESM12HR ssp126

% load ocean data
load ../process_CMIP/MPIESM12HR_ssp126.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
cnrmcm6_baseline_inds = find(ismember(year,baseline));
TF0_model = nanmean(squeeze(TF_basins(:,cnrmcm6_baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% bias
bias = TF0_model - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0_model(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).MPIESM12HR.ssp126.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).MPIESM12HR.ssp126.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).MPIESM12HR.ssp126.bias_TF = bias(glaciers(ii).sectornum);
end

%% associate future runoff from MARv3.12/MPIESM12HR ssp126
% load runoff
load ../runoff/MARv3.12-MPI-ESM1-2-HR-ssp126-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_mod = runoff;

% sort time variable
for ii=1:length(runoff_mod),
    runoff_mod(ii).time = [1950:2100];
end

% assign to glaciers
for ii=1:length(glaciers),
    % ssp585
    for jj=1:length(runoff_mod),
        if ismember(glaciers(ii).rignotid,runoff_mod(jj).rignotGlacierID),
            glaciers(ii).MPIESM12HR.ssp126.tJJA = runoff_mod(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).MPIESM12HR.ssp126.tJJA),
                yr = floor(glaciers(ii).MPIESM12HR.ssp126.tJJA(kk));
                glaciers(ii).MPIESM12HR.ssp126.QJJA(kk) = 3*runoff_mod(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).MPIESM12HR.ssp126.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).MPIESM12HR.ssp126.bias_QJJA = mean(glaciers(ii).MPIESM12HR.ssp126.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).MPIESM12HR.ssp126.QJJA = glaciers(ii).MPIESM12HR.ssp126.QJJA - glaciers(ii).MPIESM12HR.ssp126.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).MPIESM12HR.ssp126.QJJA(find(glaciers(ii).MPIESM12HR.ssp126.QJJA<0)) = 0;
end

%% calculate future MPIESM12HR melt ssp126
% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    glaciers(ii).MPIESM12HR.ssp126.tmelt = glaciers(ii).MPIESM12HR.ssp126.tJJA;
    glaciers(ii).MPIESM12HR.ssp126.melt = (glaciers(ii).MPIESM12HR.ssp126.QJJA.^0.4).*...
        interp1(glaciers(ii).MPIESM12HR.ssp126.tTF,glaciers(ii).MPIESM12HR.ssp126.TF,glaciers(ii).MPIESM12HR.ssp126.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).MPIESM12HR.ssp126.melt(find(glaciers(ii).MPIESM12HR.ssp126.melt<0)) = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(modelid == 8)
disp('Adding UKESM1-0-LL-Robin ssp585')

% associate future thermal forcing from UKESM1-0-LL-Robin ssp585

% load UKESM1-0-LL-Robin ssp585 ocean
load ../process_CMIP/UKESM1-0-LL-Robin_ssp585.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
ukesm1_baseline_inds = find(ismember(year,baseline));
TF0_ukesm1 = nanmean(squeeze(TF_basins(:,ukesm1_baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% ukesm1 bias
bias = TF0_ukesm1 - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0_ukesm1(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).UKESM1Robin.ssp585.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).UKESM1Robin.ssp585.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).UKESM1Robin.ssp585.bias_TF = bias(glaciers(ii).sectornum);
end

%% associate future runoff from MAR/UKESM1-0-LL-Robin (ssp 585 only)

% load UKESM1-0-LL-Robin/MAR ssp585
load ../runoff/MARv3.12-UKESM1-0-LL-Robin-ssp585-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_mod = runoff;

% sort time variable
for ii=1:length(runoff_mod),
    runoff_mod(ii).time = [1960:2100]; % 10 years shorter than the others!
end

% assign to glaciers
for ii=1:length(glaciers),
    % ssp585
    for jj=1:length(runoff_mod),
        if ismember(glaciers(ii).rignotid,runoff_mod(jj).rignotGlacierID),
            glaciers(ii).UKESM1Robin.ssp585.tJJA = runoff_mod(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).UKESM1Robin.ssp585.tJJA),
                yr = floor(glaciers(ii).UKESM1Robin.ssp585.tJJA(kk));
                glaciers(ii).UKESM1Robin.ssp585.QJJA(kk) = 3*runoff_mod(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).UKESM1Robin.ssp585.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).UKESM1Robin.ssp585.bias_QJJA = mean(glaciers(ii).UKESM1Robin.ssp585.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).UKESM1Robin.ssp585.QJJA = glaciers(ii).UKESM1Robin.ssp585.QJJA - glaciers(ii).UKESM1Robin.ssp585.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).UKESM1Robin.ssp585.QJJA(find(glaciers(ii).UKESM1Robin.ssp585.QJJA<0)) = 0;
end
%% calculate future UKESM1-0-LL-Robin melt (ssp585 only)

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % ssp585
    glaciers(ii).UKESM1Robin.ssp585.tmelt = glaciers(ii).UKESM1Robin.ssp585.tJJA;
    glaciers(ii).UKESM1Robin.ssp585.melt = (glaciers(ii).UKESM1Robin.ssp585.QJJA.^0.4).*...
        interp1(glaciers(ii).UKESM1Robin.ssp585.tTF,glaciers(ii).UKESM1Robin.ssp585.TF,glaciers(ii).UKESM1Robin.ssp585.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).UKESM1Robin.ssp585.melt(find(glaciers(ii).UKESM1Robin.ssp585.melt<0)) = 0;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(modelid == 9)
disp('Adding NorESM2 ssp585')

% associate future thermal forcing from NorESM2 ssp585

% load NorESM2 ssp585 ocean
load ../process_CMIP/NorESM2_ssp585.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
baseline_inds = find(ismember(year,baseline));
TF0 = nanmean(squeeze(TF_basins(:,baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% bias
bias = TF0 - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).NorESM2.ssp585.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).NorESM2.ssp585.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).NorESM2.ssp585.bias_TF = bias(glaciers(ii).sectornum);
end

%% associate future runoff from MAR/NorESM2 (ssp 585 only)

% load NorESM2/MAR ssp585
load ../runoff/MARv3.12-NorESM2-ssp585-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_mod = runoff;

% sort time variable
for ii=1:length(runoff_mod),
    runoff_mod(ii).time = [1950:2100]; 
end

% assign to glaciers
for ii=1:length(glaciers),
    % ssp585
    for jj=1:length(runoff_mod),
        if ismember(glaciers(ii).rignotid,runoff_mod(jj).rignotGlacierID),
            glaciers(ii).NorESM2.ssp585.tJJA = runoff_mod(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).NorESM2.ssp585.tJJA),
                yr = floor(glaciers(ii).NorESM2.ssp585.tJJA(kk));
                glaciers(ii).NorESM2.ssp585.QJJA(kk) = 3*runoff_mod(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).NorESM2.ssp585.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).NorESM2.ssp585.bias_QJJA = mean(glaciers(ii).NorESM2.ssp585.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).NorESM2.ssp585.QJJA = glaciers(ii).NorESM2.ssp585.QJJA - glaciers(ii).NorESM2.ssp585.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).NorESM2.ssp585.QJJA(find(glaciers(ii).NorESM2.ssp585.QJJA<0)) = 0;
end
%% calculate future NorESM2 melt (ssp585 only)

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % ssp585
    glaciers(ii).NorESM2.ssp585.tmelt = glaciers(ii).NorESM2.ssp585.tJJA;
    glaciers(ii).NorESM2.ssp585.melt = (glaciers(ii).NorESM2.ssp585.QJJA.^0.4).*...
        interp1(glaciers(ii).NorESM2.ssp585.tTF,glaciers(ii).NorESM2.ssp585.TF,glaciers(ii).NorESM2.ssp585.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).NorESM2.ssp585.melt(find(glaciers(ii).NorESM2.ssp585.melt<0)) = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(modelid == 10)
disp('Adding NorESM2 ssp245')

% associate future thermal forcing from NorESM2 ssp245

% load NorESM2 ssp245 ocean
load ../process_CMIP/NorESM2_ssp245.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
baseline_inds = find(ismember(year,baseline));
TF0 = nanmean(squeeze(TF_basins(:,baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% bias
bias = TF0 - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).NorESM2.ssp245.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).NorESM2.ssp245.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).NorESM2.ssp245.bias_TF = bias(glaciers(ii).sectornum);
end

%% associate future runoff from MAR/NorESM2 (ssp 585 only)

% load NorESM2/MAR ssp245
load ../runoff/MARv3.12-NorESM2-ssp245-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_mod = runoff;

% sort time variable
for ii=1:length(runoff_mod),
    runoff_mod(ii).time = [1950:2100]; 
end

% assign to glaciers
for ii=1:length(glaciers),
    % ssp245
    for jj=1:length(runoff_mod),
        if ismember(glaciers(ii).rignotid,runoff_mod(jj).rignotGlacierID),
            glaciers(ii).NorESM2.ssp245.tJJA = runoff_mod(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).NorESM2.ssp245.tJJA),
                yr = floor(glaciers(ii).NorESM2.ssp245.tJJA(kk));
                glaciers(ii).NorESM2.ssp245.QJJA(kk) = 3*runoff_mod(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).NorESM2.ssp245.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).NorESM2.ssp245.bias_QJJA = mean(glaciers(ii).NorESM2.ssp245.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).NorESM2.ssp245.QJJA = glaciers(ii).NorESM2.ssp245.QJJA - glaciers(ii).NorESM2.ssp245.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).NorESM2.ssp245.QJJA(find(glaciers(ii).NorESM2.ssp245.QJJA<0)) = 0;
end
%% calculate future NorESM2 melt 

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % ssp245
    glaciers(ii).NorESM2.ssp245.tmelt = glaciers(ii).NorESM2.ssp245.tJJA;
    glaciers(ii).NorESM2.ssp245.melt = (glaciers(ii).NorESM2.ssp245.QJJA.^0.4).*...
        interp1(glaciers(ii).NorESM2.ssp245.tTF,glaciers(ii).NorESM2.ssp245.TF,glaciers(ii).NorESM2.ssp245.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).NorESM2.ssp245.melt(find(glaciers(ii).NorESM2.ssp245.melt<0)) = 0;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(modelid == 11)
disp('Adding CESM2-CMIP6 ssp585')

% associate future thermal forcing from CESM2-CMIP6 ssp585

% load CESM2-CMIP6 ssp585 ocean
load ../process_CMIP/CESM2-CM6_ssp585.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
baseline_inds = find(ismember(year,baseline));
TF0 = nanmean(squeeze(TF_basins(:,baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% bias
bias = TF0 - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).CESM2CMIP6.ssp585.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).CESM2CMIP6.ssp585.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).CESM2CMIP6.ssp585.bias_TF = bias(glaciers(ii).sectornum);
end

%% associate future runoff from MAR/CESM2-CMIP6 (ssp 585 only)

% load CESM2-CMIP6/MAR ssp585
load ../runoff/MARv3.12-CESM2-CMIP6-ssp585-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_mod = runoff;

% sort time variable
for ii=1:length(runoff_mod),
    runoff_mod(ii).time = [1950:2100]; 
end

% assign to glaciers
for ii=1:length(glaciers),
    % ssp585
    for jj=1:length(runoff_mod),
        if ismember(glaciers(ii).rignotid,runoff_mod(jj).rignotGlacierID),
            glaciers(ii).CESM2CMIP6.ssp585.tJJA = runoff_mod(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).CESM2CMIP6.ssp585.tJJA),
                yr = floor(glaciers(ii).CESM2CMIP6.ssp585.tJJA(kk));
                glaciers(ii).CESM2CMIP6.ssp585.QJJA(kk) = 3*runoff_mod(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).CESM2CMIP6.ssp585.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).CESM2CMIP6.ssp585.bias_QJJA = mean(glaciers(ii).CESM2CMIP6.ssp585.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).CESM2CMIP6.ssp585.QJJA = glaciers(ii).CESM2CMIP6.ssp585.QJJA - glaciers(ii).CESM2CMIP6.ssp585.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).CESM2CMIP6.ssp585.QJJA(find(glaciers(ii).CESM2CMIP6.ssp585.QJJA<0)) = 0;
end
%% calculate future CESM2CMIP6 melt (ssp585 only)

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % ssp585
    glaciers(ii).CESM2CMIP6.ssp585.tmelt = glaciers(ii).CESM2CMIP6.ssp585.tJJA;
    glaciers(ii).CESM2CMIP6.ssp585.melt = (glaciers(ii).CESM2CMIP6.ssp585.QJJA.^0.4).*...
        interp1(glaciers(ii).CESM2CMIP6.ssp585.tTF,glaciers(ii).CESM2CMIP6.ssp585.TF,glaciers(ii).CESM2CMIP6.ssp585.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).CESM2CMIP6.ssp585.melt(find(glaciers(ii).CESM2CMIP6.ssp585.melt<0)) = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(modelid == 12)
disp('Adding CESM2-CMIP6 ssp245')

% associate future thermal forcing from CESM2-CMIP6 ssp245

% load CESM2-CMIP6 ssp245 ocean
load ../process_CMIP/CESM2-CM6_ssp245.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
baseline_inds = find(ismember(year,baseline));
TF0 = nanmean(squeeze(TF_basins(:,baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% bias
bias = TF0 - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).CESM2CMIP6.ssp245.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).CESM2CMIP6.ssp245.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).CESM2CMIP6.ssp245.bias_TF = bias(glaciers(ii).sectornum);
end

%% associate future runoff from MAR/CESM2-CMIP6 (ssp 585 only)

% load CESM2-CMIP6/MAR ssp245
load ../runoff/MARv3.12-CESM2-CMIP6-ssp245-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_mod = runoff;

% sort time variable
for ii=1:length(runoff_mod),
    runoff_mod(ii).time = [1950:2100]; 
end

% assign to glaciers
for ii=1:length(glaciers),
    % ssp245
    for jj=1:length(runoff_mod),
        if ismember(glaciers(ii).rignotid,runoff_mod(jj).rignotGlacierID),
            glaciers(ii).CESM2CMIP6.ssp245.tJJA = runoff_mod(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).CESM2CMIP6.ssp245.tJJA),
                yr = floor(glaciers(ii).CESM2CMIP6.ssp245.tJJA(kk));
                glaciers(ii).CESM2CMIP6.ssp245.QJJA(kk) = 3*runoff_mod(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).CESM2CMIP6.ssp245.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).CESM2CMIP6.ssp245.bias_QJJA = mean(glaciers(ii).CESM2CMIP6.ssp245.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).CESM2CMIP6.ssp245.QJJA = glaciers(ii).CESM2CMIP6.ssp245.QJJA - glaciers(ii).CESM2CMIP6.ssp245.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).CESM2CMIP6.ssp245.QJJA(find(glaciers(ii).CESM2CMIP6.ssp245.QJJA<0)) = 0;
end
%% calculate future CESM2CMIP6 melt 

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % ssp245
    glaciers(ii).CESM2CMIP6.ssp245.tmelt = glaciers(ii).CESM2CMIP6.ssp245.tJJA;
    glaciers(ii).CESM2CMIP6.ssp245.melt = (glaciers(ii).CESM2CMIP6.ssp245.QJJA.^0.4).*...
        interp1(glaciers(ii).CESM2CMIP6.ssp245.tTF,glaciers(ii).CESM2CMIP6.ssp245.TF,glaciers(ii).CESM2CMIP6.ssp245.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).CESM2CMIP6.ssp245.melt(find(glaciers(ii).CESM2CMIP6.ssp245.melt<0)) = 0;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(modelid == 13)
disp('Adding CESM2-CMIP6 ssp126')

% associate future thermal forcing from CESM2-CMIP6 ssp126

% load CESM2-CMIP6 ssp126 ocean
load ../process_CMIP/CESM2-CM6_ssp126.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
baseline_inds = find(ismember(year,baseline));
TF0 = nanmean(squeeze(TF_basins(:,baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% bias
bias = TF0 - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).CESM2CMIP6.ssp126.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).CESM2CMIP6.ssp126.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).CESM2CMIP6.ssp126.bias_TF = bias(glaciers(ii).sectornum);
end

%% associate future runoff from MAR/CESM2-CMIP6 (ssp 585 only)

% load CESM2-CMIP6/MAR ssp126
load ../runoff/MARv3.12-CESM2-CMIP6-ssp126-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_mod = runoff;

% sort time variable
for ii=1:length(runoff_mod),
    runoff_mod(ii).time = [1950:2100]; 
end

% assign to glaciers
for ii=1:length(glaciers),
    % ssp126
    for jj=1:length(runoff_mod),
        if ismember(glaciers(ii).rignotid,runoff_mod(jj).rignotGlacierID),
            glaciers(ii).CESM2CMIP6.ssp126.tJJA = runoff_mod(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).CESM2CMIP6.ssp126.tJJA),
                yr = floor(glaciers(ii).CESM2CMIP6.ssp126.tJJA(kk));
                glaciers(ii).CESM2CMIP6.ssp126.QJJA(kk) = 3*runoff_mod(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).CESM2CMIP6.ssp126.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).CESM2CMIP6.ssp126.bias_QJJA = mean(glaciers(ii).CESM2CMIP6.ssp126.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).CESM2CMIP6.ssp126.QJJA = glaciers(ii).CESM2CMIP6.ssp126.QJJA - glaciers(ii).CESM2CMIP6.ssp126.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).CESM2CMIP6.ssp126.QJJA(find(glaciers(ii).CESM2CMIP6.ssp126.QJJA<0)) = 0;
end
%% calculate future CESM2CMIP6 melt 

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % ssp126
    glaciers(ii).CESM2CMIP6.ssp126.tmelt = glaciers(ii).CESM2CMIP6.ssp126.tJJA;
    glaciers(ii).CESM2CMIP6.ssp126.melt = (glaciers(ii).CESM2CMIP6.ssp126.QJJA.^0.4).*...
        interp1(glaciers(ii).CESM2CMIP6.ssp126.tTF,glaciers(ii).CESM2CMIP6.ssp126.TF,glaciers(ii).CESM2CMIP6.ssp126.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).CESM2CMIP6.ssp126.melt(find(glaciers(ii).CESM2CMIP6.ssp126.melt<0)) = 0;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(modelid == 14)
disp('Adding UKESM1-0-LL-CMIP6 ssp585')

% associate future thermal forcing from UKESM1-0-LL-CMIP6 ssp585

% load UKESM1-0-LL-CMIP6 ssp585 ocean
load ../process_CMIP/UKESM1-0-LL-CM6_ssp585.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
baseline_inds = find(ismember(year,baseline));
TF0 = nanmean(squeeze(TF_basins(:,baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% bias
bias = TF0 - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).UKESM1CMIP6.ssp585.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).UKESM1CMIP6.ssp585.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).UKESM1CMIP6.ssp585.bias_TF = bias(glaciers(ii).sectornum);
end

%% associate future runoff from MAR/UKESM1-0-LL-CMIP6 (ssp 585 only)

% load UKESM1-0-LL-CMIP6/MAR ssp585
load ../runoff/MARv3.12-UKESM1-0-LL-CMIP6-ssp585-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_mod = runoff;

% sort time variable
for ii=1:length(runoff_mod),
    runoff_mod(ii).time = [1950:2100]; 
end

% assign to glaciers
for ii=1:length(glaciers),
    % ssp585
    for jj=1:length(runoff_mod),
        if ismember(glaciers(ii).rignotid,runoff_mod(jj).rignotGlacierID),
            glaciers(ii).UKESM1CMIP6.ssp585.tJJA = runoff_mod(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).UKESM1CMIP6.ssp585.tJJA),
                yr = floor(glaciers(ii).UKESM1CMIP6.ssp585.tJJA(kk));
                glaciers(ii).UKESM1CMIP6.ssp585.QJJA(kk) = 3*runoff_mod(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).UKESM1CMIP6.ssp585.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).UKESM1CMIP6.ssp585.bias_QJJA = mean(glaciers(ii).UKESM1CMIP6.ssp585.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).UKESM1CMIP6.ssp585.QJJA = glaciers(ii).UKESM1CMIP6.ssp585.QJJA - glaciers(ii).UKESM1CMIP6.ssp585.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).UKESM1CMIP6.ssp585.QJJA(find(glaciers(ii).UKESM1CMIP6.ssp585.QJJA<0)) = 0;
end
%% calculate future UKESM1CMIP6 melt (ssp585 only)

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % ssp585
    glaciers(ii).UKESM1CMIP6.ssp585.tmelt = glaciers(ii).UKESM1CMIP6.ssp585.tJJA;
    glaciers(ii).UKESM1CMIP6.ssp585.melt = (glaciers(ii).UKESM1CMIP6.ssp585.QJJA.^0.4).*...
        interp1(glaciers(ii).UKESM1CMIP6.ssp585.tTF,glaciers(ii).UKESM1CMIP6.ssp585.TF,glaciers(ii).UKESM1CMIP6.ssp585.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).UKESM1CMIP6.ssp585.melt(find(glaciers(ii).UKESM1CMIP6.ssp585.melt<0)) = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(modelid == 15)
disp('Adding UKESM1-0-LL-CMIP6 ssp245')

% associate future thermal forcing from UKESM1-0-LL-CMIP6 ssp245

% load UKESM1-0-LL-CMIP6 ssp245 ocean
load ../process_CMIP/UKESM1-0-LL-CM6_ssp245.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
baseline_inds = find(ismember(year,baseline));
TF0 = nanmean(squeeze(TF_basins(:,baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% bias
bias = TF0 - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).UKESM1CMIP6.ssp245.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).UKESM1CMIP6.ssp245.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).UKESM1CMIP6.ssp245.bias_TF = bias(glaciers(ii).sectornum);
end

%% associate future runoff from MAR/UKESM1-0-LL-CMIP6 (ssp 585 only)

% load UKESM1-0-LL-CMIP6/MAR ssp245
load ../runoff/MARv3.12-UKESM1-0-LL-CMIP6-ssp245-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_mod = runoff;

% sort time variable
for ii=1:length(runoff_mod),
    runoff_mod(ii).time = [1950:2100]; 
end

% assign to glaciers
for ii=1:length(glaciers),
    % ssp245
    for jj=1:length(runoff_mod),
        if ismember(glaciers(ii).rignotid,runoff_mod(jj).rignotGlacierID),
            glaciers(ii).UKESM1CMIP6.ssp245.tJJA = runoff_mod(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).UKESM1CMIP6.ssp245.tJJA),
                yr = floor(glaciers(ii).UKESM1CMIP6.ssp245.tJJA(kk));
                glaciers(ii).UKESM1CMIP6.ssp245.QJJA(kk) = 3*runoff_mod(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).UKESM1CMIP6.ssp245.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).UKESM1CMIP6.ssp245.bias_QJJA = mean(glaciers(ii).UKESM1CMIP6.ssp245.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).UKESM1CMIP6.ssp245.QJJA = glaciers(ii).UKESM1CMIP6.ssp245.QJJA - glaciers(ii).UKESM1CMIP6.ssp245.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).UKESM1CMIP6.ssp245.QJJA(find(glaciers(ii).UKESM1CMIP6.ssp245.QJJA<0)) = 0;
end
%% calculate future UKESM1CMIP6 melt 

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % ssp245
    glaciers(ii).UKESM1CMIP6.ssp245.tmelt = glaciers(ii).UKESM1CMIP6.ssp245.tJJA;
    glaciers(ii).UKESM1CMIP6.ssp245.melt = (glaciers(ii).UKESM1CMIP6.ssp245.QJJA.^0.4).*...
        interp1(glaciers(ii).UKESM1CMIP6.ssp245.tTF,glaciers(ii).UKESM1CMIP6.ssp245.TF,glaciers(ii).UKESM1CMIP6.ssp245.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).UKESM1CMIP6.ssp245.melt(find(glaciers(ii).UKESM1CMIP6.ssp245.melt<0)) = 0;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(modelid == 16)
disp('Adding IPSL-CM6A-LR - ssp585')

% associate future thermal forcing from IPSL-CM6A-LR ssp585

% load IPSL-CM6A-LR ssp585 ocean
load ../process_CMIP/IPSL-CM6A-LR_ssp585.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
baseline_inds = find(ismember(year,baseline));
TF0 = nanmean(squeeze(TF_basins(:,baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% bias
bias = TF0 - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).IPSLCM6ALR.ssp585.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).IPSLCM6ALR.ssp585.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).IPSLCM6ALR.ssp585.bias_TF = bias(glaciers(ii).sectornum);
end

%% associate future runoff from MAR/IPSL-CM6A-LR

% load IPSL-CM6A-LR/MAR ssp585
load ../runoff/MARv3.12-IPSL-CM6A-LR-ssp585-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_mod = runoff;

% sort time variable
for ii=1:length(runoff_mod),
    runoff_mod(ii).time = [1950:2100];
end

% assign to glaciers
for ii=1:length(glaciers),
    % ssp585
    for jj=1:length(runoff_mod),
        if ismember(glaciers(ii).rignotid,runoff_mod(jj).rignotGlacierID),
            glaciers(ii).IPSLCM6ALR.ssp585.tJJA = runoff_mod(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).IPSLCM6ALR.ssp585.tJJA),
                yr = floor(glaciers(ii).IPSLCM6ALR.ssp585.tJJA(kk));
                glaciers(ii).IPSLCM6ALR.ssp585.QJJA(kk) = 3*runoff_mod(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).IPSLCM6ALR.ssp585.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).IPSLCM6ALR.ssp585.bias_QJJA = mean(glaciers(ii).IPSLCM6ALR.ssp585.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).IPSLCM6ALR.ssp585.QJJA = glaciers(ii).IPSLCM6ALR.ssp585.QJJA - glaciers(ii).IPSLCM6ALR.ssp585.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).IPSLCM6ALR.ssp585.QJJA(find(glaciers(ii).IPSLCM6ALR.ssp585.QJJA<0)) = 0;
end
%% calculate future IPSLCM6ALR melt

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % ssp585
    glaciers(ii).IPSLCM6ALR.ssp585.tmelt = glaciers(ii).IPSLCM6ALR.ssp585.tJJA;
    glaciers(ii).IPSLCM6ALR.ssp585.melt = (glaciers(ii).IPSLCM6ALR.ssp585.QJJA.^0.4).*...
        interp1(glaciers(ii).IPSLCM6ALR.ssp585.tTF,glaciers(ii).IPSLCM6ALR.ssp585.TF,glaciers(ii).IPSLCM6ALR.ssp585.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).IPSLCM6ALR.ssp585.melt(find(glaciers(ii).IPSLCM6ALR.ssp585.melt<0)) = 0;
end

end

%% save
save glaciers_MARv312.mat glaciers
