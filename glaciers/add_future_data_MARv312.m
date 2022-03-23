clear; close all;
% script to add future runoff, ocean TF and melt to glaciers structure
load glaciers_past.mat

% modelid  1 ACCESS1.3 - rcp85
% modelid  2 CESM2 - ssp585
% modelid  3 CNRM-CM6 - ssp585
% modelid  4 CNRM-ESM2 - ssp585 
% modelid  5 MPIESM12HR - ssp585
% modelid  6 MPIESM12HR - ssp245
% modelid  7 MPIESM12HR - ssp126
% modelid  8 UKESM1-CM6, using UKESM1-0-LL - ssp585]

%modelid = [2];
modelid = [1, 2, 3, 5, 6, 7];

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
disp('Adding CESM2')
%% associate future thermal forcing from CESM2 ssp585

% load CESM2 ssp585 ocean
load ../process_CMIP/CESM2_ssp585.mat
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
    glaciers(ii).CESM2.ssp585.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).CESM2.ssp585.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).CESM2.ssp585.bias_TF = bias(glaciers(ii).sectornum);
end

%% associate future runoff from MAR/CESM2 ssp585

% load CESM2/MAR ssp585
load ../runoff/MARv3.12-CESM2-ssp585-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
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
            glaciers(ii).CESM2.ssp585.tJJA = runoff_mod(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).CESM2.ssp585.tJJA),
                yr = floor(glaciers(ii).CESM2.ssp585.tJJA(kk));
                glaciers(ii).CESM2.ssp585.QJJA(kk) = 3*runoff_mod(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).CESM2.ssp585.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).CESM2.ssp585.bias_QJJA = mean(glaciers(ii).CESM2.ssp585.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).CESM2.ssp585.QJJA = glaciers(ii).CESM2.ssp585.QJJA - glaciers(ii).CESM2.ssp585.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).CESM2.ssp585.QJJA(find(glaciers(ii).CESM2.ssp585.QJJA<0)) = 0;
end

%% calculate future CESM2 melt (ssp585 only)

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % ssp585
    glaciers(ii).CESM2.ssp585.tmelt = glaciers(ii).CESM2.ssp585.tJJA;
    glaciers(ii).CESM2.ssp585.melt = (glaciers(ii).CESM2.ssp585.QJJA.^0.4).*...
        interp1(glaciers(ii).CESM2.ssp585.tTF,glaciers(ii).CESM2.ssp585.TF,glaciers(ii).CESM2.ssp585.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).CESM2.ssp585.melt(find(glaciers(ii).CESM2.ssp585.melt<0)) = 0;
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


%% associate future thermal forcing from UKESM1-0-LL ssp585

%% load UKESM1-0-LL ssp585 ocean
%load ../process_CMIP/UKESM1-0-LL_ssp585.mat
%% load EN4
%load EN4_ISMIP6.mat
%
%% present day baseline period
%baseline = [1995:2014];
%ukesm1_baseline_inds = find(ismember(year,baseline));
%TF0_ukesm1 = nanmean(squeeze(TF_basins(:,ukesm1_baseline_inds)),2);
%
%% en4 baseline
%EN4_baseline_inds = find(ismember(regions(1).t,baseline));
%for l=1:7,
%    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
%end
%
%% ukesm1 bias
%bias = TF0_ukesm1 - TF0_EN4';
%
%% future T for forcing
%clearvars Tfuture
%for l=1:7,
%Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0_ukesm1(l));
%end
%
%% assign to glaciers
%for ii=1:length(glaciers),
%    glaciers(ii).UKESM1.ssp585.tTF = year+0.5; % follow convention of annual means
%    glaciers(ii).UKESM1.ssp585.TF = Tfuture(glaciers(ii).sectornum,:);
%    glaciers(ii).UKESM1.ssp585.bias_TF = bias(glaciers(ii).sectornum);
%end
%
%%% associate future thermal forcing from CNRM-ESM2-1 ssp585
%
%% load CNRM-ESM2-1 ssp585 ocean
%load ../process_CMIP/CNRM-ESM2-1_ssp585.mat
%% load EN4
%load EN4_ISMIP6.mat
%
%% present day baseline period
%baseline = [1995:2014];
%cnrmesm2_baseline_inds = find(ismember(year,baseline));
%TF0_cnrmesm2 = nanmean(squeeze(TF_basins(:,cnrmesm2_baseline_inds)),2);
%
%% en4 baseline
%EN4_baseline_inds = find(ismember(regions(1).t,baseline));
%for l=1:7,
%    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
%end
%
%% cnrmesm2 bias
%bias = TF0_cnrmesm2 - TF0_EN4';
%
%% future T for forcing
%clearvars Tfuture
%for l=1:7,
%Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0_cnrmesm2(l));
%end
%
%% assign to glaciers
%for ii=1:length(glaciers),
%    glaciers(ii).CNRMESM2.ssp585.tTF = year+0.5; % follow convention of annual means
%    glaciers(ii).CNRMESM2.ssp585.TF = Tfuture(glaciers(ii).sectornum,:);
%    glaciers(ii).CNRMESM2.ssp585.bias_TF = bias(glaciers(ii).sectornum);
%end
%%% associate future runoff from MAR/CNRM-ESM2-1 (ssp 585 only)
%
%% load CNRM-ESM2-1/MAR ssp585
%load ../runoff/MARv3.12-CNRM-ESM2-ssp585-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
%runoff_mod = runoff;
%
%% sort time variable
%for ii=1:length(runoff_mod),
%    runoff_mod(ii).time = [1950:2100];
%end
%
%% assign to glaciers
%for ii=1:length(glaciers),
%    % ssp585
%    for jj=1:length(runoff_mod),
%        if ismember(glaciers(ii).rignotid,runoff_mod(jj).rignotGlacierID),
%            glaciers(ii).CNRMESM2.ssp585.tJJA = runoff_mod(jj).time+0.5; % annual mean convention
%            for kk=1:length(glaciers(ii).CNRMESM2.ssp585.tJJA),
%                yr = floor(glaciers(ii).CNRMESM2.ssp585.tJJA(kk));
%                glaciers(ii).CNRMESM2.ssp585.QJJA(kk) = 3*runoff_mod(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
%            end
%        end
%    end
%end
%
%% present day baseline period
%baseline = [1995:2014];
%MARbaselineinds = find(ismember(floor(glaciers(1).CNRMESM2.ssp585.tJJA),baseline));
%RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));
%
%% calculate and account for bias over baseline and add Q baseline
%for ii=1:length(glaciers),
%    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
%    glaciers(ii).CNRMESM2.ssp585.bias_QJJA = mean(glaciers(ii).CNRMESM2.ssp585.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
%    glaciers(ii).CNRMESM2.ssp585.QJJA = glaciers(ii).CNRMESM2.ssp585.QJJA - glaciers(ii).CNRMESM2.ssp585.bias_QJJA;
%    % make sure no runoff values less than 0
%    glaciers(ii).CNRMESM2.ssp585.QJJA(find(glaciers(ii).CNRMESM2.ssp585.QJJA<0)) = 0;
%end
%%% associate future runoff from MAR/UKESM1-0-LL (ssp 585 only)
%
%% load UKESM1-0-LL/MAR ssp585
%load ../runoff/MARv3.12-UKESM1-CM6-ssp585-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
%runoff_mod = runoff;
%
%% sort time variable
%for ii=1:length(runoff_mod),
%    runoff_mod(ii).time = [1950:2100];
%end
%
%% assign to glaciers
%for ii=1:length(glaciers),
%    % ssp585
%    for jj=1:length(runoff_mod),
%        if ismember(glaciers(ii).rignotid,runoff_mod(jj).rignotGlacierID),
%            glaciers(ii).UKESM1.ssp585.tJJA = runoff_mod(jj).time+0.5; % annual mean convention
%            for kk=1:length(glaciers(ii).UKESM1.ssp585.tJJA),
%                yr = floor(glaciers(ii).UKESM1.ssp585.tJJA(kk));
%                glaciers(ii).UKESM1.ssp585.QJJA(kk) = 3*runoff_mod(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
%            end
%        end
%    end
%end
%
%% present day baseline period
%baseline = [1995:2014];
%MARbaselineinds = find(ismember(floor(glaciers(1).UKESM1.ssp585.tJJA),baseline));
%RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));
%
%% calculate and account for bias over baseline and add Q baseline
%for ii=1:length(glaciers),
%    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
%    glaciers(ii).UKESM1.ssp585.bias_QJJA = mean(glaciers(ii).UKESM1.ssp585.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
%    glaciers(ii).UKESM1.ssp585.QJJA = glaciers(ii).UKESM1.ssp585.QJJA - glaciers(ii).UKESM1.ssp585.bias_QJJA;
%    % make sure no runoff values less than 0
%    glaciers(ii).UKESM1.ssp585.QJJA(find(glaciers(ii).UKESM1.ssp585.QJJA<0)) = 0;
%end
%%% calculate future UKESM1-0-LL melt (ssp585 only)
%
%% just have to combine runoff and thermal forcing
%for ii=1:length(glaciers),
%    % ssp585
%    glaciers(ii).UKESM1.ssp585.tmelt = glaciers(ii).UKESM1.ssp585.tJJA;
%    glaciers(ii).UKESM1.ssp585.melt = (glaciers(ii).UKESM1.ssp585.QJJA.^0.4).*...
%        interp1(glaciers(ii).UKESM1.ssp585.tTF,glaciers(ii).UKESM1.ssp585.TF,glaciers(ii).UKESM1.ssp585.tJJA);
%    % make sure no melt values less than 0
%    glaciers(ii).UKESM1.ssp585.melt(find(glaciers(ii).UKESM1.ssp585.melt<0)) = 0;
%end
%%% calculate future CNRM-ESM2-1 melt (ssp585 only)
%
%% just have to combine runoff and thermal forcing
%for ii=1:length(glaciers),
%    % ssp585
%    glaciers(ii).CNRMESM2.ssp585.tmelt = glaciers(ii).CNRMESM2.ssp585.tJJA;
%    glaciers(ii).CNRMESM2.ssp585.melt = (glaciers(ii).CNRMESM2.ssp585.QJJA.^0.4).*...
%        interp1(glaciers(ii).CNRMESM2.ssp585.tTF,glaciers(ii).CNRMESM2.ssp585.TF,glaciers(ii).CNRMESM2.ssp585.tJJA);
%    % make sure no melt values less than 0
%    glaciers(ii).CNRMESM2.ssp585.melt(find(glaciers(ii).CNRMESM2.ssp585.melt<0)) = 0;
%end
%%% save

save glaciers_MARv312.mat glaciers
