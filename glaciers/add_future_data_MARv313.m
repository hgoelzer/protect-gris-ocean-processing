clear; close all;
% script to add future runoff, ocean TF and melt to glaciers structure
load glaciers_past.mat

% modelid  1 IPSL-CM6A-LR - ssp585 - e04-LWC7_2
% modelid  2 IPSL-CM6A-LR - ssp585 - e54-LWC7_2
% modelid  3 IPSL-CM6A-LR - ssp585 - e05
% modelid  4 IPSL-CM6A-LR - ssp585 - e55
% modelid  5 IPSL-CM6A-LR - ssp585 - e63

%modelid = [1,2,3,4,5];
modelid = [1,3,4,5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MARv3.13 variants

if any(modelid == 1)
disp('Adding IPSL-CM6A-LR - ssp585 - e04')

% associate future thermal forcing

% load ocean
load ../process_CMIP/IPSL-CM6A-LR_ssp585_2300.mat
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
    glaciers(ii).IPSLCM6ALRe04.ssp585.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).IPSLCM6ALRe04.ssp585.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).IPSLCM6ALRe04.ssp585.bias_TF = bias(glaciers(ii).sectornum);
end

%% associate future runoff

% load runoff
load ../runoff/MARv3.13-e04-LWC7_2-IPSL-CM6A-LR-ssp585-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_mod = runoff;

% sort time variable
for ii=1:length(runoff_mod),
    runoff_mod(ii).time = [1950:2300];
end

% assign to glaciers
for ii=1:length(glaciers),
    % ssp585
    for jj=1:length(runoff_mod),
        if ismember(glaciers(ii).rignotid,runoff_mod(jj).rignotGlacierID),
            glaciers(ii).IPSLCM6ALRe04.ssp585.tJJA = runoff_mod(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).IPSLCM6ALRe04.ssp585.tJJA),
                yr = floor(glaciers(ii).IPSLCM6ALRe04.ssp585.tJJA(kk));
                glaciers(ii).IPSLCM6ALRe04.ssp585.QJJA(kk) = 3*runoff_mod(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).IPSLCM6ALRe04.ssp585.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).IPSLCM6ALRe04.ssp585.bias_QJJA = mean(glaciers(ii).IPSLCM6ALRe04.ssp585.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).IPSLCM6ALRe04.ssp585.QJJA = glaciers(ii).IPSLCM6ALRe04.ssp585.QJJA - glaciers(ii).IPSLCM6ALRe04.ssp585.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).IPSLCM6ALRe04.ssp585.QJJA(find(glaciers(ii).IPSLCM6ALRe04.ssp585.QJJA<0)) = 0;
end
%% calculate future melt

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % ssp585
    glaciers(ii).IPSLCM6ALRe04.ssp585.tmelt = glaciers(ii).IPSLCM6ALRe04.ssp585.tJJA;
    glaciers(ii).IPSLCM6ALRe04.ssp585.melt = (glaciers(ii).IPSLCM6ALRe04.ssp585.QJJA.^0.4).*...
        interp1(glaciers(ii).IPSLCM6ALRe04.ssp585.tTF,glaciers(ii).IPSLCM6ALRe04.ssp585.TF,glaciers(ii).IPSLCM6ALRe04.ssp585.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).IPSLCM6ALRe04.ssp585.melt(find(glaciers(ii).IPSLCM6ALRe04.ssp585.melt<0)) = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(modelid == 2)
disp('Adding IPSL-CM6A-LR - ssp585 - e54')

% associate future thermal forcing

% load ocean
load ../process_CMIP/IPSL-CM6A-LR_ssp585_2300.mat
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
    glaciers(ii).IPSLCM6ALRe54.ssp585.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).IPSLCM6ALRe54.ssp585.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).IPSLCM6ALRe54.ssp585.bias_TF = bias(glaciers(ii).sectornum);
end

%% associate future runoff

% load runoff
load ../runoff/MARv3.13-e54-LWC7_2-IPSL-CM6A-LR-ssp585-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_mod = runoff;

% sort time variable
for ii=1:length(runoff_mod),
    runoff_mod(ii).time = [1950:2300];
end

% assign to glaciers
for ii=1:length(glaciers),
    % ssp585
    for jj=1:length(runoff_mod),
        if ismember(glaciers(ii).rignotid,runoff_mod(jj).rignotGlacierID),
            glaciers(ii).IPSLCM6ALRe54.ssp585.tJJA = runoff_mod(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).IPSLCM6ALRe54.ssp585.tJJA),
                yr = floor(glaciers(ii).IPSLCM6ALRe54.ssp585.tJJA(kk));
                glaciers(ii).IPSLCM6ALRe54.ssp585.QJJA(kk) = 3*runoff_mod(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).IPSLCM6ALRe54.ssp585.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).IPSLCM6ALRe54.ssp585.bias_QJJA = mean(glaciers(ii).IPSLCM6ALRe54.ssp585.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).IPSLCM6ALRe54.ssp585.QJJA = glaciers(ii).IPSLCM6ALRe54.ssp585.QJJA - glaciers(ii).IPSLCM6ALRe54.ssp585.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).IPSLCM6ALRe54.ssp585.QJJA(find(glaciers(ii).IPSLCM6ALRe54.ssp585.QJJA<0)) = 0;
end
%% calculate future 

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % ssp585
    glaciers(ii).IPSLCM6ALRe54.ssp585.tmelt = glaciers(ii).IPSLCM6ALRe54.ssp585.tJJA;
    glaciers(ii).IPSLCM6ALRe54.ssp585.melt = (glaciers(ii).IPSLCM6ALRe54.ssp585.QJJA.^0.4).*...
        interp1(glaciers(ii).IPSLCM6ALRe54.ssp585.tTF,glaciers(ii).IPSLCM6ALRe54.ssp585.TF,glaciers(ii).IPSLCM6ALRe54.ssp585.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).IPSLCM6ALRe54.ssp585.melt(find(glaciers(ii).IPSLCM6ALRe54.ssp585.melt<0)) = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(modelid == 3)
disp('Adding IPSL-CM6A-LR - ssp585 - e05')

% associate future thermal forcing

% load ocean
load ../process_CMIP/IPSL-CM6A-LR_ssp585_2300.mat
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
    glaciers(ii).IPSLCM6ALRe05.ssp585.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).IPSLCM6ALRe05.ssp585.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).IPSLCM6ALRe05.ssp585.bias_TF = bias(glaciers(ii).sectornum);
end

%% associate future runoff

% load runoff
load ../runoff/MARv3.13-e05-IPSL-CM6A-LR-ssp585-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_mod = runoff;

% sort time variable
for ii=1:length(runoff_mod),
    runoff_mod(ii).time = [1950:2300];
end

% assign to glaciers
for ii=1:length(glaciers),
    % ssp585
    for jj=1:length(runoff_mod),
        if ismember(glaciers(ii).rignotid,runoff_mod(jj).rignotGlacierID),
            glaciers(ii).IPSLCM6ALRe05.ssp585.tJJA = runoff_mod(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).IPSLCM6ALRe05.ssp585.tJJA),
                yr = floor(glaciers(ii).IPSLCM6ALRe05.ssp585.tJJA(kk));
                glaciers(ii).IPSLCM6ALRe05.ssp585.QJJA(kk) = 3*runoff_mod(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).IPSLCM6ALRe05.ssp585.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).IPSLCM6ALRe05.ssp585.bias_QJJA = mean(glaciers(ii).IPSLCM6ALRe05.ssp585.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).IPSLCM6ALRe05.ssp585.QJJA = glaciers(ii).IPSLCM6ALRe05.ssp585.QJJA - glaciers(ii).IPSLCM6ALRe05.ssp585.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).IPSLCM6ALRe05.ssp585.QJJA(find(glaciers(ii).IPSLCM6ALRe05.ssp585.QJJA<0)) = 0;
end
%% calculate future melt

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % ssp585
    glaciers(ii).IPSLCM6ALRe05.ssp585.tmelt = glaciers(ii).IPSLCM6ALRe05.ssp585.tJJA;
    glaciers(ii).IPSLCM6ALRe05.ssp585.melt = (glaciers(ii).IPSLCM6ALRe05.ssp585.QJJA.^0.4).*...
        interp1(glaciers(ii).IPSLCM6ALRe05.ssp585.tTF,glaciers(ii).IPSLCM6ALRe05.ssp585.TF,glaciers(ii).IPSLCM6ALRe05.ssp585.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).IPSLCM6ALRe05.ssp585.melt(find(glaciers(ii).IPSLCM6ALRe05.ssp585.melt<0)) = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(modelid == 4)
disp('Adding IPSL-CM6A-LR - ssp585 - e55')

% associate future thermal forcing

% load ocean
load ../process_CMIP/IPSL-CM6A-LR_ssp585_2300.mat
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
    glaciers(ii).IPSLCM6ALRe55.ssp585.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).IPSLCM6ALRe55.ssp585.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).IPSLCM6ALRe55.ssp585.bias_TF = bias(glaciers(ii).sectornum);
end

%% associate future runoff

% load runoff
load ../runoff/MARv3.13-e55-IPSL-CM6A-LR-ssp585-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_mod = runoff;

% sort time variable
for ii=1:length(runoff_mod),
    runoff_mod(ii).time = [1950:2300];
end

% assign to glaciers
for ii=1:length(glaciers),
    % ssp585
    for jj=1:length(runoff_mod),
        if ismember(glaciers(ii).rignotid,runoff_mod(jj).rignotGlacierID),
            glaciers(ii).IPSLCM6ALRe55.ssp585.tJJA = runoff_mod(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).IPSLCM6ALRe55.ssp585.tJJA),
                yr = floor(glaciers(ii).IPSLCM6ALRe55.ssp585.tJJA(kk));
                glaciers(ii).IPSLCM6ALRe55.ssp585.QJJA(kk) = 3*runoff_mod(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).IPSLCM6ALRe55.ssp585.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).IPSLCM6ALRe55.ssp585.bias_QJJA = mean(glaciers(ii).IPSLCM6ALRe55.ssp585.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).IPSLCM6ALRe55.ssp585.QJJA = glaciers(ii).IPSLCM6ALRe55.ssp585.QJJA - glaciers(ii).IPSLCM6ALRe55.ssp585.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).IPSLCM6ALRe55.ssp585.QJJA(find(glaciers(ii).IPSLCM6ALRe55.ssp585.QJJA<0)) = 0;
end
%% calculate future melt

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % ssp585
    glaciers(ii).IPSLCM6ALRe55.ssp585.tmelt = glaciers(ii).IPSLCM6ALRe55.ssp585.tJJA;
    glaciers(ii).IPSLCM6ALRe55.ssp585.melt = (glaciers(ii).IPSLCM6ALRe55.ssp585.QJJA.^0.4).*...
        interp1(glaciers(ii).IPSLCM6ALRe55.ssp585.tTF,glaciers(ii).IPSLCM6ALRe55.ssp585.TF,glaciers(ii).IPSLCM6ALRe55.ssp585.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).IPSLCM6ALRe55.ssp585.melt(find(glaciers(ii).IPSLCM6ALRe55.ssp585.melt<0)) = 0;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(modelid == 5)
disp('Adding IPSL-CM6A-LR - ssp585 - e63')

% associate future thermal forcing

% load ocean
load ../process_CMIP/IPSL-CM6A-LR_ssp585_2300.mat
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
    glaciers(ii).IPSLCM6ALRe63.ssp585.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).IPSLCM6ALRe63.ssp585.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).IPSLCM6ALRe63.ssp585.bias_TF = bias(glaciers(ii).sectornum);
end

%% associate future runoff

% load runoff
load ../runoff/MARv3.13-e63-IPSL-CM6A-LR-ssp585-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_mod = runoff;

% sort time variable
for ii=1:length(runoff_mod),
    runoff_mod(ii).time = [1950:2300];
end

% assign to glaciers
for ii=1:length(glaciers),
    % ssp585
    for jj=1:length(runoff_mod),
        if ismember(glaciers(ii).rignotid,runoff_mod(jj).rignotGlacierID),
            glaciers(ii).IPSLCM6ALRe63.ssp585.tJJA = runoff_mod(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).IPSLCM6ALRe63.ssp585.tJJA),
                yr = floor(glaciers(ii).IPSLCM6ALRe63.ssp585.tJJA(kk));
                glaciers(ii).IPSLCM6ALRe63.ssp585.QJJA(kk) = 3*runoff_mod(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).IPSLCM6ALRe63.ssp585.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).IPSLCM6ALRe63.ssp585.bias_QJJA = mean(glaciers(ii).IPSLCM6ALRe63.ssp585.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).IPSLCM6ALRe63.ssp585.QJJA = glaciers(ii).IPSLCM6ALRe63.ssp585.QJJA - glaciers(ii).IPSLCM6ALRe63.ssp585.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).IPSLCM6ALRe63.ssp585.QJJA(find(glaciers(ii).IPSLCM6ALRe63.ssp585.QJJA<0)) = 0;
end
%% calculate future melt

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % ssp585
    glaciers(ii).IPSLCM6ALRe63.ssp585.tmelt = glaciers(ii).IPSLCM6ALRe63.ssp585.tJJA;
    glaciers(ii).IPSLCM6ALRe63.ssp585.melt = (glaciers(ii).IPSLCM6ALRe63.ssp585.QJJA.^0.4).*...
        interp1(glaciers(ii).IPSLCM6ALRe63.ssp585.tTF,glaciers(ii).IPSLCM6ALRe63.ssp585.TF,glaciers(ii).IPSLCM6ALRe63.ssp585.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).IPSLCM6ALRe63.ssp585.melt(find(glaciers(ii).IPSLCM6ALRe63.ssp585.melt<0)) = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% save
save glaciers_MARv313.mat glaciers
