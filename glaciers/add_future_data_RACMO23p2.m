clear; close all;
% script to add future runoff, ocean TF and melt to glaciers structure
load glaciers_past.mat

% modelid  1 CESM2-Leo - ssp585
% modelid  2 CESM2-CMIP6 - ssp245
% modelid  3 CESM2-CMIP6 - ssp126

%modelid = [1];
%modelid = [2];
modelid = [1, 2, 3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(modelid == 1)
disp('Adding CESM2-Leo - ssp585')

%% associate future thermal forcing from CESM2-Leo ssp585

% load CESM2 ssp585 ocean
load ../process_CMIP/CESM2-Leo_ssp585.mat
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
    glaciers(ii).CESM2Leo.ssp585.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).CESM2Leo.ssp585.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).CESM2Leo.ssp585.bias_TF = bias(glaciers(ii).sectornum);
end

%% associate future runoff from CESM2-Leo ssp585

% load CESM2-Leo/RACMO ssp585
load ../runoff/RACMO2.3p2-CESM2-Leo-ssp585-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_mod = runoff;


% sort time variable
for ii=1:length(runoff_mod),
    runoff_mod(ii).time = [1950:2100];
    % duplicate last year for CESM2-Leo which misses year 2100 
    runoff_mod(ii).runoff(1,151)=runoff_mod(ii).runoff(1,150);
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

%% calculate future CESM2-Leo melt

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

if any(modelid == 2)
disp('Adding CESM2CMIP6 - ssp245')

%% associate future thermal forcing from CESM2-CMIP6

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

%% associate future runoff from CESM2-CMIP6/RACMO ssp245

% load CESM2-CMIP6/RACMO ssp245
load ../runoff/RACMO2.3p2-CESM2-CMIP6-ssp245-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_mod = runoff;


% sort time variable
for ii=1:length(runoff_mod),
    runoff_mod(ii).time = [1950:2100];
    % duplicate last year for CESM2_CMIP6 which misses year 2100 
    runoff_mod(ii).runoff(1,151)=runoff_mod(ii).runoff(1,150);
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

%% calculate future CESM2-CMIP6 melt

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

if any(modelid == 3)
disp('Adding CESM2CMIP6 - ssp126')

%% associate future thermal forcing from CESM2-CMIP6 ssp126

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

%% associate future runoff from CESM2-CMIP6/RACMO ssp126

% load CESM2-CMIP6/RACMO ssp126
load ../runoff/RACMO2.3p2-CESM2-CMIP6-ssp126-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_mod = runoff;


% sort time variable
for ii=1:length(runoff_mod),
    runoff_mod(ii).time = [1950:2100];
    % duplicate last year for CESM which misses year 2100 
    runoff_mod(ii).runoff(1,151)=runoff_mod(ii).runoff(1,150);
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

%% calculate future CESM2-CMIP6 melt

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

save glaciers_RACMO23p2.mat glaciers
