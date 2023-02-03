clear; close all;
% script to add future runoff, ocean TF and melt to glaciers structure
load glaciers_past.mat

% modelid  1 NorESM2LM B1500r02
% modelid  1 NorESM2LM B2500r02

modelid = [1 2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(modelid == 1)
disp('Adding NorESM2LM B1500r02')

%% associate future thermal forcing 

% load ocean
load ../process_CMIP/NorESM2LM_B1500r02.mat
% load baseline
load ../process_CMIP/NorESM2LM_ctrl.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1820:1849];
baseline_inds = find(ismember(year_historical,baseline));
TF0 = nanmean(squeeze(TF_basins_historical(:,baseline_inds)),2);

% en4 baseline
EN4_baseline = [1995:2014];
EN4_baseline_inds = find(ismember(regions(1).t,EN4_baseline));
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
    glaciers(ii).NorESM2LM.B1500r02.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).NorESM2LM.B1500r02.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).NorESM2LM.B1500r02.bias_TF = bias(glaciers(ii).sectornum);
end

%% associate future runoff 

% load runoff
load ../runoff/MEC1-NorESM2LM-B1500r02-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_mod = runoff;

% workaround for constant runoff
for ii=1:length(runoff_mod),
    runoff_mod(ii).time = [1850:2249];
    % fill  
    runoff_mod(ii).runoff(1,1:400)=runoff_mod(ii).runoff(1,1);
end
% sort time variable
%for ii=1:length(runoff_mod),
%    runoff_mod(ii).time = [1850:2249];
%    % duplicate last year 
%    %runoff_mod(ii).runoff(1,151)=runoff_mod(ii).runoff(1,150);
%end

% assign to glaciers
for ii=1:length(glaciers),
    % scen
    for jj=1:length(runoff_mod),
        if ismember(glaciers(ii).rignotid,runoff_mod(jj).rignotGlacierID),
            glaciers(ii).NorESM2LM.B1500r02.tJJA = runoff_mod(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).NorESM2LM.B1500r02.tJJA),
                yr = floor(glaciers(ii).NorESM2LM.B1500r02.tJJA(kk));
                glaciers(ii).NorESM2LM.B1500r02.QJJA(kk) = 3*runoff_mod(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1850:1879];
MARbaselineinds = find(ismember(floor(glaciers(1).NorESM2LM.B1500r02.tJJA),baseline));
RACMObaseline = [1995:2014];
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),RACMObaseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).NorESM2LM.B1500r02.bias_QJJA = mean(glaciers(ii).NorESM2LM.B1500r02.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).NorESM2LM.B1500r02.QJJA = glaciers(ii).NorESM2LM.B1500r02.QJJA - glaciers(ii).NorESM2LM.B1500r02.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).NorESM2LM.B1500r02.QJJA(find(glaciers(ii).NorESM2LM.B1500r02.QJJA<0)) = 0;
end

%% calculate future melt

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % ssp585
    glaciers(ii).NorESM2LM.B1500r02.tmelt = glaciers(ii).NorESM2LM.B1500r02.tJJA;
    glaciers(ii).NorESM2LM.B1500r02.melt = (glaciers(ii).NorESM2LM.B1500r02.QJJA.^0.4).*...
        interp1(glaciers(ii).NorESM2LM.B1500r02.tTF,glaciers(ii).NorESM2LM.B1500r02.TF,glaciers(ii).NorESM2LM.B1500r02.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).NorESM2LM.B1500r02.melt(find(glaciers(ii).NorESM2LM.B1500r02.melt<0)) = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(modelid == 2)
disp('Adding NorESM2LM B2500r02')

%% associate future thermal forcing 

% load ocean
load ../process_CMIP/NorESM2LM_B2500r02.mat
% load baseline
load ../process_CMIP/NorESM2LM_ctrl.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1820:1849];
baseline_inds = find(ismember(year_historical,baseline));
TF0 = nanmean(squeeze(TF_basins_historical(:,baseline_inds)),2);

% en4 baseline
EN4_baseline = [1995:2014];
EN4_baseline_inds = find(ismember(regions(1).t,EN4_baseline));
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
    glaciers(ii).NorESM2LM.B2500r02.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).NorESM2LM.B2500r02.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).NorESM2LM.B2500r02.bias_TF = bias(glaciers(ii).sectornum);
end

%% associate future runoff 

% load runoff
load ../runoff/MEC1-NorESM2LM-B2500r02-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_mod = runoff;

%% workaround for constant runoff
%for ii=1:length(runoff_mod),
%    runoff_mod(ii).time = [1850:2249];
%    % fill  
%    runoff_mod(ii).runoff(1,150:400)=runoff_mod(ii).runoff(1,150);
%end
% sort time variable
%for ii=1:length(runoff_mod),
%    runoff_mod(ii).time = [1850:2249];
%    % duplicate last year 
%    %runoff_mod(ii).runoff(1,151)=runoff_mod(ii).runoff(1,150);
%end

% assign to glaciers
for ii=1:length(glaciers),
    % scen
    for jj=1:length(runoff_mod),
        if ismember(glaciers(ii).rignotid,runoff_mod(jj).rignotGlacierID),
            glaciers(ii).NorESM2LM.B2500r02.tJJA = runoff_mod(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).NorESM2LM.B2500r02.tJJA),
                yr = floor(glaciers(ii).NorESM2LM.B2500r02.tJJA(kk));
                glaciers(ii).NorESM2LM.B2500r02.QJJA(kk) = 3*runoff_mod(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1850:1879];
MARbaselineinds = find(ismember(floor(glaciers(1).NorESM2LM.B2500r02.tJJA),baseline));
RACMObaseline = [1995:2014];
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),RACMObaseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).NorESM2LM.B2500r02.bias_QJJA = mean(glaciers(ii).NorESM2LM.B2500r02.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).NorESM2LM.B2500r02.QJJA = glaciers(ii).NorESM2LM.B2500r02.QJJA - glaciers(ii).NorESM2LM.B2500r02.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).NorESM2LM.B2500r02.QJJA(find(glaciers(ii).NorESM2LM.B2500r02.QJJA<0)) = 0;
end

%% calculate future melt

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % ssp585
    glaciers(ii).NorESM2LM.B2500r02.tmelt = glaciers(ii).NorESM2LM.B2500r02.tJJA;
    glaciers(ii).NorESM2LM.B2500r02.melt = (glaciers(ii).NorESM2LM.B2500r02.QJJA.^0.4).*...
        interp1(glaciers(ii).NorESM2LM.B2500r02.tTF,glaciers(ii).NorESM2LM.B2500r02.TF,glaciers(ii).NorESM2LM.B2500r02.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).NorESM2LM.B2500r02.melt(find(glaciers(ii).NorESM2LM.B2500r02.melt<0)) = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


save glaciers_MEC.mat glaciers
