% Analyse relative importance of runoff and tf

df = load('runoff_tf_MARv312.mat');
rtf = df.retreat;

load('projected_retreat_MARv312.mat');

dataset = zeros(7,16,3);

t0 = 65;
tend = 151;

for l = 1:7;

% ACCESS RCP8.5
dataset(l,1,1) = rtf.ACCESS.RCP85.Qmed(l,tend)-rtf.ACCESS.RCP85.Qmed(l,t0);
dataset(l,1,2) = rtf.ACCESS.RCP85.TF(l,tend)-rtf.ACCESS.RCP85.TF(l,t0);
dataset(l,1,3) = retreat.ACCESS.RCP85.med(l,tend)-retreat.ACCESS.RCP85.med(l,t0);
% CESM2 ssp585
dataset(l,2,1) = rtf.CESM2Leo.ssp585.Qmed(l,tend)-rtf.CESM2Leo.ssp585.Qmed(l,t0);
dataset(l,2,2) = rtf.CESM2Leo.ssp585.TF(l,tend)-rtf.CESM2Leo.ssp585.TF(l,t0);
dataset(l,2,3) = retreat.CESM2Leo.ssp585.med(l,tend)-retreat.CESM2Leo.ssp585.med(l,t0);
% CNRM-CM6-1 ssp585
dataset(l,3,1) = rtf.CNRMCM6.ssp585.Qmed(l,tend)-rtf.CNRMCM6.ssp585.Qmed(l,t0);
dataset(l,3,2) = rtf.CNRMCM6.ssp585.TF(l,tend)-rtf.CNRMCM6.ssp585.TF(l,t0);
dataset(l,3,3) = retreat.CNRMCM6.ssp585.med(l,tend)-retreat.CNRMCM6.ssp585.med(l,t0);
% CNRM-ESM2-1 ssp585
dataset(l,4,1) = rtf.CNRMESM2.ssp585.Qmed(l,tend)-rtf.CNRMESM2.ssp585.Qmed(l,t0);
dataset(l,4,2) = rtf.CNRMESM2.ssp585.TF(l,tend)-rtf.CNRMESM2.ssp585.TF(l,t0);
dataset(l,4,3) = retreat.CNRMESM2.ssp585.med(l,tend)-retreat.CNRMESM2.ssp585.med(l,t0);
% MPIESM12HR 
dataset(l,5,1) = rtf.MPIESM12HR.ssp585.Qmed(l,tend)-rtf.MPIESM12HR.ssp585.Qmed(l,t0);
dataset(l,5,2) = rtf.MPIESM12HR.ssp585.TF(l,tend)-rtf.MPIESM12HR.ssp585.TF(l,t0);
dataset(l,5,3) = retreat.MPIESM12HR.ssp585.med(l,tend)-retreat.MPIESM12HR.ssp585.med(l,t0);
dataset(l,6,1) = rtf.MPIESM12HR.ssp245.Qmed(l,tend)-rtf.MPIESM12HR.ssp245.Qmed(l,t0);
dataset(l,6,2) = rtf.MPIESM12HR.ssp245.TF(l,tend)-rtf.MPIESM12HR.ssp245.TF(l,t0);
dataset(l,6,3) = retreat.MPIESM12HR.ssp245.med(l,tend)-retreat.MPIESM12HR.ssp245.med(l,t0);
dataset(l,7,1) = rtf.MPIESM12HR.ssp126.Qmed(l,tend)-rtf.MPIESM12HR.ssp126.Qmed(l,t0);
dataset(l,7,2) = rtf.MPIESM12HR.ssp126.TF(l,tend)-rtf.MPIESM12HR.ssp126.TF(l,t0);
dataset(l,7,3) = retreat.MPIESM12HR.ssp126.med(l,tend)-retreat.MPIESM12HR.ssp126.med(l,t0);
% UKESM1-0-LL ssp585
dataset(l,8,1) = rtf.UKESM1Robin.ssp585.Qmed(l,tend)-rtf.UKESM1Robin.ssp585.Qmed(l,t0);
dataset(l,8,2) = rtf.UKESM1Robin.ssp585.TF(l,tend)-rtf.UKESM1Robin.ssp585.TF(l,t0);
dataset(l,8,3) = retreat.UKESM1Robin.ssp585.med(l,tend)-retreat.UKESM1Robin.ssp585.med(l,t0);
% NorESM2
dataset(l,9,1) = rtf.NorESM2.ssp585.Qmed(l,tend)-rtf.NorESM2.ssp585.Qmed(l,t0);
dataset(l,9,2) = rtf.NorESM2.ssp585.TF(l,tend)-rtf.NorESM2.ssp585.TF(l,t0);
dataset(l,9,3) = retreat.NorESM2.ssp585.med(l,tend)-retreat.NorESM2.ssp585.med(l,t0);
dataset(l,10,1) = rtf.NorESM2.ssp245.Qmed(l,tend)-rtf.NorESM2.ssp245.Qmed(l,t0);
dataset(l,10,2) = rtf.NorESM2.ssp245.TF(l,tend)-rtf.NorESM2.ssp245.TF(l,t0);
dataset(l,10,3) = retreat.NorESM2.ssp245.med(l,tend)-retreat.NorESM2.ssp245.med(l,t0);
% CESM2CMIP6
dataset(l,11,1) = rtf.CESM2CMIP6.ssp585.Qmed(l,tend)-rtf.CESM2CMIP6.ssp585.Qmed(l,t0);
dataset(l,11,2) = rtf.CESM2CMIP6.ssp585.TF(l,tend)-rtf.CESM2CMIP6.ssp585.TF(l,t0);
dataset(l,11,3) = retreat.CESM2CMIP6.ssp585.med(l,tend)-retreat.CESM2CMIP6.ssp585.med(l,t0);
dataset(l,12,1) = rtf.CESM2CMIP6.ssp245.Qmed(l,tend)-rtf.CESM2CMIP6.ssp245.Qmed(l,t0);
dataset(l,12,2) = rtf.CESM2CMIP6.ssp245.TF(l,tend)-rtf.CESM2CMIP6.ssp245.TF(l,t0);
dataset(l,12,3) = retreat.CESM2CMIP6.ssp245.med(l,tend)-retreat.CESM2CMIP6.ssp245.med(l,t0);
dataset(l,13,1) = rtf.CESM2CMIP6.ssp126.Qmed(l,tend)-rtf.CESM2CMIP6.ssp126.Qmed(l,t0);
dataset(l,13,2) = rtf.CESM2CMIP6.ssp126.TF(l,tend)-rtf.CESM2CMIP6.ssp126.TF(l,t0);
dataset(l,13,3) = retreat.CESM2CMIP6.ssp126.med(l,tend)-retreat.CESM2CMIP6.ssp126.med(l,t0);
% UKESM1CMIP6
dataset(l,14,1) = rtf.UKESM1CMIP6.ssp585.Qmed(l,tend)-rtf.UKESM1CMIP6.ssp585.Qmed(l,t0);
dataset(l,14,2) = rtf.UKESM1CMIP6.ssp585.TF(l,tend)-rtf.UKESM1CMIP6.ssp585.TF(l,t0);
dataset(l,14,3) = retreat.UKESM1CMIP6.ssp585.med(l,tend)-retreat.UKESM1CMIP6.ssp585.med(l,t0);
dataset(l,15,1) = rtf.UKESM1CMIP6.ssp245.Qmed(l,tend)-rtf.UKESM1CMIP6.ssp245.Qmed(l,t0);
dataset(l,15,2) = rtf.UKESM1CMIP6.ssp245.TF(l,tend)-rtf.UKESM1CMIP6.ssp245.TF(l,t0);
dataset(l,15,3) = retreat.UKESM1CMIP6.ssp245.med(l,tend)-retreat.UKESM1CMIP6.ssp245.med(l,t0);
% IPSLCM6ALR ssp585
dataset(l,16,1) = rtf.IPSLCM6ALR.ssp585.Qmed(l,tend)-rtf.IPSLCM6ALR.ssp585.Qmed(l,t0);
dataset(l,16,2) = rtf.IPSLCM6ALR.ssp585.TF(l,tend)-rtf.IPSLCM6ALR.ssp585.TF(l,t0);
dataset(l,16,3) = retreat.IPSLCM6ALR.ssp585.med(l,tend)-retreat.IPSLCM6ALR.ssp585.med(l,t0);

end

% sign change for better interpretation
for l = 1:7
    for m = 1:16
        dataset(l,m,1) = -dataset(l,m,1);
        dataset(l,m,3) = -dataset(l,m,3);
    end
end

% model selection
% mm =  1 for ACCESS RCP8.5
% mm =  2 for CESM2Leo ssp585
% mm =  3 for CNRM-CM6-1 ssp585
% mm =  4 for CNRM-ESM2-1 ssp585
% mm =  5 for MPI-ESM1-2-H ssp585
% mm =  6 for MPI-ESM1-2-H ssp245
% mm =  7 for MPI-ESM1-2-H ssp126
% mm =  8 for UKESM1Robin ssp585
% mm =  9 for NorESM2 - ssp585
% mm = 10 for NorESM2 - ssp245
% mm = 11 for CESM2-CMIP6 - ssp585
% mm = 12 for CESM2-CMIP6 - ssp245
% mm = 13 for CESM2-CMIP6 - ssp126
% mm = 14 for UKESM1-0-LL-CMIP6 - ssp585
% mm = 15 for UKESM1-0-LL-CMIP6 - ssp245
% mm = 16 for IPSL-CM6A-LR - ssp585

%size(dataset)
%% only keep ssp585 ; start deleting from end
%%for m = [15, 13, 12, 10, 7, 6]
%% only keep ssp585 and remove doubles; start deleting from end
%for m = [15, 13, 12, 10, 8, 7, 6, 2]
%%for m = [15, 13, 12, 11, 10, 8, 7, 6]
%% remove ssp585; start deleting from end
%%for m = [16, 14, 11, 9, 8, 5, 4, 3, 2, 1]
%    dataset(:,m,:) = [];
%end
%size(dataset)

% only keep some regions 1 - 7; start deleting from end
% 'SE' 'SW' 'CE' 'CW' 'NE' 'NW' 'NO'
size(dataset)
%for l = [7, 5, 4, 3, 2, 1]
for l = [7, 5, 3, 2, 1]
    dataset(l,:,:) = [];
end
size(dataset)

nl = size(dataset,1)

% correlation
cc = corrcoef(squeeze(dataset(:,:,1)),squeeze(dataset(:,:,3))); cc1 = cc(1,2)
cc = corrcoef(squeeze(dataset(:,:,2)),squeeze(dataset(:,:,3))); cc2 = cc(1,2)
cc = corrcoef(squeeze(dataset(:,:,1)),squeeze(dataset(:,:,2))); cc12 = cc(1,2)

figure
hold on; box on
for l=1:nl
    plot(dataset(l,:,1),dataset(l,:,3),'o')
end
xlabel('Runoff')
ylabel('Retreat')
title(num2str(cc1))

figure
hold on; box on
for l=1:nl
    plot(dataset(l,:,2),dataset(l,:,3),'o')
end
xlabel('TF')
ylabel('Retreat')
title(num2str(cc2))

figure
hold on; box on
for l=1:nl
    plot(dataset(l,:,1),dataset(l,:,2),'o')
end
xlabel('Runoff')
ylabel('TF')
title(num2str(cc12))
