% script to generate CMIP5 ocean forcings
%clear; close all;

clearvars -except modelid 

% modelid = 1 for NorESM2 esm-piControl aka NorESM2LM_ctrl
% modelid = 2 for NorESM2 esm-bell-1500PgC r2i1p1f1 aka NorESM2LM_B1500r02
% modelid = 3 for NorESM2 esm-bell-2500PgC r2i1p1f1 aka NorESM2LM_B2500r02

%esm-piControl
%esm-bell-1500PgC
%esm-bell-1750PgC
%esm-bell-2000PgC
%esm-bell-2500PgC
%esm-brch-bell-1750PgC-000y-cdr0250PgC
%esm-brch-bell-1750PgC-100y-cdr0250PgC
%esm-brch-bell-2000PgC-000y-cdr0500PgC
%esm-brch-bell-2000PgC-100y-cdr0500PgC
%esm-brch-bell-2500PgC-000y-cdr1000PgC
%esm-brch-bell-2500PgC-100y-cdr1000PgC

modelid = 3;

% freezing point parameters
l1 = -5.73e-2;
l2 = 8.32e-2;
l3 = -7.53e-4;

% regular grid
dd = 50; % grid resolution in km
xreg = [-3e6:dd*10^3:3e6];
yreg = [-4e6:dd*10^3:0e6];
zreg = [200:25:500];
[Xreg,Yreg] = meshgrid(xreg,yreg);


%% NorESM2 esm-N1850
if modelid == 1,

    % path to temperature output
    Tfiletoload = '/projects/NS9708K/heig/NorESM2-LM/esm-N1850/thetao_Oyr_NorESM2-LM_esm-N1850_r2i1p1f1_gr_1880-1909_part.nc';
    % path to salinity output
    Sfiletoload = '/projects/NS9708K/heig/NorESM2-LM/esm-N1850/so_Oyr_NorESM2-LM_esm-N1850_r2i1p1f1_gr_1880-1909_part.nc';
    % problem with lat lon
    latlonname= '/projects/NS9708K/heig/NorESM2-LM/esm-N1850/latlon_part.nc';
    % output filename
    matname ='NorESM2LM_ctrl.mat';
end

%%% IMOPSE
%% NorESM2-LM esm-bell-1500PgC r2i1p1f1
if modelid == 2,

    % path to temperature output
    Tfiletoload = '/projects/NS9708K/heig/NorESM2-LM/esm-bell-1500PgC/thetao_Oyr_NorESM2-LM_esm-bell-1500PgC_r2i1p1f1_gr_1850-2249_part.nc';
    % path to salinity output
    Sfiletoload = '/projects/NS9708K/heig/NorESM2-LM/esm-bell-1500PgC/so_Oyr_NorESM2-LM_esm-bell-1500PgC_r2i1p1f1_gr_1850-2249_part.nc';
    % output filename
    matname ='NorESM2LM_B1500r02.mat';
end

%% NorESM2-LM esm-bell-2500PgC r2i1p1f1
if modelid == 3,

    % path to temperature output
    Tfiletoload = '/projects/NS9708K/heig/NorESM2-LM/esm-bell-2500PgC/thetao_Oyr_NorESM2-LM_esm-bell-2500PgC_r2i1p1f1_gr_1850-2249_part.nc';
    % path to salinity output
    Sfiletoload = '/projects/NS9708K/heig/NorESM2-LM/esm-bell-2500PgC/so_Oyr_NorESM2-LM_esm-bell-2500PgC_r2i1p1f1_gr_1850-2249_part.nc';
    % output filename
    matname ='NorESM2LM_B2500r02.mat';
end

if modelid == 1,
    % read fields from separate latlon file
    lon=ncread(latlonname,'longitude');
    lat=ncread(latlonname,'latitude');
    z=ncread(latlonname,'lev');
    lon_vert=ncread(latlonname,'vertices_longitude');
    lat_vert=ncread(latlonname,'vertices_latitude');

else
    % read fields
    lon=ncread(Tfiletoload,'longitude');
    lat=ncread(Tfiletoload,'latitude');
    z=ncread(Tfiletoload,'lev');
    %time=ncread(Tfiletoload,'time');
    lon_vert=ncread(Tfiletoload,'vertices_longitude');
    lat_vert=ncread(Tfiletoload,'vertices_latitude');
end

% read only TS between 200 and 500 m and one point either side
depthinds = [max(find(z<200)):min(find(z>500))];
z = z(depthinds);
T=ncread(Tfiletoload,'thetao',[1,1,depthinds(1),1],[Inf,Inf,length(depthinds),Inf]);
S=ncread(Sfiletoload,'so',[1,1,depthinds(1),1],[Inf,Inf,length(depthinds),Inf]);

% time field is a bit wierd - make own array
if modelid == 1
year=[1820:1849];
else
year=[1850:2249];
end

% reshape
T = reshape(T,size(T,1)*size(T,2),size(T,3),size(T,4));
S = reshape(S,size(S,1)*size(S,2),size(S,3),size(S,4));
lon = lon(:);
lat = lat(:);

% trim to area of interest
% spatial trim
inds1 = find(lon>270 & lat>55);
inds2 = find(lon<10 & lat>55);
spatialinds = [inds1;inds2];
T = squeeze(T(spatialinds,:,:));
S = squeeze(S(spatialinds,:,:));
lon = lon(spatialinds);
lat = lat(spatialinds);

% check all variables are double
T = double(T);
S = double(S);
lon = double(lon);
lat = double(lat);

% temperature is in kelvin with NaNs set to 0
% change 0s to NaNs
T(find(T==0)) = NaN;
% already in celsius
%T = T-273.15;
% salinity has NaNs set to 0 so change these
S(find(S==0)) = NaN;

% calculate thermal forcing
depth = permute(repmat(z,1,size(S,1),length(year)),[2,1,3]);
TF = T - (l1*S+l2+l3*depth);

% get coords of NorESM points on BedMachine grid
% gridcell vertices
[x_vert,y_vert] = latlon2utm(lat_vert,lon_vert);
% gridcell centres
[x,y] = latlon2utm(lat,lon);

% interpolate TF onto regular grid
for jj=1:size(TF,2),
    for kk=1:size(TF,3),
        TF0 = squeeze(TF(:,jj,kk));
        f = scatteredInterpolant(x,y,TF0,'linear','none');
        TFreg(:,:,jj,kk) = f(Xreg,Yreg);
    end
end

% load ice-ocean basin definitions
load ../final_region_def/ice_ocean_sectors.mat

% find regular grid points inside ice-ocean basins
basin = NaN(size(Xreg,1),size(Xreg,2));
for k=1:7, % since there are 7 basins
    regions(k).inds = find(inpolygon(Xreg,Yreg,regions(k).ocean.x,regions(k).ocean.y));
    basin(regions(k).inds) = k;
end

% get thermal forcing profile per basin per year
TF_basins_z = NaN(7,length(z),length(year));
TFreg_vec = reshape(TFreg,size(TFreg,1)*size(TFreg,2),size(TFreg,3),size(TFreg,4));
for jj=1:7,
    for kk=1:length(z),
        TF_basins_z(jj,kk,:) = nanmean(TFreg_vec(regions(jj).inds,kk,:),1);
    end
end

% take mean of thermal forcing over 200-500 m
% first interpolate onto regular grid then take mean
for jj=1:7,
    for kk=1:length(year),
        TF_basins(jj,kk) = nanmean(interp1(z,TF_basins_z(jj,:,kk),zreg,'linear',NaN));
    end
end

% create rough spatial map of TF by taking naive depth mean
TF = nanmean(TF,2);
TF = squeeze(TF(:,1,:));

% save
save(matname, 'x', 'y', 'TF', 'TF_basins', 'year')


% save historical/reference part of run in slightly different format
if modelid == 1,

%% NorESM2 esm-N1850 aka ctrl
year_historical = year;
TF_historical = TF;
TF_basins_historical = TF_basins;
save(matname, 'x', 'y', 'TF_historical', 'TF_basins_historical', 'year_historical')

end
