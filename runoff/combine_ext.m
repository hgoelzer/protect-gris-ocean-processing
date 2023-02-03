% Recombine with ext
% For runs that are too long to fit into memory, segment processing
% and combine here. 200 years OK, 400 not

clear

% 2 files
infile1='MEC1-NorESM2LM-B2500-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat';
infile2='MEC1-NorESM2LM-B2500ext-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat';
outfile='MEC1-NorESM2LM-B2500r02-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat';

% load
load(infile1)
ext=load(infile2)

% combine
len = length(runoff);
for i=1:len
    runoff(i).runoff(201:400) = ext.runoff(i).runoff(1:200);
    % adjust timer
    runoff(i).time = 1850:2249;
end

save(outfile,'runoff')
