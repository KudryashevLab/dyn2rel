% This script generates tomograms file (.star) in RELION-4.0 format,
% based on provided particles file (.tbl) in Dynamo format
% and corresponding tomograms file (.vll) in Dynamo format.
%
% The script was prepared by:
% Vasilii Mikirtumov, author (https://github.com/vamikirt) 
% Xiaofeng Chu, contributor (https://github.com/xchu3123)
%

%% Activate Dynamo
run /path/to/Dynamo/dynamo_activate.m

%% Set input/output paths and parameters
particle_tbl_filepath = '/path/to/particles.tbl';
tomogram_vll_filepath = '/path/to/vll/tomo_list.vll';
order_list_csv_filepath = '/path/to/order_list.csv';

dose_per_tilt = 4.75; %in e-/A^2

% path to output tomograms.star file
tomogram_star_filepath = '/path/to/tomograms.star'; 

%% Load input files
tbl = dread(particle_tbl_filepath);
vll_file = readlines(tomogram_vll_filepath);
tomo_tbl_id_list = sort(unique(tbl(:,20),'stable')); 

%% Generate .star file for tomograms to be imported in RELION-4.0
fp_star = fopen(tomogram_star_filepath, 'w+');
i = 1;
fprintf(fp_star, '\ndata_\n\nloop_\n');
fprintf(fp_star, '_rlnTomoName                  #%d\n',i); i = i+1;
fprintf(fp_star, '_rlnTomoTiltSeriesName        #%d\n',i); i = i+1;
fprintf(fp_star, '_rlnTomoImportCtfPlotterFile     #%d\n',i); i = i+1;
fprintf(fp_star, '_rlnTomoImportImodDir         #%d\n',i); i = i+1; 
fprintf(fp_star, '_rlnTomoImportFractionalDose  #%d\n',i); i = i+1; % dose per tilt
fprintf(fp_star, '_rlnTomoImportOrderList       #%d\n',i); i = i+1; % order list is just a text file with 2 columns: acquisition order, tilt angle
for j = 1:length(tomo_tbl_id_list)
    vll_tomo_name = vll_file(j); 
    vll_tomo_name = regexp(vll_tomo_name, 'tomogram_\d\d\d', 'match');
    vll_tomo_name = vll_tomo_name(1);
    fprintf(fp_star, '%s\t%s\t%s\t%s\t%.1f\t%s\n',vll_tomo_name,...
        sprintf('./imod_folders/%s/%s.st',vll_tomo_name, vll_tomo_name),...
        sprintf('./imod_folders/%s/%s.defocus', vll_tomo_name, vll_tomo_name),...
        sprintf('./imod_folders/%s', vll_tomo_name),...
        dose_per_tilt, order_list_csv_filepath);  
end
fclose(fp_star);