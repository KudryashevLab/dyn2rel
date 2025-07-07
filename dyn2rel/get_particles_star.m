% This script generates particles file (.star) in RELION-4.0 format,
% based on provided particles file (.tbl) in Dynamo format
% and corresponding tomograms file (.vll) in Dynamo format.
%
% The script was prepared by:
% Vasilii Mikirtumov, author (https://github.com/vamikirt) 
% Xiaofeng Chu, contributor (https://github.com/xchu3123)
%

%% Activate Dynamo
run /path/to/Dynamo/dynamo_activate.m

%% Activate utils (to use rot_M2eZYZ)
addpath /path/to/dyn2rel/utils

%% Link BatchRunTomo folder from TomoBEAR output
tomobear_BRT_path = '/path/to/tomobear/output/X_BatchRunTomo_1';
imod_files_link_path = './imod_folders';

linking_com = sprintf('ln -s %s %s', tomobear_BRT_path, imod_files_link_path);
system(linking_com);

%% Set input/output paths and parameters
particle_tbl_filepath = '/path/to/particles.tbl';
tomogram_vll_filepath = '/path/to/vll/tomo_list.vll';

% adjust accordingly to your input particles.tbl Dynamo table
bin = 8; 

% path to output particles.star file
particle_star_filepath = '/path/to/particles.star';

%% Load input files
tbl = dread(particle_tbl_filepath);
vll_file = readlines(tomogram_vll_filepath);
tomo_tbl_id_list = sort(unique(tbl(:,20),'stable')); 

%% Generate .star file for particles to be imported in RELION-4.0
TomoManifoldIndex = 1; 
fp_star = fopen(particle_star_filepath, 'w+');
i = 1;

fprintf(fp_star, '\ndata_particles\n\nloop_\n');
fprintf(fp_star, '_rlnTomoName          #%d\n',i); i = i+1;
fprintf(fp_star, '_rlnTomoParticleId    #%d\n',i); i = i+1;
fprintf(fp_star, '_rlnTomoManifoldIndex #%d\n',i); i = i+1;  
fprintf(fp_star, '_rlnCoordinateX       #%d\n',i); i = i+1;  
fprintf(fp_star, '_rlnCoordinateY       #%d\n',i); i = i+1; 
fprintf(fp_star, '_rlnCoordinateZ       #%d\n',i); i = i+1;  
fprintf(fp_star, '_rlnOriginXAngst      #%d\n',i); i = i+1;  
fprintf(fp_star, '_rlnOriginYAngst      #%d\n',i); i = i+1;  
fprintf(fp_star, '_rlnOriginZAngst      #%d\n',i); i = i+1;   
fprintf(fp_star, '_rlnAngleRot          #%d\n',i); i = i+1;  
fprintf(fp_star, '_rlnAngleTilt         #%d\n',i); i = i+1;  
fprintf(fp_star, '_rlnAnglePsi          #%d\n',i); i = i+1;  
fprintf(fp_star, '_rlnClassNumber       #%d\n',i); i = i+1;  
fprintf(fp_star, '_rlnRandomSubset      #%d\n',i); i = i+1; 

for j = 1:length(tomo_tbl_id_list)
    tomo_tbl_id = tomo_tbl_id_list(j) ;
    sel = find(tbl(:,20) == (tomo_tbl_id));
    tbl_tomo = tbl(sel,:);
    if isempty(tbl_tomo)
        continue
        disp(sprintf('tomogram with index %d has no particle', tomo_id));
    end
    vll_tomo_name = vll_file(tomo_tbl_id); 
    vll_tomo_name = regexp(vll_tomo_name, 'tomogram_\d\d\d', 'match');
    vll_tomo_name = vll_tomo_name(1);
    
    for k = 1:size(tbl_tomo,1)
        
        % particle orientation angles - convert to required convention:
        % 1. get rotation matrix from ZXZ angles, used by Dynamo
        R = dynamo_euler2matrix(tbl_tomo(k,7:9));
        % 2. generate from rotation matrix ZYZ angles, used by RELION-4.0
        euler = rot_M2eZYZ(R);
        
        % particle coordinates - recenter and rescale to bin1 (i.e. unbinned)
        pos = bin * ( tbl_tomo(k,[24 25 26]) + tbl_tomo(k,[4 5 6]) );
        
        fprintf(fp_star, '%s\t%06d\t%d',vll_tomo_name,tbl_tomo(k,1),TomoManifoldIndex);
        fprintf(fp_star, '\t%.1f\t%.1f\t%.1f',round(pos(1)),round(pos(2)),round(pos(3)));
        fprintf(fp_star, '\t%.4f\t%.4f\t%.4f',0,0,0); 
        fprintf(fp_star, '\t%12.6f\t%12.6f\t%12.6f',euler(1),euler(2),euler(3));
        fprintf(fp_star, '\t%d\t%12d',1,mod(tbl_tomo(k,1),2)+1); 
        fprintf(fp_star, '\n');
    end
end
fclose(fp_star);