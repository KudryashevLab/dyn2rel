% 
% the script takes .tbl file from dynamo and corresponding .vll file
% generate tomo.star and particle.star for relion ver4.0 
% the script was originally made by Vasilii Mikirtumov (https://github.com/vamikirt)
% modified by Xiaofeng Chu (https://github.com/xchu3123)
%
%% activate auxiliarly utils
addpath /path/to/dyn2rel/utils

%% preparation of ptcls.star
bear_path = 'path/to/tomobear/output'; 
tbl = dread('path/to/particle/tbl//class1.tbl');
vll_file = readlines('path/to/vll/tomo_list.vll');
tomo_tbl_id_list = sort(unique(tbl(:,20),'stable')); 

linking_com = sprintf('ln -s %s/11_BatchRunTomo_1 ./imod_folders', bear_path);
system(linking_com);

unbin_apix = 0.84;
bin = 8; 
TomoManifoldIndex = 1; 
fp_star = fopen('./particles.star', 'w+');
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
        disp(sprintf('tomogram with index %d has no particle', tomo_id))
    end
    vll_tomo_name = vll_file(tomo_tbl_id); 
    vll_tomo_name = regexp(vll_tomo_name, 'tomogram_\d\d\d', 'match');
    vll_tomo_name = vll_tomo_name(1);
    
    for k = 1:size(tbl_tomo,1) 
        R = dynamo_euler2matrix(tbl_tomo(k,7:9)); 
        euler = rot_M2eZYZ(R); 
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

%% preparin tomos.star

dose_per_tilt = 4.75;
fp_star = fopen('tomo.star', 'w+');
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
        dose_per_tilt, 'order_list.csv');  
end
fclose(fp_star);
