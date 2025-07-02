%
% the script was originally made by Vasilii Mikirtumov (https://github.com/vamikirt)
% modified by Xiaofeng Chu (https://github.com/xchu3123)
%
%%

clear all;

cc_std_fact = 2;
excl_r_fact = 4;

rec = dread('path/to/template/binnedTemplate.mrc');
mask_erase = dynamo_ellipsoid((size(rec))/excl_r_fact - 1, size(rec,1), (size(rec,1))/2, 0); %mask_erase = dynamo_ellipsoid(size(rec_resampled)*obj.configuration.exclusion_radius_box_size_ratio, size(rec_resampled,1), size(rec_resampled,1)/2, mask_gaussian_fall_off);


tab_all   = dynamo_table_blank(1);
parfor (i = 1:63, 12)

    tomo_name = char(sprintf('tomo_%03d.TM',i));

    cc_n   = sprintf('%s/cc.mrc',tomo_name);
    tomo_mask = sprintf('%s/tomo_%03d_mask.mrc',tomo_name,i);
    if ~isfile(cc_n)
        continue
    end
    system(sprintf('mkdir %s_peak',tomo_name));
    cc = dread(cc_n);

    if isfile(tomo_mask)
        cc = cc .* dread(tomo_mask);
    end
    tlt = [-60:3:60];
    tilt_max = dread(sprintf('%s/tilt.mrc',tomo_name));
    tdrot_max = dread(sprintf('%s/tdrot.mrc',tomo_name));
    narot_max = dread(sprintf('%s/narot.mrc',tomo_name));

    tab_tomo = dynamo_table_blank(1);

    mean_val = mean(cc(:));
    std_val = std(cc(:));

    half_mask_size = round(size(mask_erase,1)/2);

    disp(sprintf('Extracting peaks from %s',tomo_name));

    L = round(size(mask_erase,1)*0.5); 
    point_file_name = sprintf('%s_peak/%s.point',tomo_name,tomo_name(1:8));
    fid = fopen(point_file_name, "w+");
    precision = 4;
    j = 0;
    cutoff = mean_val + cc_std_fact * std_val;
    while length(tab_tomo(:,1)) < 2000
        j = j + 1;
        [a,b] = dynamo_peak_subpixel(real(cc));
        if b < cutoff || ~any(cc(:)) || j > 1000
            break;
        end
        ar = round([0.5,0.5,0.5]+a);
        if ar(3)+L-1 > 250 || ar(2)+L-1 > 480 || ar(1)+L-1 > 464
            disp(sprintf('skip %d',j))
            cc(ar) = 0;
            continue
        end

        tab_tomo(j,1)  = j;
        tab_tomo(j,24:26) = [0.5,0.5,0.5] + a;
        tab_tomo(j,10) = b;
        tab_tomo(j,7)  = tdrot_max(round(a(1)),round(a(2)), round(a(3)));
        tab_tomo(j,8)  =  tilt_max(round(a(1)),round(a(2)), round(a(3)));
        tab_tomo(j,9)  = narot_max(round(a(1)),round(a(2)), round(a(3)));
        tab_tomo(j,13) = 1;
        tab_tomo(j,14) = tlt(1);
        tab_tomo(j,15) = tlt(end);
        tab_tomo(j,20) = i;
        tab_tomo(j,32) = 1;
        tab_tomo(j,2)  = 1;
        tab_tomo(j,3)  = 1;
      
        cc(ar(1)-L:ar(1)+L-1, ar(2)-L:ar(2)+L-1,ar(3)-L:ar(3)+L-1) =  cc(ar(1)-L : ar(1)+L-1, ar(2)-L : ar(2)+L-1 , ar(3) - L: ar(3)+L-1) .* (1-mask_erase);


        fprintf(fid, "%d %." + precision + "f %." + precision...
            + "f %." + precision...
            + "f \n", 1, round(tab_tomo(j,24),precision),...
            round(tab_tomo(j,25),precision),...
            round(tab_tomo(j,26),precision));

    end

    model_file_name = sprintf('%s_peak/%s.mod',tomo_name,tomo_name(1:8));
    if ~(fseek(fid, 1, 'bof') == -1)
        system("point2model " + point_file_name + " " + model_file_name);
    end
    dwrite(tab_tomo, sprintf("%s_peak/%s.tbl",tomo_name,tomo_name(1:11)));
    fclose(fid); 
    cc_name = sprintf("%s_peak/%s_cc.mrc",tomo_name,tomo_name(1:11));
    dwrite(cc,cc_name)
end


%% 

for i = 1:63
    tomo_name = char(sprintf('tomo_%03d.TM',i));
    if ~isfile(sprintf('%s_peak/%s.tbl',tomo_name,tomo_name))
        continue
    end
    tab_tomo = dread(sprintf('%s_peak/%s.tbl',tomo_name,tomo_name));
    if i ==1 
        tab_all = tab_tomo;
    else
        
        tab_all = cat(1,tab_all,tab_tomo);
    end
end
tab_all(:,1) = [1:length(tab_all)];

%%
vllfile = fopen('bin16.vll', 'w');
for i = 1:63
   fprintf(vllfile,'path/to/tomograms/tomo_%02d.rec\n',i); 
end
fclose(vllfile);
