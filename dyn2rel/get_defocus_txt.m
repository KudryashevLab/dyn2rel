% This script generates defocus (.txt) file in IMOD format,
% based on the output log files (_gctf.log) from GCTF,
% with the following tab-separated parameters:
% - micrograph number
% - defocus u
% - defocus v
% - astigmatism angle
% - phase shift
% - ( = 0)
% - cross-correlation
% - resolution limit
%
% The script was prepared by:
% Vasilii Mikirtumov (https://github.com/vamikirt)

%% Set input/output paths
tomobear_GCTF_path = '/path/to/tomobear/output/X_GCTFCtfphaseflipCTFCorrection_1';
particle_tbl_filepath = '/path/to/particles.tbl';

% path to folder with defocus.txt files per each tomogram
defocus_files_dirpath = '/path/to/defocus_files/folder';

%% Load input files
tbl = dread(particle_tbl_filepath);

%% Generate defocus.txt file in IMOD style
tid = unique(tbl(:,20));
for i=1:length(tid)
	[~,lst] = system(['ls ' sprintf('%s/tomogram_%03d/slices/*_gctf.log', tomobear_GCTF_path, tid(i))]);
	lst = split(strtrim(lst));
	def = fopen(sprintf('%s/tomogram_%03d_defocus.txt', defocus_files_dirpath, tid(i)), 'w');
	for j=1:length(lst)
    	fid = fopen(cell2mat(lst(j)), 'r');
    	line_divided_text = textscan(fid, "%s", "delimiter", "\n");
    	final_values = line_divided_text{1}{contains(line_divided_text{1}, "Final Values")};
    	res_lim = line_divided_text{1}{contains(line_divided_text{1}, "Resolution limit")};
    	fin_val_splt = strsplit(final_values);res_lim_splt = strsplit(res_lim);
    	fclose(fid);
    	fprintf(def, '%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n',j,str2double(cell2mat(fin_val_splt(1))),str2double(cell2mat(fin_val_splt(2))), str2double(cell2mat(fin_val_splt(3))), 0, str2double(cell2mat(fin_val_splt(4))),str2double(cell2mat(res_lim_splt(7))));
	end
	fclose(def);
end