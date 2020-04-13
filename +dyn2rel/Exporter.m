classdef Exporter < handle
    %   Creates a RELION stack from a DYNAMO tbl and an Array of Tomograms.
    %   For each entry in the table, it calculates its location in the
    % projection according to:
    %     P_position = tbl_mul*(table(i,[24 25 26])+table(i,[4 5 6]));
    % after solving rotations and orientations, it crops a patch of size
    % out_siz and it bins it according to out_bin before saving it.
    %
    % Exporter Properties:
    %   xcor_sel - fraction of particles to extract, according to CC (0 to 1).
    %   invert   - (bool) invert contrast. Default: true.
    %   out_siz  - size of the unbinned .mrcs file.
    %   out_bin  - binning of the .mrcs.
    %   tbl_mul  - scale of the table.
    %   split_h  - split particle in random halves.
    %
    % Exporter Methods:
    %    exec              - exports the dynamo project as a relion project
    
    properties
        xcor_sel % fraction of particles to extract, according to CC (0 to 1).
        invert   % (bool) invert contrast. Default: true.
        out_siz  % size of the unbinned .mrcs file.
        out_bin  % binning of the .mrcs.
        tbl_mul  % binning of the table.
        split_h  % split halves. < 1: no random subset, 1: even-odd, > 1: table column.
        skip_gb  % skip goldbeads factor.
    end
    
    properties(Access=private)
        use_ctf
        detector_size
        magnification
    end
    
    methods
        function obj = Exporter
            % Creates an empty Tomogram class.
            obj.xcor_sel = 0.25;
            obj.invert   = true;
            obj.out_siz  = 200;
            obj.out_bin  = 0;
            obj.tbl_mul  = 1;
            obj.split_h  = 0;
            obj.skip_gb  = 0;
            
            obj.use_ctf       = false;
            obj.detector_size = 14;
            obj.magnification = 1;
        end
        
        function out_tbl = exec(obj,out_pfx,tomo_list,table)
            % Exports the dynamo project as a relion project.
            %   .exec(OUT_PFX,TOMO_LIST.TABLE) Creates a OUT_PFX.star and a
            %   OUT_PFX.mrcs file from TOMO_LIST, an array of class
            %   Tomogram, according to the dynamo table TABLE. TABLE can be
            %   a variable or a filename.
            %   
            % See also dyn2rel.Tomogram
            
            if( nargout > 0 )
                out_tbl = [];
            end
            
            fprintf('Exporting DYNAMO project to RELION\n');
            
            %%%%%%% Get best particles.
            if( ischar(table) )
                table_w = obj.filter_table(dread(table),obj.xcor_sel);
            else
                table_w = obj.filter_table(table,obj.xcor_sel);
            end
            
            %%%%%%% Inform about the tomograms on the list.
            tomo_idxs = unique(table_w(:,20));
            fprintf('%d tomograms available (from 1 to %d)\n',length(tomo_list),length(tomo_list));
            fprintf('%d tomograms used in table:',length(tomo_idxs));
            fprintf(' %d',tomo_idxs);
            fprintf('\n');
            
            %%%%%%% Check all tomograms have CTF information:
            obj.use_ctf = true;
            for i = 1:length(tomo_list)
                if( isempty(tomo_list(i).defocus) )
                    obj.use_ctf = false;
                end
            end
            if( ~obj.use_ctf )
                fprintf('At least one tomograms does not have CTF info, ignoring all CTF info.\n');
            end
            
            %%%%%%% Inform random halves policy:
            if( obj.split_h < 1 )
                fprintf('No random subsets.\n');
            elseif( obj.split_h == 1 )
                fprintf('Using EVEN-ODD random subsets.\n');
            else
                fprintf('Random subset according to column %d of the table.\n',obj.split_h);
            end
            
            %%%%%%% Init .STAR file:
            fp_star = fopen([out_pfx '.star'],'w');
            obj.write_star_header(fp_star);
            
            %%%%%%% Init .MRCS file:
            fp_mrcs = fopen([out_pfx '.mrcs'],'wb');
            obj.init_mrcs(fp_mrcs)
            
            %%%%%%% Iterate over each tomogram:
            fprintf('Extracting %d particles:\n',size(table_w,1));
            particle_counter = 1;
            skipped_particles = 0;
            skipped_goldbeads = 0;
            mrcs_min = +inf;
            mrcs_max = -inf;
            mrcs_avg = 0;
            for tomo_ix = 1:length(tomo_idxs)
                tomo_id = tomo_idxs(tomo_ix);
                tomo_ix_list = obj.get_tomo_index(tomo_id,tomo_list);
                if( tomo_ix_list>0 )
                    table_w_t = table_w( table_w(:,20)==tomo_id, :);
                    fprintf('    Tomogram %3d (%6d particles): %5.1f%%',tomo_id,size(table_w_t,1),0);
                    [mrc, mrc_size] = obj.read_mrc(tomo_list(tomo_ix_list));
                    if( obj.skip_gb )
                        std_mrc = std(mrc(:));
                    end
                    pos = obj.project_positions(tomo_list(tomo_ix_list),table_w_t,mrc_size);
                    obj.magnification = obj.detector_size*10000/( tomo_list(tomo_ix_list).pix_siz *(2^(obj.out_bin)));
                    
                    %%%%%%% Iterate over each particle:
                    for part_ix = 1:size(table_w_t,1)
                        fprintf('\b\b\b\b\b\b%5.1f%%',100*part_ix/size(table_w_t,1));

                        %%%%%%% Crop projections:
                        proj = dynamo_crop(mrc,obj.out_siz,round(pos(part_ix,[1 2]))-1,1);
                        if( ~isempty(proj) )
                            
                            %%%%%%% Check for goldbeads:
                            process_particle = true;                            
                            if( obj.skip_gb > 0 )
                                tmp_proj = medfilt2(proj,[5 5],'symmetric');
                                if( (min(tmp_proj(:))-mean(tmp_proj(:))) < -obj.skip_gb*std_mrc )
                                    skipped_goldbeads = skipped_goldbeads + 1;
                                    process_particle = false;
                                end
                            end
                            
                            %%%%%%% Skip if it is a golbead:
                            if( process_particle )
                                
                                if( nargout > 0 )
                                    out_tbl = [out_tbl; table_w_t(part_ix,:)];
                                end

                                %%%%%%% Bin, Invert and Normalize:
                                if( obj.out_bin > 0 )
                                    proj = dbin( proj, obj.out_bin );
                                end
                                proj = (proj-mean(proj(:)))/std(proj(:));
                                if( obj.invert > 0 )
                                    proj = -proj;
                                end

                                %%%%%%% Gather statistics:
                                mrcs_min = min(mrcs_min,min(proj(:)));
                                mrcs_max = max(mrcs_max,max(proj(:)));
                                mrcs_avg = mrcs_avg + mean(proj(:));

                                %%%%%%% Save proj in stack:
                                fwrite(fp_mrcs,proj,'single');

                                %%%%%%% Save entry in STAR file
                                obj.write_star_entry(fp_star,tomo_list(tomo_ix_list),table_w_t(part_ix,:),pos(part_ix,:),particle_counter,[out_pfx '.mrcs']);
                                particle_counter = particle_counter + 1;
                            end
                        else
                            skipped_particles = skipped_particles + 1;
                        end
                    end

                    fprintf('\n');
                else
                    fprintf('    No tomogram with ID %d found on tomogram list. Skipping.\n',tomo_id);
                end                
                
                
            end
            
            %%%%%%% Final report:
            fprintf('Finished. Files %s.star and %s.mrcs were created.\n',out_pfx,out_pfx);
            fprintf('    Particles cropped: %d.\n',particle_counter-1);
            if( skipped_particles > 0 )
                fprintf('    Particles skipped: %d.\n',skipped_particles);
            end
            if( skipped_goldbeads > 0 )
                fprintf('    Goldbeads skipped: %d.\n',skipped_goldbeads);
            end
            
            %%%%%%% Close .STAR file:
            fclose(fp_star);
            
            %%%%%%% Close .MRCS file:
            out_size = obj.out_siz / (2^obj.out_bin);
            obj.finish_mrcs(fp_mrcs,[out_size out_size particle_counter-1], mrcs_min, mrcs_max, mrcs_avg, tomo_list(end).pix_siz);
            fclose(fp_mrcs);
        end
    end
    
    methods(Access=private)
        function write_star_header(obj,fp_star)
            i = 1;
            fprintf(fp_star, '\ndata_\n\nloop_\n');
            fprintf(fp_star, '_rlnImageName           #%d\n',i); i = i+1;
            fprintf(fp_star, '_rlnOriginX             #%d\n',i); i = i+1;
            fprintf(fp_star, '_rlnOriginY             #%d\n',i); i = i+1;
            fprintf(fp_star, '_rlnAngleRot            #%d\n',i); i = i+1; % First Euler Angle
            fprintf(fp_star, '_rlnAngleTilt           #%d\n',i); i = i+1; % Second Euler Angle
            fprintf(fp_star, '_rlnAnglePsi            #%d\n',i); i = i+1; % Third Euler Angle
            fprintf(fp_star, '_rlnCoordinateX         #%d\n',i); i = i+1; % Position of particle in micrograph
            fprintf(fp_star, '_rlnCoordinateY         #%d\n',i); i = i+1; % Position of particle in micrograph
            fprintf(fp_star, '_rlnCoordinateZ         #%d\n',i); i = i+1; % Position (delta) of particle in micrograph
            
            if( obj.use_ctf )
                fprintf(fp_star, '_rlnDefocusU            #%d\n',i); i = i+1;
                fprintf(fp_star, '_rlnDefocusV            #%d\n',i); i = i+1;
                fprintf(fp_star, '_rlnDefocusAngle        #%d\n',i); i = i+1;
                fprintf(fp_star, '_rlnVoltage             #%d\n',i); i = i+1;
                fprintf(fp_star, '_rlnSphericalAberration #%d\n',i); i = i+1;
                fprintf(fp_star, '_rlnAmplitudeContrast   #%d\n',i); i = i+1;
                fprintf(fp_star, '_rlnMagnification       #%d\n',i); i = i+1;
                fprintf(fp_star, '_rlnDetectorPixelSize   #%d\n',i); i = i+1;
                fprintf(fp_star, '_rlnCtfBfactor          #%d\n',i); i = i+1;
            end
            
            if( obj.split_h > 0 )
                fprintf(fp_star, '_rlnRandomSubset        #%d\n',i);
            end
        end
        
        function write_star_entry(obj,fp_star,tomo_class,table_entry,pos,particle_ix,mrcs_name)
            
            off = [0;0];
            
            Ttot = tomo_class.Ttlt*tomo_class.Txf;
            R = dynamo_euler2matrix(table_entry(7:9));
            euler = obj.rot_M2eZYZ(R*Ttot(1:3,1:3));
            
            fprintf(fp_star, '%06d@%s',particle_ix,mrcs_name);
            fprintf(fp_star, '\t%f\t%f',off(1),off(2));
            fprintf(fp_star, '\t%12.6f\t%12.6f\t%12.6f',euler(1),euler(2),euler(3));
            fprintf(fp_star, '\t%.1f\t%.1f\t%.1f',round(pos(1)),round(pos(2)),round(pos(3)));
            
            if( obj.use_ctf )
                dz = tomo_class.pix_siz * pos(3);
                fprintf(fp_star, '\t%12.6f',tomo_class.defocus.U+dz);
                fprintf(fp_star, '\t%12.6f',tomo_class.defocus.V+dz);
                fprintf(fp_star, '\t%f',tomo_class.defocus.A);
                fprintf(fp_star, '\t%f',tomo_class.defocus.volt);
                fprintf(fp_star, '\t%f',tomo_class.defocus.sph_ab);
                fprintf(fp_star, '\t%f',tomo_class.defocus.ac);
                fprintf(fp_star, '\t%f',obj.magnification);
                fprintf(fp_star, '\t%f',obj.detector_size);
                fprintf(fp_star, '\t%f',tomo_class.defocus.bfactor);
            end
            
            if( obj.split_h == 1 )
                fprintf(fp_star, '\t%12d',mod(particle_ix,2)+1);
            elseif( obj.split_h > 1 )
                fprintf(fp_star, '\t%12d',table_entry(obj.split_h));
            end
            
            fprintf(fp_star, '\n');
        end
        
        function pos = project_positions(obj,tomo_class,table,mrc_size)
            
            tomo_center = tomo_class.tomo_size/2 + 1;
            stack_center = mrc_size/2 + 1;
            
            if( (stack_center(1)>stack_center(2)) && (tomo_center(1)<tomo_center(2)) )
                tomo_center([1 2 3]) = tomo_center([2 1 3]);
            elseif( (stack_center(1)<stack_center(2)) && (tomo_center(1)>tomo_center(2)) )
                tomo_center([1 2 3]) = tomo_center([2 1 3]);
            end
            
            pos = obj.tbl_mul*( table(:,[24 25 26]) + table(:,[4 5 6]) );
            
            T = inv(tomo_class.Ttlt*tomo_class.Txf);

            pos(:,1) = pos(:,1) - tomo_center(1);
            pos(:,2) = pos(:,2) - tomo_center(2);
            pos(:,3) = pos(:,3) - tomo_center(3);
            
            pos = [pos ones(size(pos(:,1)))]*(T');
            
            pos(:,1) = pos(:,1) + stack_center(1);
            pos(:,2) = pos(:,2) + stack_center(2);
            pos = pos(:,1:3);
        end
    end
    
    methods(Access=private,Static)
        function table_w = filter_table(table,th)
            if( th > 0 && th < 1 )
                table_w = table(table(:,10)>quantile(table(:,10),1-th),:);
            else
                table_w = table;
            end
        end
        
        function [mrc, mrc_size] = read_mrc(tomo_class)
            fp = fopen(tomo_class.file_stack,'r');
            fseek(fp, 0, 'bof');
            stack_size = fread(fp,3,'uint32');
            stack_mode = fread(fp,1,'uint32');
            
            % Check mode
            if( stack_mode ~= 2 )
                error('Stack mode is %d, the mode supported is 2 (32-bit float)\n',stack_mode);
            end
            
            % Check size
            if( (stack_size(3) > 1) && (stack_size(3) < tomo_class.HD_ix ) )
                fclose(fp);
                error('Stack has %d frames, requested frame number %d\n',stack_size(3), tomo_class.HD_ix);
            end
            
            % Skip to frame to be read
            if( stack_size(3) == 1 )
                fseek(fp, 1024, 'bof');
            else
                fseek(fp, 1024 + 4*(tomo_class.HD_ix-1)*stack_size(1)*stack_size(2), 'bof');
            end
            
            % Read frame:
            mrc = fread(fp,[stack_size(1) stack_size(2)],'single');
            mrc_size = size(mrc);
        end
        
        function data = read_tlt(filename)
            fp = fopen(filename,'r+');
            data = fscanf(fp,'%f');
            fclose(fp);
            data = data(:);
        end
        
        function data = read_xf(filename)
            fp = fopen(filename,'r+');
            data = fscanf(fp,'%f',[6 inf]);
            fclose(fp);
        end
        
        function data = read_defocus(filename)
            fp = fopen(filename,'r+');
            data = [];
            tline = fgetl(fp);
            while( length(tline) > 1 )
                splitted = sscanf(tline,'%f',inf);
                data(end+1) = 10*splitted(5);
                tline = fgetl(fp);
            end
            fclose(fp);
        end
    
        function init_mrcs(fp_mrcs)
            fseek(fp_mrcs, 0, 'bof');
            buffer = uint8(zeros(1024,1));
            fwrite(fp_mrcs,buffer,'uint8');
        end
        
        function finish_mrcs(fp_mrcs, st_size, mrcs_min, mrcs_max, mrcs_avg, pix_size)
            fseek(fp_mrcs, 0, 'bof');
            first_values = uint32([st_size 2]);
            map_values   = uint32([1 2 3]);
            stat_values  = single([mrcs_min,mrcs_max,mrcs_avg]);
            xyz_len      = single(st_size);
            xyz_len      = pix_size*xyz_len;

            fseek(fp_mrcs, 0, 'bof');
            fwrite(fp_mrcs,first_values,'uint32');

            fseek(fp_mrcs, 28, 'bof');
            fwrite(fp_mrcs,first_values,'uint32');
            
            fseek(fp_mrcs, 40, 'bof');
            fwrite(fp_mrcs,xyz_len,'float32');
            
            fseek(fp_mrcs, 64, 'bof');
            fwrite(fp_mrcs,map_values,'uint32');

            fseek(fp_mrcs, 76, 'bof');
            fwrite(fp_mrcs,stat_values,'single');
        end
        
        function eu = rot_M2eZYZ(R)
            
            eu = zeros(1,3);
            tol = 5e-5;

            if( R(3,3) < (1-tol))

                if( R(3,3) > (tol-1) )

                    % GENERAL CASE
                    eu(2) = acos ( R(3,3) )*180/pi;
                    eu(1) = atan2( R(2,3), R(1,3) )*180/pi;
                    eu(3) = atan2( R(3,2),-R(3,1) )*180/pi;

                else

                    % r22 <= -1
                    eu(2) = 180;
                    eu(1) = -atan2( R(2,1), R(2,2) )*180/pi;
                    eu(3) = 0;

                end
            else

                % r22 <= -1
                eu(2) = 0;
                eu(1) = atan2( R(2,1), R(2,2) )*180/pi;
                eu(3) = 0;

            end

        end
        
        function ix = get_tomo_index(tomo_id,tomo_list)
            ix = 0;
            for i = 1:length(tomo_list)
                if( tomo_id == tomo_list(i).index )
                    ix = i;
                end
            end
        end
    end
end 
