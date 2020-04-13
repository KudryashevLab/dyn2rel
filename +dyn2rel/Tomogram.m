classdef Tomogram < handle
    % Holds a Tomogram's info: ST, XF, ALI, DEFOCUS files, and its size.
    %
    % Tomogram Properties:
    %    index      - Tomogram index.
    %    file_stack - .st, .ali or .mrc unbinned stack.
    %    tomo_size  - Size of the unbinned reconstructed tomogram.
    %    HD_ix      - Index of the High-Dose projection.
    %    pix_siz    - Pixel size in angstroms.
    %    Txf        - Affine matrix, read from .xf file.
    %    Ttlt       - Affine matrix, read from .tlt file.
    %    defocus    - Defocus information.
    %
    % Tomogram Methods:
    %    info              - Shows the info of the tomogram.
    %    set_stack         - Sets the stack used to generate the tomogram.
    %    read_ali_files    - Reads the alignment information of the stack.
    %    set_tomogram_size - Set the 3D size of the reconstructed tomogram.
    %    set_defocus       - Set the defocus values for the tomogram.
    
    properties
        index      % Tomogram index.
    end
    
    properties(SetAccess=private)
        file_stack % .st, .ali or .mrc unbinned stack.
        tomo_size  % Size of the unbinned reconstructed tomogram.
        HD_ix      % Index of the High-Dose projection.
        pix_siz    % Pixel size in angstroms.
        Txf        % Affine matrix, read from .xf file.
        Ttlt       % Affine matrix, read from .tlt file.
        defocus    % Defocus information.
    end
    
    methods
        function obj = Tomogram
            % Creates an empty Tomogram class.
            obj.index = 0;
            obj.tomo_size = [];
            obj.defocus = [];
            obj.file_stack = [];
            obj.Txf = eye(4,4);
            obj.Ttlt = eye(4,4);
            obj.HD_ix = 1;
            obj.pix_siz = 1;
        end
        
        function info(obj)
            % Shows the info of the tomogram.
            fprintf(['Tomogram Index: %3d\n'],obj.index);
            if( isempty(obj.file_stack) )
                fprintf('No stack file set\n');
            else
                fprintf('Stack file: %s\n',obj.file_stack);
            end
            if( isempty(obj.tomo_size) )
                fprintf('Unknown tomogram size\n');
            else
                fprintf('Tomogram size: %d %d %d\n',obj.tomo_size(1),obj.tomo_size(2),obj.tomo_size(3));
            end
            fprintf(['Pixel size: %.3f' 197 '\n'],obj.pix_siz);
            fprintf('High-Dose projection: %d\n',obj.HD_ix);
            fprintf('Ttlt: %6.3f %6.3f %6.3f | %5.1f\n',obj.Ttlt(1,1),obj.Ttlt(1,2),obj.Ttlt(1,3),obj.Ttlt(1,4));
            fprintf('      %6.3f %6.3f %6.3f | %5.1f\n',obj.Ttlt(2,1),obj.Ttlt(2,2),obj.Ttlt(2,3),obj.Ttlt(2,4));
            fprintf('      %6.3f %6.3f %6.3f | %5.1f\n',obj.Ttlt(3,1),obj.Ttlt(3,2),obj.Ttlt(3,3),obj.Ttlt(3,4));
            fprintf('Txf:  %6.3f %6.3f %6.3f | %5.1f\n',obj.Txf(1,1),obj.Txf(1,2),obj.Txf(1,3),obj.Txf(1,4));
            fprintf('      %6.3f %6.3f %6.3f | %5.1f\n',obj.Txf(2,1),obj.Txf(2,2),obj.Txf(2,3),obj.Txf(2,4));
            fprintf('      %6.3f %6.3f %6.3f | %5.1f\n',obj.Txf(3,1),obj.Txf(3,2),obj.Txf(3,3),obj.Txf(3,4));
            if( isempty(obj.defocus) )
                fprintf('No defocus information\n');
            else
                obj.defocus.info;
            end
        end
        
        function obj = set_stack(obj,filename,pix_siz)
            % Sets the stack used to generate the tomogram.
            %   .set_stack(FILENAME) FILENAME is a MRC file with the
            %   following extensions: .st, .ali, .mrc, .mrcs
            %   .set_stack(FILENAME,PIX_SIZ) uses PIX_SIZ, in angstroms,
            %   instead of the one read from the MRC header.
            %
            if(nargin > 1)
                if( exist(filename,'file') )
                    obj.file_stack = filename;
                else
                    error('%s does not exists.',filename);
                end
            end
            
            if(nargin > 2)
                obj.pix_siz = pix_siz;
            else
                obj.pix_siz = obj.read_apix_mrc(obj.file_stack);
            end
            
        end
        
        function obj = read_ali_files(obj,tlt_file,xf_file,hd_proj)
            % Reads the alignment information of the stack.
            %   .read_ali_files(TLT_FILE) Reads a .tlt file Assumes the 
            %   High-dose projection is at 0 degrees.
            %   .read_ali_files(TLT_FILE,XF_FILE) Reads a .tlt and .xf 
            %   files. Assumes the High-dose projection is at 0 degrees.
            %   .read_ali_files(TLT_FILE,XF_FILE,HD_PROJ) Reads a .tlt and
            %   .xf files and uses the tilt angle at the index HD_PROJ, 
            %   which corresponds to the High-dose projection.
            %   .read_ali_files(TLT_FILE,[].HD_PROJ) The same as the
            %   previous one, but it ignores the .xf file (the stack is a
            %   .ali file).
            
            if( nargin < 4 )
                hd_proj = 0;
            end
            
            if( nargin < 3 )
                xf_file= [];
            end
            
            % Read TLT:
            tlt_angles = obj.read_tlt(tlt_file);
            if( hd_proj < 1 )
                obj.HD_ix = round(length(tlt_angles)/2);
            else
                obj.HD_ix = hd_proj;
            end
            ct = cosd(tlt_angles(obj.HD_ix));
            st = sind(tlt_angles(obj.HD_ix));
            obj.Ttlt = [ ct 0 -st 0; 0 1 0 0; st 0 ct 0; 0 0 0 1];
            
            % Read XF:
            if( isempty(xf_file) )
                obj.Txf = eye(4,4);
            else
                xf_raw = obj.read_xf(xf_file);
                obj.Txf = eye(4,4);
                %[xf(1) xf(2) 0 xf(5); xf(3) xf(4) 0 xf(6); 0 0 1 0; 0 0 0 1];
                obj.Txf(1,:) = [xf_raw(1,obj.HD_ix) xf_raw(2,obj.HD_ix) 0 xf_raw(5,obj.HD_ix)];
                obj.Txf(2,:) = [xf_raw(3,obj.HD_ix) xf_raw(4,obj.HD_ix) 0 xf_raw(6,obj.HD_ix)];
            end
        end
    
        function obj = set_tomogram_size(obj,tomo_size)
            % Set the 3D size of the reconstructed tomogram.
            %   .set_tomogram_size(TOMO_SIZE) TOMO_SIZE must be a 3-element
            %   vector.
            
            if( numel(tomo_size) ~= 3 )
                error('The tomogram size must have three values');
            else
                obj.tomo_size = reshape(tomo_size,[1 3]);
            end
        end
        
        function obj = set_defocus(obj,def_U,def_V,def_A)
            % Set the defocus values for the tomogram.
            %   .set_defocus(DEF_VAL) set DEF_VAL as the defocus without
            %   astigmatism. In angstroms.
            %   .set_defocus(DEFOCUS_FILE) read the defocus value from the
            %   .defocus file DEFOCUS_FILE.
            %   .set_defocus(DEF_U,DEF_V,DEF_A) sets the defocus with
            %   astigmatism.
            if( nargin == 2 )
                if( isnumeric( def_U ) )
                    obj.defocus = dyn2rel.CTF;
                    obj.defocus.U = def_U;
                    obj.defocus.V = def_U;
                    obj.defocus.A = 0;
                elseif( ischar( def_U ) )
                    data = obj.read_defocus(def_U);
                    obj.defocus = dyn2rel.CTF;
                    obj.defocus.U = data(obj.HD_ix);
                    obj.defocus.V = data(obj.HD_ix);
                    obj.defocus.A = 0;
                end
            end
            
            if( nargin == 4 )
                obj.defocus = dyn2rel.CTF;
                obj.defocus.U = def_U;
                obj.defocus.V = def_V;
                obj.defocus.A = def_A;
            end
            
        end
        
        function obj = set_bfactor(obj,bfactor)
            % Set the bfactor values for the tomogram.
            %   .set_bfactor(BFACTOR) set BFACTOR as the BFactor.
            if( ~isempty( obj.defocus ) )
                obj.defocus.bfactor = bfactor;
            end
            
        end
    end
    
    
    methods(Static)
        function tomos_list = create_tomos_list(N_TOMOS)
            % Creates an array with Tomograms.
            %   TOMOS_LIST = .CREATE_TOMOS_LIST(N_TOMOS) creates an array with N_TOMOS
            %   empty tomograms (dyn2rel.Tomogram class).
            % See also dyn2rel.Tomogram

            tomos_list(1:N_TOMOS) = dyn2rel.Tomogram;

            for i = 1:N_TOMOS
                tomos_list(i) = dyn2rel.Tomogram;
                tomos_list(i).index = i;
            end

        end
    end

    
    methods(Access=private,Static)
        function apix = read_apix_mrc(filename)
            fp = fopen(filename,'r');
            fseek(fp, 28, 'bof');
            mx   = fread(fp,1,'uint32');
            fseek(fp, 40, 'bof');
            xlen = fread(fp,1,'float32');
            fclose(fp);
            if( xlen == 0 ) 
                apix = 1;
            else
                apix = xlen/mx;
            end
            
        end
        
        function data = read_tlt(filename)
            fp = fopen(filename,'r');
            data = fscanf(fp,'%f');
            fclose(fp);
            data = data(:);
        end
        
        function data = read_xf(filename)
            fp = fopen(filename,'r');
            data = fscanf(fp,'%f',[6 inf]);
            fclose(fp);
        end
        
        function data = read_defocus(filename)
            fp = fopen(filename,'r');
            data = [];
            tline = fgetl(fp);
            while( length(tline) > 1 )
                splitted = sscanf(tline,'%f',inf);
                data(end+1) = 10*splitted(5);
                tline = fgetl(fp);
            end
            fclose(fp);
        end
    end
end 
