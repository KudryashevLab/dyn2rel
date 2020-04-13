classdef CTF < handle
    % Holds the CTF's information.
    
    properties
        U
        V
        A
        volt
        sph_ab
        ac
        bfactor
    end
    
    methods
        function obj = CTF
            % Creates CTF class with the default values.
            obj.U       = 0;
            obj.V       = 0;
            obj.A       = 0;
            obj.volt    = 300;
            obj.sph_ab  = 2.7;
            obj.ac      = 0.07;
            obj.bfactor = 10.97;
        end
        
        function info(obj)
            % Shows the info of the CTF.
            fprintf(['Defocus U: %.0f' 197 '\n'],obj.U);
            fprintf(['Defocus V: %.0f' 197 '\n'],obj.V);
            fprintf('Defocus angle: %.2f\n',obj.A);
            fprintf('Voltage: %d\n',obj.volt);
            fprintf('Spherical aberration: %.2f\n',obj.sph_ab);
            fprintf('Amplitud Contrast: %.2f\n',obj.ac);
            fprintf('B-Factor: %.3f\n',obj.bfactor);
        end
    end
    
    methods(Access=private,Static)
        
    end
end 
