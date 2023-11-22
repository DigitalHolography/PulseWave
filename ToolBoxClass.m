classdef ToolBoxClass
     properties
        %Path of the PW dir and the output dir inside
        PW_path_main char 
        PW_path_dir char
        PW_path_png char 
        PW_path_eps char 
        PW_path_txt char 
        PW_path_avi char 
        PW_path_mp4 char
        PW_path_pulswave char 
        main_foldername char
        PW_folder_name char
        stride double
        fs double
        type char
        f1 double
        f2 double
        minPCA double
        maxPCA double
        x_barycentre double
        y_barycentre double
        ScalingFactorVelocityInPlane double
        ScalingFactorVelocityCRA_AVG double
        ScalingFactorVelocityCRA_RMS double
        ARI_hue_max double
        ARI_hue_min double
        ARI_inflexion_point_hue double
        ARI_slope_hue double
        ARI_val_max double
        ARI_val_min double
        ARI_inflexion_point_val double
        ARI_slope_val double

     end


     methods
         function obj = ToolBoxClass(path)

             PW_params = Parameters(path);

             %% Creating paths
            idx = 0 ;
            split_path = strsplit(path, '\');
            obj.main_foldername = split_path{end-1};
            obj.PW_path_main = fullfile(path,'pulsewave');

            if ~exist(obj.PW_path_main,'dir')
                mkdir(obj.PW_path_main);
            end

            PW_folder_name = strcat( obj.main_foldername, '_PW');

            while (exist(fullfile(obj.PW_path_main, sprintf('%s_%d', PW_folder_name,idx)), 'dir'))
                idx = idx + 1 ;
            end
            obj.PW_folder_name = sprintf('%s_%d', PW_folder_name,idx);
            obj.PW_path_dir = fullfile(obj.PW_path_main, obj.PW_folder_name) ;
            obj.PW_path_png = fullfile(obj.PW_path_dir, 'png');
            obj.PW_path_eps = fullfile(obj.PW_path_dir, 'eps');
            obj.PW_path_txt = fullfile(obj.PW_path_dir, 'txt');
            obj.PW_path_avi = fullfile(obj.PW_path_dir, 'avi');
            obj.PW_path_mp4 = fullfile(obj.PW_path_dir, 'mp4');

            %% Reading Cache Parameters from .mat
            dir_path_mat = fullfile(path,'mat');
            file_path_mat = fullfile(dir_path_mat,[obj.main_foldername,'.mat']);
            %[~,filename_mat,~] = fileparts(file_path_mat);

            if exist(file_path_mat) % .mat with cache from holowaves is present, timeline can be computed
                disp('reading cache parameters');
                load(file_path_mat,'cache') ;
                obj.stride = cache.batch_stride ;
                obj.fs = (cache.Fs)/1000;
                obj.type = cache.time_transform.type;
                obj.f1 = cache.time_transform.f1;
                obj.f2 = cache.time_transform.f2;
                obj.minPCA = cache.time_transform.min_PCA;
                obj.maxPCA = cache.time_transform.max_PCA;
                disp('done.')

            elseif exist(fullfile(path,[obj.main_foldername,'.mat']))
                disp('reading cache parameters');
                load(fullfile(path,[obj.main_foldername,'.mat']),'cache') ;
                obj.stride = cache.batch_stride ;
                obj.fs = (cache.Fs)/1000; %conversion in kHz
                obj.type = cache.time_transform.type;
                obj.f1 = cache.time_transform.f1;
                obj.f2 = cache.time_transform.f2;
                obj.minPCA = cache.time_transform.min_PCA;
                obj.maxPCA = cache.time_transform.max_PCA;
                disp('done.')
            else
                disp('no mat file found');
                obj.stride = 0 ;
                obj.fs = 0;
                obj.type = 'None';
                obj.f1 = 0;
                obj.f2 = 0;
                obj.minPCA = 0;
                obj.maxPCA = 0;
            end

            %% Calculation of the Sacling Factors

%             obj.ScalingFactorVelocityInPlane = 1000 * 1000 * PW_params.lambda / PW_params.opticalIndex * (3/PW_params.theta)^(1/2); % 1000 for kHz -> Hz and 1000 for m -> mm
            obj.ScalingFactorVelocityInPlane = 30;
            obj.ScalingFactorVelocityCRA_AVG  = 1000 * 1000 * PW_params.lambda / PW_params.opticalIndex * (PW_params.theta/2); % 1000 for kHz -> Hz and 1000 for m -> mm
            obj.ScalingFactorVelocityCRA_RMS  = 1000 * 1000 * PW_params.lambda / PW_params.opticalIndex * (1/(2+2*(PW_params.theta^3)/3))^(1/2); % 1000 for kHz -> Hz and 1000 for m -> mm
            %% Parameters the color maps 

            obj.ARI_hue_max = 0;
            obj.ARI_hue_min = 0.15;
            obj.ARI_inflexion_point_hue = 0.7;
            obj.ARI_slope_hue = 10;
            obj.ARI_val_max = 0.4;
            obj.ARI_val_min = 1;
            obj.ARI_inflexion_point_val = 1;
            obj.ARI_slope_val = 10;

         end
        



     end
end