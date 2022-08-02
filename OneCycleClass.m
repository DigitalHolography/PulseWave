classdef OneCycleClass
    properties
        data (1,:) cell
        directory char
        filenames (1,:) cell
        nbFiles {mustBeNumeric , mustBePositive}  
    end

    methods
        function obj = OneCycleClass(files,path,refvideosize)
            arguments
                files
                path
                refvideosize = []
            end
            obj.directory = path;
            obj.filenames = files;
            obj.nbFiles = length(files) ; 
            obj.data = cell(1,obj.nbFiles) ;
            for ii = 1 : obj.nbFiles
                currentFilePath = fullfile(obj.directory,obj.filenames{ii});
                [filepath,name,ext] = fileparts(currentFilePath);
                if (ext == '.avi')
                    disp(['reading : ',filepath,name,ext]);
                    V = VideoReader(fullfile(path,files{ii}));
                    video = zeros(V.Height, V.Width, V.NumFrames);
                    for n = 1 : V.NumFrames
                        video(:,:,n) = rgb2gray(read(V, n));
                    end
                    obj.data{ii} = video;
                elseif (ext == '.raw')
                    disp(['reading : ',filepath,name,ext]);
                    fileID = fopen(currentFilePath);
                    video = fread(fileID,'float32');
                    fclose(fileID)
                    obj.data{ii} = reshape(video,refvideosize);
                else
                    disp([filepath,name,ext,' : non recognized video format']);
                end
            end
        end

        function [sys_index_list_cell, mask] = getSystole(obj)
            sys_index_list_cell = cell(obj.nbFiles) ;
            for i = 1:obj.nbFiles
                [sys_index_list_cell{i}, mask] = find_systole_index(obj.data{i});
            end
        end

        function writeCut(obj, mask, sys_index_list_cell, Ninterp, add_infos)
            idx = 0 ;
            while (exist(fullfile(obj.directory, sprintf("one_cycle_%d",idx)), 'dir'))
                idx = idx + 1 ;
            end
            one_cycle_dir = fullfile(obj.directory, sprintf("one_cycle_%d",idx)) ; 
            mkdir(one_cycle_dir);
            for n = 1 : obj.nbFiles
                one_cycle_video = create_one_cycle(obj.data{n}, mask, sys_index_list_cell{n}, Ninterp) ;
                [~,name,ext]=fileparts(obj.filenames{n}) ;
                w = VideoWriter(fullfile(one_cycle_dir,strcat(name,'_one_cycle',ext))) ;
                open(w)
                for j = 1:size(one_cycle_video,3)
                    writeVideo(w,one_cycle_video(:,:,j)) ;
                end
                close(w);
                disp(fullfile('./',strcat(name,'.mat'))) ; 
%                 cache_exists = exist(fullfile('./',strcat(name,'.mat'))) ; 
%                 disp(cache_exists)
                if add_infos
                    pulseAnalysis(one_cycle_video,one_cycle_dir, name,sys_index_list_cell{n});
                end 
                
            end



        end
    end
end