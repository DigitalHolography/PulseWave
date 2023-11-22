classdef OneCycleClass
    properties
        reference (1,:) cell % M0 AVI
        reference_norm (1,:) cell % M0 norm AVI
        reference_interp (1,:) cell % M0 AVI
        reference_norm_interp (1,:) cell % M0 norm AVI
        dataM2M0 (1,:) cell % RMS M2/M0
        dataM1M0 (1,:) cell % AVG M1/M0
        dataM2M0_interp (1,:) cell % RMS M2/M0
        dataM1M0_interp (1,:) cell % AVG M1/M0
        dataM0 (1,:) cell % M0 raw
        dataM0_interp (1,:) cell % M0 raw
        dataM1 (1,:) cell % M1 raw
        dataM2 (1,:) cell % M2 raw
        dataSH (1,:) cell % SH raw
        dataSH_interp (1,:) cell % SH raw
        directory char
        nbFiles {mustBeNumeric , mustBePositive}
        filenames (1,:) cell
        k double %interpolaton parameter
        % For sectioning
        video_loaded
        pictureSection
        videoSection
        flag_FlatField
        flag_SH_analysis
        flag_PulseWave_analysis
        isdone_flatfield 
        ToolBoxmaster ToolBoxClass
        




    end

    methods
        function obj = OneCycleClass(path)

            arguments
                path
            end
            
            obj.directory = path;

            tmp_idx = regexp(path, '\');
            obj.filenames = {path(tmp_idx(end-1)+1:end-1)};

            obj.nbFiles = 1 ;
            obj.reference = cell(1,obj.nbFiles) ;
            obj.reference_norm = cell(1,obj.nbFiles);
            obj.dataM2M0 = cell(1,obj.nbFiles) ;
            obj.dataM1M0 = cell(1,obj.nbFiles) ;
            obj.dataM2M0_interp = cell(1,obj.nbFiles) ;
            obj.dataM1M0_interp = cell(1,obj.nbFiles) ;
            obj.dataM0 = cell(1,obj.nbFiles) ;
            obj.dataM0_interp = cell(1,obj.nbFiles) ;
            obj.dataM1 = cell(1,obj.nbFiles) ;
            obj.dataM2 = cell(1,obj.nbFiles) ;
            obj.dataSH = cell(1,obj.nbFiles) ;
            obj.dataSH_interp = cell(1,obj.nbFiles) ;
            obj.k = 0;
            
         


            try
                checkPulsewaveParamsFromTxt(obj.directory);
                PW_params = Parameters(obj.directory);
            catch
                dir_path_txt = fullfile(path,'txt');
                delete(fullfile(dir_path_txt,'InputPulsewaveParams.txt'));
                checkPulsewaveParamsFromTxt(obj.directory);
                PW_params = Parameters(obj.directory);
            end


            obj.k = PW_params.k;
            obj.isdone_flatfield = 0;
            obj.ToolBoxmaster =  ToolBoxClass(obj.directory);

            for ii = 1 : obj.nbFiles

                %% Ref loading
                dir_path_avi = fullfile(obj.directory, 'avi');
                NameRefAviFile = strcat(obj.filenames{ii}, '_M0');
                RefAviFilePath = fullfile(dir_path_avi, NameRefAviFile);
                ext = '.avi';

                disp(['reading : ',RefAviFilePath]);
                V = VideoReader(fullfile(dir_path_avi,[NameRefAviFile,ext]));
                video = zeros(V.Height, V.Width, V.NumFrames);
                for n = 1 : V.NumFrames
                    video(:,:,n) = rgb2gray(read(V, n));
                end
                obj.reference{ii} = video;
                refvideosize = size(obj.reference{ii});

                %% File loading
                dir_path_raw = fullfile(obj.directory ,'raw');
                RawFilePath = fullfile(dir_path_raw,strcat(obj.filenames{ii}, '_DopplerRMS.raw'));
                NameRawFile = strcat(obj.filenames{ii}, '_DopplerRMS');
                ext = '.raw';

                % sqrt M2/M0 : DopplerRMS
                disp(['reading : ',fullfile(dir_path_raw,NameRawFile)]);
                fileID = fopen(RawFilePath);
                video = fread(fileID,'float32');
                fclose(fileID);
                obj.dataM2M0{ii} = reshape(video,refvideosize);
                % M1/M0 : DopplerAvg
                tmpname = NameRawFile;
                tmpname(end-2:end) = 'AVG';
                disp(['reading : ',fullfile(dir_path_raw,[tmpname,ext])]);
                fileID = fopen(fullfile(dir_path_raw,[tmpname,ext]));
                videoM1M0 = fread(fileID,'float32');
                fclose(fileID);
                obj.dataM1M0{ii} = reshape(videoM1M0,refvideosize);
                % Import Moment 0
                tmpname = strcat(NameRawFile(1:end-10), 'moment0');
                disp(['reading : ',fullfile(dir_path_raw,[tmpname,ext])]);
                fileID = fopen(fullfile(dir_path_raw,[tmpname,ext]));
                videoM0 = fread(fileID,'float32');
                fclose(fileID);
                obj.dataM0{ii} = reshape(videoM0,refvideosize);
                % Import Moment 1
                tmpname = strcat(NameRawFile(1:end-10), 'moment1');
                disp(['reading : ',fullfile(dir_path_raw,[tmpname,ext])]);
                fileID = fopen(fullfile(dir_path_raw,[tmpname,ext]));
                videoM1 = fread(fileID,'float32');
                fclose(fileID);
                obj.dataM1{ii} = reshape(videoM1,refvideosize);
                % Import Moment 2
                tmpname = strcat(NameRawFile(1:end-10), 'moment2');
                disp(['reading : ',fullfile(dir_path_raw,[tmpname,ext])]);
                fileID = fopen(fullfile(dir_path_raw,[tmpname,ext]));
                videoM2 = fread(fileID,'float32');
                fclose(fileID);
                obj.dataM2{ii} = reshape(videoM2,refvideosize);
%                 % Import SH
% 
%                 tmpname = strcat(NameRawFile(1:end-10), 'SH');
%                 disp(['reading : ',fullfile(dir_path_raw,[tmpname,ext])]);
%                 fileID = fopen(fullfile(dir_path_raw,[tmpname,ext]));
%                 videoSH = fread(fileID,'float32');
%                 fclose(fileID);
%                 obj.dataSH{ii} = reshape(videoSH,refvideosize(1),refvideosize(2),[]);
%   


            end
        end

        function obj = MomentNormalize(obj)
            for ii = 1 : obj.nbFiles
                avgM0 = mean(obj.dataM0{ii},[1 2]);
                avgRef = mean(obj.reference{ii},[1 2]);


                %donnée temporaire
                tmp_M2M0 = obj.dataM2M0{ii};
                tmp_M2 = obj.dataM2{ii};
                tmp_M1 = obj.dataM1{ii};
                tmp_M1M0 = obj.dataM1M0{ii};
                tmp_Ref = obj.reference{ii};

                for jj = 1 : size(obj.dataM2M0{ii}, 3)
                    tmp_M2M0(:,:,jj) = sqrt((tmp_M2(:,:,jj))./avgM0(:,:,jj));
                    tmp_M1M0(:,:,jj) = (tmp_M1(:,:,jj))./avgM0(:,:,jj);
                    tmp_Ref(:,:,jj) = tmp_Ref(:,:,jj)./avgRef(:,:,jj);

                end
                obj.dataM2M0{ii} = tmp_M2M0;
                obj.dataM1M0{ii} = tmp_M1M0;
                obj.reference_norm{ii} = tmp_Ref;
            end
        end

        function obj = MomentFlatField(obj)
            PW_params = Parameters(obj.directory);
            for ii = 1 : obj.nbFiles
                height = size(obj.dataM2M0_interp{1}, 1);
                width = size(obj.dataM2M0_interp{1}, 2);
                length = size(obj.dataM2M0_interp{1}, 3);
                num_frames = size(obj.reference{1}, 3);
                gwRatio = PW_params.flatField_gwRatio;
                flatField_border = PW_params.flatField_border;


                tmp_dataM0 = zeros(height, width, length);
                tmp_dataM2M0 = zeros(height, width, length);
                tmp_calc_data = obj.dataM2M0_interp{1};
                tmp_calc_data_M0 = obj.dataM0_interp{1};
                parfor i = 1 : num_frames % loop over frames
                    tmp_dataM0(:,:,i) =  flat_field_correction_old(tmp_calc_data_M0(:,:,i), gwRatio*height, flatField_border);
                    tmp_dataM2M0(:,:,i) = flat_field_correction_old(tmp_calc_data(:,:,i), gwRatio*height, flatField_border);

                end

                obj.dataM2M0_interp{1} = tmp_dataM2M0;
                obj.dataM0_interp{1} = tmp_dataM0;

            end
        end

        function obj = cropAllVideo(obj)
            %Crop a video (matrix dim 3)
            PW_params = Parameters(obj.directory);
            firstFrame = PW_params.videoStartFrameIndex;
            lastFrame = PW_params.videoEndFrameIndex;

            if firstFrame > 0 && firstFrame < size(obj.dataM2M0{1},3) && lastFrame > firstFrame && lastFrame <= size(obj.dataM2M0{1},3)


                obj.reference{1} = obj.reference{1}(:,:,firstFrame:lastFrame);
                obj.reference_norm{1} = obj.reference{1}(:,:,firstFrame:lastFrame);
                obj.dataM2M0{1} = obj.dataM2M0{1}(:,:,firstFrame:lastFrame);
                obj.dataM1M0{1} = obj.dataM1M0{1}(:,:,firstFrame:lastFrame);
                obj.dataM0{1} = obj.dataM0{1}(:,:,firstFrame:lastFrame);
                obj.dataM1{1} = obj.dataM1{1}(:,:,firstFrame:lastFrame);
                obj.dataM2{1} = obj.dataM2{1}(:,:,firstFrame:lastFrame);


                disp(['Data cube frame: ',num2str(firstFrame), '/', num2str(size(obj.dataM2M0{1},3)), ' to ',num2str(lastFrame), '/', num2str(size(obj.dataM2M0{1},3))])
            else
                disp('Wrong value for the first frame. Set as 1')
                disp('Wrong value for the last frame. Set as the end.')
                disp(['Data cube frame: 1/', num2str(size(obj.dataM2M0{1},3)), ' to ',num2str(size(obj.dataM2M0{1},3)), '/', num2str(size(obj.dataM2M0{1},3))])
            end
        end

        function obj = Interpolate(obj) %ref = TRUE indicates the object is the reference
            height = size(obj.reference{1}, 1);
            width = size(obj.reference{1}, 2);
            num_frames = size(obj.reference{1}, 3);
            k_interp = obj.k;
            height = (height-1)*(2^k_interp-1)+height;
            width = (width-1)*(2^k_interp-1)+width;


            tmp_ref = zeros(height,width, size(obj.dataM2M0{1}, 3));
            tmp_ref_norm = zeros(height,width, size(obj.dataM2M0{1}, 3));
            tmp_dataM0 = zeros(height, width, size(obj.dataM2M0{1}, 3));
            tmp_dataM2M0 = zeros(height, width, size(obj.dataM2M0{1}, 3));
            tmp_data_M1M0 = zeros(height, width, size(obj.dataM2M0{1}, 3));

            tmp_calc_ref = obj.reference{1};
            tmp_calc_ref_norm = obj.reference_norm{1};
            tmp_calc_data = obj.dataM2M0{1};
            tmp_calc_data_M0 = obj.dataM0{1};
            tmp_calc_data_M1M0 = obj.dataM1M0{1};

            parfor i = 1 : num_frames % loop over frames
                tmp_dataM0(:,:,i) = interp2(tmp_calc_data_M0(:,:,i), k_interp);
                tmp_dataM2M0(:,:,i) = interp2(tmp_calc_data(:,:,i), k_interp);
                tmp_data_M1M0(:,:,i) = interp2(tmp_calc_data_M1M0(:,:,i), k_interp);
                tmp_ref(:,:,i) = interp2( tmp_calc_ref(:,:,i), k_interp);
                tmp_ref_norm(:,:,i) = interp2(tmp_calc_ref_norm(:,:,i), k_interp);
            end

            obj.reference_interp{1} = tmp_ref;
            obj.reference_norm_interp{1} = tmp_ref_norm;
            obj.dataM2M0_interp{1} = tmp_dataM2M0;
            obj.dataM1M0_interp{1} = tmp_data_M1M0;
            obj.dataM0_interp{1} = tmp_dataM0;

%             if obj.flag_SH == 1
%                 tmp_data_SH = zeros(height, width, size(obj.dataSH{1}, 3));
%                 tmp_calc_data_SH = obj.dataSH{1};
%                 for i = 1 : size(obj.dataSH{1}, 3)
%                     tmp_data_SH(:,:,i) = interp2(tmp_calc_data_SH(:,:,i), k_interp);
%                 end
%                 obj.dataSH_interp{1} = tmp_data_SH;
% 
%             end



        end

        function onePulse(obj, Ninterp)
            PW_params = Parameters(obj.directory);
            ToolBox = obj.ToolBoxmaster;
           

            meanIm = squeeze(mean(obj.reference_interp{1}, 3));
            blurred_mask = imgaussfilt(double(meanIm),0.05*size(meanIm,1),'Padding',0);
            [ToolBox.y_barycentre,ToolBox.x_barycentre] = find(blurred_mask == max(blurred_mask,[],'all'));

            mkdir(ToolBox.PW_path_dir);
            mkdir(ToolBox.PW_path_png);
            mkdir(ToolBox.PW_path_eps);
            mkdir(ToolBox.PW_path_txt);
            mkdir(ToolBox.PW_path_avi);
            mkdir(ToolBox.PW_path_mp4);

            path_dir_txt = fullfile(obj.directory,'txt');
            path_file_txt_params = fullfile(path_dir_txt,'InputPulsewaveParams.txt');
            copyfile(path_file_txt_params,ToolBox.PW_path_txt );
            
            %saving times
            path_file_txt_exe_times = fullfile(path_dir_txt, 'ExecutionTimes.txt'));

            %% FlatField 

            if obj.flag_FlatField && (obj.isdone_flatfield == 0) 

                disp("Correcting data with flat field")
                for ii = 1 : obj.nbFiles
                height = size(obj.dataM2M0_interp{1}, 1);
                width = size(obj.dataM2M0_interp{1}, 2);
                length = size(obj.dataM2M0_interp{1}, 3);
                num_frames = size(obj.reference{1}, 3);
                gwRatio = PW_params.flatField_gwRatio;
                flatField_border = PW_params.flatField_border;


                tmp_dataM0 = zeros(height, width, length);
                tmp_dataM2M0 = zeros(height, width, length);
                tmp_calc_data = obj.dataM2M0_interp{1};
                tmp_calc_data_M0 = obj.dataM0_interp{1};
                parfor i = 1 : num_frames % loop over frames
                    tmp_dataM0(:,:,i) =  flat_field_correction_old(tmp_calc_data_M0(:,:,i), gwRatio*height, flatField_border);
                    tmp_dataM2M0(:,:,i) = flat_field_correction_old(tmp_calc_data(:,:,i), gwRatio*height, flatField_border);

                end

                obj.dataM2M0_interp{1} = tmp_dataM2M0;
                obj.dataM0_interp{1} = tmp_dataM0;

                end
                obj.isdone_flatfield = 1;
            else
                disp("no flatfield applied")
            end

            %% Creating Masks

            tic
          [maskArtery, maskVein, maskVessel,maskBackground,maskCRA,maskCRV,maskSectionArtery] = createMasks(obj.reference_norm_interp{1} ,obj.dataM1M0_interp{1}, obj.directory, ToolBox);
            disp('CreatMasks timing :')
            time = toc;
            disp(time)
            save_time(path_file_txt_exe_times, 'createMasks', time)

            tic
          % [maskArtery, maskVein, maskVessel,maskBackground,maskCRA,maskCRV,maskSectionArtery]= createMasksNew(obj.reference_interp{1} ,obj.dataM1M0_interp{1}, obj.directory, ToolBox);
            % [mask_artery, mask_vein, mask_vessel,mask_background,mask_CRA,maskCRV] = createMasksNew(obj.reference_interp{1} ,obj.dataM1M0_interp{1}, obj.directory, ToolBox);

            disp('CreatMasks New timing :')
            time = toc;
            disp(time)
            %save_time(path_file_txt_exe_times, 'createMasksNew', time)

            %% PulseWave Analysis 

            if obj.flag_PulseWave_analysis
                close all

                sys_index_list_cell = cell(obj.nbFiles) ;
                fullPulseWave_cell = cell(obj.nbFiles) ;
                for n = 1:obj.nbFiles
                    tic
                    % [sys_index_list_cell{n}, fullPulseWave_cell{n}] = find_systole_index(obj.reference_norm_interp{n}, obj.directory,maskArtery);
                    [sys_index_list_cell{n}, fullPulseWave_cell{n}] = find_systole_index(obj.reference_interp{n}, obj.directory,maskArtery);
                    disp('FindSystoleIndex timing :')
                    time = toc;
                    disp(time)
                    save_time(path_file_txt_exe_times, 'find_systole_index', time)

                    [v_RMS_one_cycle,v_RMS_all] = pulseAnalysis(Ninterp,obj.dataM2M0_interp{n},obj.dataM1M0_interp{n},sys_index_list_cell{n},maskArtery,maskVein,maskBackground ,ToolBox,obj.directory);
                    disp('PulseAnalysis timing :')
                    time = toc;
                    disp(time)
                    save_time(path_file_txt_exe_times, 'pulseAnalysis', time)

                    tic
                    velocity_map(maskArtery, maskVein, v_RMS_one_cycle, ToolBox); %FIXME Histo en trop, jsute pour faire la FLOWMAP ? 
                    disp('Velocity map timing :')
                    time = toc;
                    disp(time)
                    save_time(path_file_txt_exe_times, 'velocity_map', time)

                    tic
                    velocityHistogramm(v_RMS_all, maskArtery,maskVein ,ToolBox)
                    disp('Velocity Histogramm timing :')
                    toc


                    tic
                    ArterialResistivityIndex(v_RMS_one_cycle,obj.dataM2M0_interp{n}, maskArtery,  ToolBox);
                    disp('ArterialResistivityIndex timing :')
                    time = toc;
                    disp(time)
                    save_time(path_file_txt_exe_times, 'ArterialResistivityIndex', time)

                    %[v] = pulseAnalysisTest(Ninterp,obj.dataM2M0_interp{n},obj.dataM1M0_interp{n},sys_index_list_cell{n},maskArtery,maskVessel,maskVein,maskBackground ,ToolBox,obj.directory);
                    %                     detectElasticWave(datacube, maskArtery, maskCRA);
                    tic
                    try
                        flow_rate(maskArtery, maskVein, maskCRA,maskSectionArtery, v_RMS_all,obj.dataM0_interp{1}, ToolBox, obj.k,obj.directory);
                    catch
                        disp('FlowRate timing :')
                    
                    time = toc;
                    disp(time)
                    save_time(path_file_txt_exe_times, 'flow_rate', time)
                    %                     try
                    %                         flow_rate_old(maskArtery, maskVein, maskCRA, v_RMS, ToolBox, obj.k,obj.directory);
                    %                     catch
                    %                     end

                   

%                         tic
%                         spectrogram(maskArtery,maskBackground,  obj.dataSH_interp{n}, ToolBox);
%                         disp('Spectrogram timing :')
%                         time = toc;
%                         disp(time)
%                         save_time(path_file_txt_exe_times, 'spectrogram', time)
%                     else

                    end
                end
            end

            %% Spectrum Analysis

            if obj.flag_SH_analysis

                %% Import SH
                cubeSize = 256 ;
                cubeFreqLength = 32 ;
                tmpname  = strcat(ToolBox.main_foldername, '_SH');
                ext = '.raw';
                disp(['reading : ',fullfile(obj.directory,'raw',[tmpname,ext])]);
                fileID = fopen(fullfile(obj.directory,'raw',[tmpname,ext]));
                videoSH = fread(fileID,'float32');
                fclose(fileID);
                SH_cube=  reshape(videoSH,cubeSize,cubeSize,cubeFreqLength,[]);

                tic
                spectrum_analysis(maskArtery,maskBackground ,SH_cube ,ToolBox,obj.dataM0_interp{1});
                disp('Spectrum analysis :')
                time = toc;
                disp(time)
                save_time(path_file_txt_exe_times, 'spectrum_analysis', time)

                tic
                spectrogram_new(maskArtery,maskBackground ,SH_cube ,ToolBox);
                disp('Spectrogram timing :')
                time = toc;
                disp(time)
                save_time(path_file_txt_exe_times, 'spectrogram_new', time)
            end
            displaySuccessMsg(1);
            close all 
        end
    end
end