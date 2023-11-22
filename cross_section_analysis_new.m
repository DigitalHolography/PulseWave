function [avg_blood_volume_rate,std_blood_volume_rate, cross_section_area, avg_blood_velocity, cross_section_mask,total_avg_blood_volume_rate,total_std_blood_volume_rate] = cross_section_analysis_new(locs, width, mask, v_RMS, slice_half_thickness, k,ToolBox, path,fig)
% validate_cross_section
%   Detailed explanation goes here FIXME
PW_params = Parameters(path);
subImg_cell = cell(size(locs,1));
subVideo_cell = cell(size(locs,1));

[M,N,T_max] = size(v_RMS);
width_cross_section = zeros(size(locs,1),1);
cross_section_area = zeros(size(locs,1),1);
avg_blood_velocity = zeros(size(locs,1),1);
avg_blood_volume_rate = zeros(size(locs,1),1);
std_blood_velocity = zeros(size(locs,1),1);
std_blood_volume_rate = zeros(size(locs,1),1);
cross_section_mask = zeros(size(mask));
mask_sections = zeros(M,N,size(locs,1));
total_avg_blood_volume_rate = zeros(T_max,1);
total_std_blood_volume_rate = zeros(T_max,1);

% %% VARIABLES FOR VELOCITY PROFILE VIDEO

tilt_angle_list = zeros(1,length(locs));


img_v_artery = squeeze(mean(v_RMS,3)).* mask;
v_RMS_masked = v_RMS.* mask;
for ii = 1:size(locs,1)
    if width(ii)>2
    subImgHW = round(width(ii)*PW_params.cropSection_scaleFactorWidth);
    %FIXME bords d IMG, 

    xRange = round(-subImgHW/2) + locs(ii,2) : round(subImgHW/2) + locs(ii,2);
    yRange = round(-subImgHW/2) + locs(ii,1) : round(subImgHW/2) + locs(ii,1);
    subImg = img_v_artery(yRange ,xRange);

    %make disk mask
    %FIXME img anamorphique
    subImg = cropCircle(subImg);


    theta = linspace(0,180,181);
    projx = zeros(size(subImg,1),length(theta));
    projy = zeros(size(subImg,2),length(theta));
    Video_subIm_rotate = zeros(size(subImg,1),size(subImg,2),length(theta));
    for tt = 1:length(theta)
        tmpImg = imrotate(subImg,theta(tt),'bilinear','crop');
        Video_subIm_rotate(:,:,tt) = tmpImg;
        projx(:,tt) = squeeze(sum(tmpImg,1));
        projy(:,tt) = squeeze(sum(tmpImg,2));
    end
    figure(3001)
    imagesc(projx)

    figure(3002)
    imagesc(projy)
    % avi
    

 

    % [max_projx,tilt_idx] = max(projx(:),[],'all','linear');
    % [row,col] = ind2sub(size(squeeze(projx)),tilt_idx);
    % tilt_angle{ii} = col;
%     [~,tilt_angle] = find(projx == max_projx); %x_max angle de rotation pour une coupe normale au vaisseau
    
    projx_bin = (projx == 0);
    list_x = squeeze(sum(projx_bin, 1));
    [~, idc] = max(list_x);
    tilt_angle_list(ii) = idc(1);
    subImg = imrotate(subImg,tilt_angle_list(ii),'bilinear','crop');
    subImg_cell{ii} = subImg;
    subVideo = v_RMS_masked(yRange, xRange, :);
    for tt = 1:size(subVideo,3)
        subVideo(:,:,tt) = imrotate(cropCircle(subVideo(:,:,tt)),tilt_angle_list(ii),'bilinear','crop');
    end
    subVideo_cell{ii} = subVideo;
    section_cut = projx(:,tilt_angle_list(ii));
    for zz = 1:length(section_cut)
        if section_cut(zz) < 0 
            section_cut(zz) = 0;
        end 
    end 

    tmp_section = (section_cut./max(section_cut))*size(section_cut,1);


    figure(1000+ii)
    xAx = linspace(0,size(projx,1),size(projx,1));
    xAy = linspace(0,size(projx,2),size(projx,2));
    imagesc(xAy,xAx,projx)
    colormap("gray")
    axis image 
    hold on 
    x = [ tilt_angle_list(ii)  tilt_angle_list(ii)];
    y = [0 size(projx,1)];
    line(x,y,'Color','red','LineStyle',':','LineWidth',3)
    axis off;
    hold off;
    set(gca,'PlotBoxAspectRatio',  [1,1.618,1]);
    f = getframe(gca);              %# Capture the current window

    imwrite(f.cdata, fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,['_Artery_Section_proj_' num2str(ii) '.png'])));

   % Video_subIm_rotate = circshift(Video_subIm_rotate,[0 0 -tilt_angle_list(ii)]);
    w = VideoWriter(fullfile(ToolBox.PW_path_avi,strcat(ToolBox.main_foldername,['_Artery_Section_' num2str(ii) '.avi'])));
    tmp_video = mat2gray( Video_subIm_rotate);
    open(w)
    for j = 1:size( Video_subIm_rotate,3)
        writeVideo(w,tmp_video(:,:,j)) ;
    end
    close(w);



    %[ ~, ~, tmp, ~] = findpeaks(section_cut,1:size(subImg,1), 'MinPeakWidth', round(PW_params.cropSection_scaleFactorSize*size(mask,1)));
    tmp = nnz(section_cut);

    % if tmp contains more than 1 element, we select the first one
    % if tmp is empty, we select 0 as width because then it's a 1-2 pixel noise peak
    if isempty(tmp) 
        width_cross_section(ii) = 0;
    else
        width_cross_section(ii) = tmp(1);
    end


    figure(fig+ii)
    xAx = linspace(0,size(section_cut,1),size(subImg,1));
    imagesc(xAx,xAx,subImg)
    colormap("gray")
    axis image 
    hold on 
    p = plot(xAx, tmp_section);
    p.LineWidth = 2;
    p.Color = 'red';
    p.LineStyle = ':';
    set(gca,'PlotBoxAspectRatio',  [1,1,1]);
    x = [round(size(subImg,1)/2)-round(width_cross_section(ii)/2) round(size(subImg,1)/2)+round(width_cross_section(ii)/2)];
    y = [round(size(subImg,1)/2) round(size(subImg,1)/2)];
    line(x,y,'Color','red','LineWidth',3)
    axis off;
    f = getframe(gca);              %# Capture the current 
    
    %bords blancs
    imwrite(f.cdata, fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,['_Artery_Section_' num2str(fig+ii) '.png'])));


    mask_slice_subImg = false(size(subImg,1),size(subImg,2));
    slice_center = round(size(subImg,1)/2);

    slice_half_thickness_tmp = min(slice_half_thickness,floor(subImgHW/2));
    mask_slice_subImg((slice_center - slice_half_thickness_tmp):(slice_center + slice_half_thickness_tmp),:) = true;
    mask_slice_subImg = imrotate(double(mask_slice_subImg),-tilt_angle_list(ii),'bilinear','crop');
    mask_slice_subImg = mask_slice_subImg>PW_params.cropSection_maskThreshold;


    %% Average the blood flow calculation over a circle before dividing by the section

    mask_current_slice = false(size(mask));
    mask_current_slice(1:size(mask_slice_subImg,1), 1:size(mask_slice_subImg,2)) = mask_slice_subImg;

    shift_x = locs(ii,2) - round(size(mask_slice_subImg,1)/2);
    shift_y = locs(ii,1) - round(size(mask_slice_subImg,2)/2);

    mask_current_slice = circshift(mask_current_slice, [shift_y shift_x]);
    mask_current_slice = mask_current_slice.*mask;
    mask_current_slice_inverse = ~mask_current_slice;
    mask_current_slice = double(mask_current_slice);
    mask_current_slice_inverse = double(mask_current_slice_inverse);
    cross_section_mask = cross_section_mask + mask_current_slice;
    mask_sections(:,:,ii) = mask_current_slice;

    end
    cross_section_area(ii) = pi*((width_cross_section(ii)/2)*(PW_params.cropSection_pixelSize/2^k))^2; % /2 because radius=d/2 - 0.0102/2^k mm = size pixel with k coef interpolation
end

local_velocity_pulses = zeros(size(locs,1),T_max);
for tt = 1:T_max
    current_frame = v_RMS(:,:,tt);
    all_velocity = zeros(M,N);
    Total_cross_section =0;
    for ii = 1:size(locs,1)

        %FIXME mean(v_RMS,3)
        tmp = current_frame.*mask_sections(:,:,ii);
        %tmp_velocity = zeros(1,size(nnz(tmp(:))));
        xRange = round(-subImgHW/2) + locs(ii,2) : round(subImgHW/2) + locs(ii,2);
        yRange = round(-subImgHW/2) + locs(ii,1) : round(subImgHW/2) + locs(ii,1);
        subFrame = tmp(yRange ,xRange);
        subFrame = cropCircle(subFrame);
        subFrame = imrotate(subFrame,tilt_angle_list(ii),'bilinear','crop');
        avg_profil = mean(subFrame,1);
        for ll = 1:size(subFrame,1)
            subFrame(ll,:) = subFrame(ll,:)-avg_profil;

        end


        %FIXME calcul std avg avec des v = 0
        %avg_blood_velocity(ii,tt) = sum(tmp(:))/nnz(tmp(:)); 
        avg_blood_velocity(ii,tt) = mean(tmp(tmp~=0)); 
        %std_blood_velocity(ii,tt) = std(tmp(tmp~=0));
        std_blood_velocity(ii,tt) = std(subFrame(subFrame~=0));
        avg_blood_volume_rate(ii,tt) = avg_blood_velocity(ii,tt)*cross_section_area(ii)*60; % microL/min
        std_blood_volume_rate(ii,tt) = std_blood_velocity(ii,tt)*cross_section_area(ii)*60; % microL/min
        all_velocity(:,:) = all_velocity(:,:) + current_frame.*mask_sections(:,:,ii);
        Total_cross_section = Total_cross_section+cross_section_area(ii);
        

        %     figure(101)
        %     plot(plot_values);
        %     findpeaks(plot_values,size(plot_values, 2),'MinPeakProminence',param_peak);
        %     % text(locs,pks,num2str((1:numel(pks))'))
        %     text(locs,pks,string(round(avg_blood_rate,3)))
        %     title("Peaks of luminosity")
        %     pbaspect([1.618 1 1]);

        %print(['-f' num2str(70+ii)],'-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,['_Artery_Section_' num2str(ii) '.png']))) ;
        %print(['-f' num2str(1000+ii)],'-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,['_Proj_Artery_Section_' num2str(ii) '.png']))) ;



    end

  
    if ~isempty(avg_blood_volume_rate)
        total_avg_blood_volume_rate(tt) = sum(avg_blood_volume_rate(:,tt));
        total_std_blood_volume_rate(tt) = mean(std_blood_volume_rate(:,tt));
%         total_std_blood_volume_rate(tt) = (std(all_velocity(all_velocity~=0))*Total_cross_section*60);% microL/min; %FIXME calculer le vrai std 
    else
        total_avg_blood_volume_rate(tt)= 0;
        total_std_blood_volume_rate(tt) = 0;
    end
end % ii


try
viscosity_new(subImg_cell , subVideo_cell,  ToolBox);
%viscosity_video = viscosity(subImg_cell, subVideo_cell, tilt_angle_list, ToolBox.PW_path_dir, ToolBox.main_foldername);
catch




end