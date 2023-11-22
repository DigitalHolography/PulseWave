function [] = flow_rate(maskArtery, maskVein, maskCRA,maskSection, v_RMS,dataM0, ToolBox, k,path)

PW_params = Parameters(path);

%maskArtery = imdilate(maskArtery,strel('disk',5));
%FIXME function velocity map

v_RMS_AVG = mean(v_RMS,3);

%% Find the locations of the sections in artery

maskSectionArtery = maskSection.*maskArtery;

figure(111)
imagesc(maskSectionArtery.*v_RMS_AVG);

maskSectionArtery = bwlabel(maskSectionArtery);

figure(222)
imagesc(maskSectionArtery)

nb_sections_artery = max(maskSectionArtery,[],'all');
masksSectionsArtery = zeros(size(maskArtery,1),size(maskArtery,2),nb_sections_artery);

for ii = 1:nb_sections_artery
    masksSectionsArtery(:,:,ii) = (maskSectionArtery == ii);
end 

SubImg_locs_artery = zeros(nb_sections_artery,2);
SubImg_width_artery = zeros(nb_sections_artery,1);

for ii = 1:nb_sections_artery
    [row,col] = find(masksSectionsArtery(:,:,ii));
    SubImg_locs_artery(ii,1) = round(mean(row));
    SubImg_locs_artery(ii,2) = round(mean(col));
    SubImg_width_artery(ii) = 0.01*size(maskArtery,1);
end

%% Find the locations of the sections in artery

maskSectionVein = maskSection.*maskVein;

figure(1111)
imagesc(maskSectionVein.*v_RMS_AVG);

maskSectionVein = bwlabel(maskSectionVein);

figure(2222)
imagesc(maskSectionVein)

nb_sections_vein = max(maskSectionVein,[],'all');
masksSectionsVein = zeros(size(maskVein,1),size(maskVein,2),nb_sections_vein);

for ii = 1:nb_sections_vein
    masksSectionsVein(:,:,ii) = (maskSectionVein == ii);
end 

SubImg_locs_vein = zeros(nb_sections_vein,2);
SubImg_width_vein = zeros(nb_sections_vein,1);

for ii = 1:nb_sections_vein
    [row,col] = find(masksSectionsVein(:,:,ii));
    SubImg_locs_vein(ii,1) = round(mean(row));
    SubImg_locs_vein(ii,2) = round(mean(col));
    SubImg_width_vein(ii) = 0.01*size(maskVein,1);
end

%% Compute blood volume rate
mask_artery = maskArtery;
mask_vein = maskVein;
[avg_blood_volume_rate_vein,std_blood_volume_rate_vein, cross_section_area_vein, avg_blood_velocity_vein, cross_section_mask_vein,total_avg_blood_volume_rate_vein,total_std_blood_volume_rate_artery] = cross_section_analysis_new(SubImg_locs_vein, SubImg_width_vein, mask_vein, v_RMS, PW_params.flowRate_sliceHalfThickness, k,ToolBox,path,80);
[avg_blood_volume_rate_artery,std_blood_volume_rate_artery, cross_section_area_artery, avg_blood_velocity_artery, cross_section_mask_artery,total_avg_blood_volume_rate_artery,total_std_blood_volume_rate_vein] = cross_section_analysis_new(SubImg_locs_artery, SubImg_width_artery, mask_artery, v_RMS, PW_params.flowRate_sliceHalfThickness, k,ToolBox,path,70);




labels_arteries = {};
labels_veins = {};
strXlabel = 'Time(s)' ;%createXlabelTime(1);
strYlabel = 'Velocity (mm.s-1)';

figure(57)
hold on
for ii=1:size(avg_blood_volume_rate_artery,1)
   plot(smoothdata(avg_blood_velocity_artery(ii,:),'lowess'))
   labels_arteries{end+1} = strcat('A',num2str(ii));
end
title('Blood velocity in artery sections');
legend(labels_arteries);
fontsize(gca,12,"points") ;
xlabel(strXlabel,'FontSize',14) ;
ylabel(strYlabel,'FontSize',14) ;
pbaspect([1.618 1 1]) ;
set(gca, 'LineWidth', 2);
axis tight;

figure(58)
hold on
for ii=1:size(avg_blood_volume_rate_vein,1)
   plot(smoothdata(avg_blood_velocity_vein(ii,:),'lowess'))
   labels_veins{end+1} = strcat('V',num2str(ii));
end
title('Blood velocity in vein sections');
legend(labels_veins);
fontsize(gca,12,"points") ;
xlabel(strXlabel,'FontSize',14) ;
ylabel(strYlabel,'FontSize',14) ;
pbaspect([1.618 1 1]) ;
set(gca, 'LineWidth', 2);
axis tight;


dataM0 = mat2gray(dataM0);
maskOnes = ones(size(maskArtery,1), size(maskArtery,2),size(dataM0,3));
fullTime = linspace(0,size(v_RMS,3)*ToolBox.stride/ToolBox.fs/1000,size(v_RMS,3));
mean_M0 = mean(dataM0,3);
ratio_etiquette = 1.2;
%% Saving Vein and artery data in txt

for ii=1:size(avg_blood_volume_rate_artery,1)
    plot2txt(fullTime,avg_blood_volume_rate_artery(ii,:), strcat('blood_volume_rate_artery_A',num2str(ii)), ToolBox)
    plot2txt(fullTime,std_blood_volume_rate_artery(ii,:),  strcat('blood_volume_rate_artery_std_A',num2str(ii)), ToolBox)
    plot2txt(fullTime,avg_blood_velocity_artery(ii,:),  strcat('avg_velocity_artery_A',num2str(ii)), ToolBox)
    
end
for ii=1:size(avg_blood_volume_rate_vein,1)
    plot2txt(fullTime,avg_blood_volume_rate_vein(ii,:), strcat('blood_volume_rate_vein_V',num2str(ii)), ToolBox)
    plot2txt(fullTime,std_blood_volume_rate_vein(ii,:),  strcat('blood_volume_rate_vein_std_V',num2str(ii)), ToolBox)
    plot2txt(fullTime,avg_blood_velocity_vein(ii,:),  strcat('avg_velocity_vein_V',num2str(ii)), ToolBox)
   
end

%% Vein and artery numerotation

maskRGB = ones(size(maskArtery,1), size(maskArtery,2),3);


maskRGB(:,:,3) = mat2gray(mean_M0).*~cross_section_mask_vein + maskRGB(:,:,3).*cross_section_mask_vein;
maskRGB(:,:,2) = mat2gray(mean_M0).*~cross_section_mask_artery;
maskRGB(:,:,1) = mat2gray(mean_M0).*~cross_section_mask_artery+maskRGB(:,:,1).*cross_section_mask_artery;



figure(652)
imshow(maskRGB);

x_center = ToolBox.x_barycentre;
y_center = ToolBox.y_barycentre;
for ii=1:size(avg_blood_volume_rate_artery,1)
    new_x = x_center + ratio_etiquette*(SubImg_locs_artery(ii,2)-x_center);
    new_y = y_center + ratio_etiquette*(SubImg_locs_artery(ii,1)-y_center);
    text(new_x, new_y, strcat("A", num2str(ii)), "FontWeight", "bold","FontSize", 15,   "Color", "white", "BackgroundColor", "black");
end
for ii=1:size(avg_blood_volume_rate_vein,1)
    new_x = x_center + ratio_etiquette*(SubImg_locs_vein(ii,2)-x_center);
    new_y = y_center + ratio_etiquette*(SubImg_locs_vein(ii,1)-y_center);
    text(new_x, new_y, strcat("V", num2str(ii)), "FontWeight", "bold","FontSize", 15,   "Color", "white", "BackgroundColor", "black");
end
drawnow
ax = gca;
ax.Units = 'pixels';
marg = 30;
pos = ax.Position;
rect = [-marg, -marg, pos(3)+2*marg, pos(4)+2*marg];
F_numerotation = getframe(ax,rect);



%% Volume Rate in Arteries 
Color_std = [0.7 0.7 0.7];

maskRGB_artery = ones(size(maskArtery,1), size(maskArtery,2),3);
volume_rate_video_artery = zeros(size(v_RMS,1),size(v_RMS,2),3,size(v_RMS,3));

mean_volume_rate_artery = ones(length(fullTime),1)*mean(total_avg_blood_volume_rate_artery);

maskRGB_artery(:,:,3) = mat2gray(mean_M0).*~cross_section_mask_artery;
maskRGB_artery(:,:,2) = mat2gray(mean_M0).*~cross_section_mask_artery;
maskRGB_artery(:,:,1) = mat2gray(mean_M0).*~cross_section_mask_artery+maskRGB_artery(:,:,1).*cross_section_mask_artery;


figure(111222)
imagesc(mean(v_RMS,3).*cross_section_mask_artery)
colormap("gray")
Total_blood_flow_rate_artery = sum(avg_blood_volume_rate_artery);
disp('Total_blood_flow_rate')
disp(Total_blood_flow_rate_artery)
title('Artery Sections' );
axis image
axis off
for tt = 1:size(v_RMS,3)

    hue = 0*(maskOnes(:,:,tt).*cross_section_mask_artery);
    sat = maskOnes(:,:,tt).*cross_section_mask_artery;
    val = dataM0(:,:,tt).*~cross_section_mask_artery+maskOnes(:,:,tt).*cross_section_mask_artery;

    tmp_frame = hsv2rgb(hue,sat,val);
    imshow(tmp_frame);
    x_center = ToolBox.x_barycentre;
    y_center = ToolBox.y_barycentre;

    for ii=1:size(avg_blood_volume_rate_artery,1)
        new_x = x_center + ratio_etiquette*(SubImg_locs_artery(ii,2)-x_center);
        new_y = y_center + ratio_etiquette*(SubImg_locs_artery(ii,1)-y_center);
        if round(avg_blood_volume_rate_artery(ii,tt),1) > 0
            text(new_x, new_y, string(round(avg_blood_volume_rate_artery(ii,tt),1)), "FontWeight", "bold","FontSize", 15,   "Color", "white", "BackgroundColor", "black");
        else
            text(new_x, new_y, string(0), "FontWeight", "bold","FontSize", 15,   "Color", "white", "BackgroundColor", "black");
        end
    end
    tmp_bvra = round(total_avg_blood_volume_rate_artery(tt));
    tmp_bvra = (tmp_bvra>0)*tmp_bvra;
    title(['Total blood volume rate : ' num2str(tmp_bvra) ' µL/min (arteries) - ']);
    set(gca,'FontSize', 16)
    drawnow
    ax = gca;
    ax.Units = 'pixels';
    marg = 30;
    pos = ax.Position;
    rect = [-marg, -marg, pos(3)+2*marg, pos(4)+2*marg];
    Total_blood_volume_rate = getframe(ax,rect);
    volume_rate_video_artery(:,:,:,tt) = imresize(Total_blood_volume_rate.cdata,[size(volume_rate_video_artery,1) size(volume_rate_video_artery,2)]);
end

w = VideoWriter(fullfile(ToolBox.PW_path_avi,strcat(ToolBox.main_foldername,'_ volume_rate_video_artery.avi')));
tmp = mat2gray( volume_rate_video_artery);
open(w)
for j = 1:size( volume_rate_video_artery,4)
    writeVideo(w,tmp(:,:,:,j)) ;  
end
close(w);

figure(121)
imshow(maskRGB_artery);
x_center = ToolBox.x_barycentre;
y_center = ToolBox.y_barycentre;
for ii=1:size(avg_blood_volume_rate_artery,1)
    new_x = x_center + ratio_etiquette*(SubImg_locs_artery(ii,2)-x_center);
    new_y = y_center + ratio_etiquette*(SubImg_locs_artery(ii,1)-y_center);
    text(new_x, new_y, string(round(mean(avg_blood_volume_rate_artery(ii,:),2))), "FontWeight", "bold","FontSize", 15,   "Color", "white", "BackgroundColor", "black");
end
title(['Total blood volume rate : ' num2str(round(mean(total_avg_blood_volume_rate_artery))) ' µL/min (arteries) - ']);
drawnow
ax = gca;
ax.Units = 'pixels';
marg = 30;
pos = ax.Position;
rect = [-marg, -marg, pos(3)+2*marg, pos(4)+2*marg];
F_Total_blood_volume_rate_artery = getframe(ax,rect);




plot_vol_rate_artery = figure(125);
plot_vol_rate_artery.Position(3:4) = [600 275];
%plot(fullTime,total_avg_blood_volume_rate_artery,'-k','LineWidth',2);
curve1 = total_avg_blood_volume_rate_artery+0.5*total_std_blood_volume_rate_artery;
curve2 = total_avg_blood_volume_rate_artery-0.5*total_std_blood_volume_rate_artery;
fullTime2 = [fullTime, fliplr(fullTime)];
inBetween = [curve1', fliplr(curve2')];
fill(fullTime2, inBetween, Color_std);
hold on;
plot(fullTime,curve1,"Color", Color_std, 'LineWidth', 2);
plot(fullTime, curve2, "Color",Color_std, 'LineWidth', 2);
plot(fullTime,total_avg_blood_volume_rate_artery,'-k','LineWidth',1);
axis tight ;
hold off


ax_vol_rate_artery = axis;
ax_vol_rate_artery(3) = 0 ;
axis([ax_vol_rate_artery(1) ax_vol_rate_artery(2) ax_vol_rate_artery(3) ax_vol_rate_artery(4)]);
yticks(sort([0 20 40 round(mean_volume_rate_artery(1))]))

Plot_volume_rate_artery = getframe(gcf);
Plot_volume_rate_video_artery = zeros(size(Plot_volume_rate_artery.cdata,1),size(Plot_volume_rate_artery.cdata,2),3,size(v_RMS,3));

for tt = 1:length(fullTime)
    figure(125);
    curve1 = total_avg_blood_volume_rate_artery(1:tt)+0.5*total_std_blood_volume_rate_artery(1:tt);
    curve2 = total_avg_blood_volume_rate_artery(1:tt)-0.5*total_std_blood_volume_rate_artery(1:tt);
    tmp_fullTime = [fullTime(1:tt), fliplr(fullTime(1:tt))];
    inBetween = [curve1', fliplr(curve2')];
    fill( tmp_fullTime, inBetween, Color_std);
    hold on;
    plot(fullTime(1:tt),curve1,"Color", Color_std, 'LineWidth', 2);
    plot(fullTime(1:tt), curve2, "Color",Color_std, 'LineWidth', 2);
    plot(fullTime(1:tt),total_avg_blood_volume_rate_artery(1:tt),'-k','LineWidth',2);
    plot(fullTime(1:tt),mean_volume_rate_artery(1:tt),':k','LineWidth',2);

%     plot(fullTime(1:tt),total_avg_blood_volume_rate_artery(1:tt),'-k',...
%         fullTime(1:tt),total_avg_blood_volume_rate_artery(1:tt)+total_std_blood_volume_rate_artery(1:tt),'-k',...
%         fullTime(1:tt),total_avg_blood_volume_rate_artery(1:tt)-total_std_blood_volume_rate_artery(1:tt),'-k',...
%         fullTime(1:tt),mean_volume_rate_artery(1:tt),':k','LineWidth',2);

    ylabel('Blood volume rate (µL/min)')
    xlabel('Time (s)')
    title("Total blood volume rate in arteries")
    fontsize(gca,14,"points") ;
    yticks(sort([0 20 40 round(mean_volume_rate_artery(1))]))
    set(gca,'PlotBoxAspectRatio',[2.5 1 1])
    axis([ax_vol_rate_artery(1) ax_vol_rate_artery(2) ax_vol_rate_artery(3) ax_vol_rate_artery(4)]);
    Plot_volume_rate_artery = getframe(gcf);
    Plot_volume_rate_video_artery(:,:,:,tt) = Plot_volume_rate_artery.cdata;

end

plot2txt(fullTime(1:tt),total_avg_blood_volume_rate_artery(1:tt), 'TotalBloodVolumeRateArtery', ToolBox)
plot2txt(fullTime(1:tt),total_std_blood_volume_rate_artery(1:tt), 'TotalBloodVolumeRateArtery', ToolBox)

w = VideoWriter(fullfile(ToolBox.PW_path_avi,strcat(ToolBox.main_foldername,'_plot_volume_rate_video_artery.avi')));
tmp = mat2gray( Plot_volume_rate_video_artery);
open(w)
for j = 1:size( Plot_volume_rate_video_artery,4)
    writeVideo(w,tmp(:,:,:,j)) ;  
end
close(w);

%% Volume Rate in Veins 

maskRGB_vein = ones(size(maskVein,1), size(maskVein,2),3);
volume_rate_video_vein = zeros(size(v_RMS,1),size(v_RMS,2),3,size(v_RMS,3));

mean_volume_rate_vein = ones(length(fullTime),1)*mean(total_avg_blood_volume_rate_vein);

maskRGB_vein(:,:,1) = mat2gray(mean_M0).*~cross_section_mask_vein;
maskRGB_vein(:,:,2) = mat2gray(mean_M0).*~cross_section_mask_vein;
maskRGB_vein(:,:,3) = mat2gray(mean_M0).*~cross_section_mask_vein + maskRGB_vein(:,:,3).*cross_section_mask_vein;


figure(111223)
imagesc(mean(v_RMS,3).*cross_section_mask_vein)
colormap("gray")
Total_blood_flow_rate_vein = sum(avg_blood_volume_rate_vein);
disp('Total_blood_flow_rate')
disp(Total_blood_flow_rate_vein)
title('Veins Sections' );
axis image
axis off
for tt = 1:size(v_RMS,3)

    hue = 0.6*(maskOnes(:,:,tt).*cross_section_mask_vein);
    sat = maskOnes(:,:,tt).*cross_section_mask_vein;
    val = dataM0(:,:,tt).*~cross_section_mask_vein + maskOnes(:,:,tt).*cross_section_mask_vein;

    tmp_frame = hsv2rgb(hue,sat,val);
    imshow(tmp_frame);
    x_center = ToolBox.x_barycentre;
    y_center = ToolBox.y_barycentre;

    for ii=1:size(avg_blood_volume_rate_vein,1)
        new_x = x_center + ratio_etiquette*(SubImg_locs_vein(ii,2)-x_center);
        new_y = y_center + ratio_etiquette*(SubImg_locs_vein(ii,1)-y_center);
        if round(avg_blood_volume_rate_vein(ii,tt),1) > 0
            text(new_x, new_y, string(round(avg_blood_volume_rate_vein(ii,tt),1)), "FontWeight", "bold","FontSize", 15,   "Color", "white", "BackgroundColor", "black");
        else
            text(new_x, new_y, string(0), "FontWeight", "bold","FontSize", 15,   "Color", "white", "BackgroundColor", "black");
        end
    end
    tmp_bvrv = round(total_avg_blood_volume_rate_vein(tt));
    tmp_bvrv = (tmp_bvrv>0)*tmp_bvrv;
    title(['Total blood volume rate : ' num2str(tmp_bvrv) ' µL/min (veins) - ']);
    set(gca,'FontSize', 16)
    drawnow
    ax = gca;
    ax.Units = 'pixels';

     drawnow
    ax = gca;
    ax.Units = 'pixels';
    marg = 30;
    pos = ax.Position;
    rect = [-marg, -marg, pos(3)+2*marg, pos(4)+2*marg];
    Total_blood_volume_rate = getframe(ax,rect);
    volume_rate_video_vein(:,:,:,tt) = imresize(Total_blood_volume_rate.cdata,[size(volume_rate_video_vein,1) size(volume_rate_video_vein,2)]);
end

w = VideoWriter(fullfile(ToolBox.PW_path_avi,strcat(ToolBox.main_foldername,'_ volume_rate_video_vein.avi')));
tmp = mat2gray( volume_rate_video_vein);
open(w)
for j = 1:size( volume_rate_video_vein,4)
    writeVideo(w,tmp(:,:,:,j)) ;  
end
close(w);

figure(122) 
imshow(maskRGB_vein);
x_center = ToolBox.x_barycentre;
y_center = ToolBox.y_barycentre;
for ii=1:size(avg_blood_volume_rate_vein,1)
    new_x = x_center + ratio_etiquette*(SubImg_locs_vein(ii,2)-x_center);
    new_y = y_center + ratio_etiquette*(SubImg_locs_vein(ii,1)-y_center);
    text(new_x, new_y, string(round(mean(avg_blood_volume_rate_vein(ii,:),2))), "FontWeight", "bold","FontSize", 15,   "Color", "white", "BackgroundColor", "black");
end
title(['Total blood volume rate : ' num2str(round(mean(total_avg_blood_volume_rate_vein))) ' µL/min (veins) - ']);
drawnow
ax = gca;
ax.Units = 'pixels';
marg = 30;
pos = ax.Position;
rect = [-marg, -marg, pos(3)+2*marg, pos(4)+2*marg];
F_Total_blood_volume_rate_vein = getframe(ax,rect);





plot_vol_rate_vein = figure(126);
plot_vol_rate_vein.Position(3:4) = [600 275];
%plot(fullTime,total_avg_blood_volume_rate_vein,'-k','LineWidth',2);
curve1 = total_avg_blood_volume_rate_vein+0.5*total_std_blood_volume_rate_vein;
curve2 = total_avg_blood_volume_rate_vein-0.5*total_std_blood_volume_rate_vein;
fullTime2 = [fullTime, fliplr(fullTime)];
inBetween = [curve1', fliplr(curve2')];
fill(fullTime2, inBetween, Color_std);
hold on;
plot(fullTime,curve1,"Color",Color_std, 'LineWidth', 2);
plot(fullTime, curve2, "Color",Color_std, 'LineWidth', 2);
plot(fullTime,total_avg_blood_volume_rate_vein,'-k','LineWidth',1);
axis tight ;
hold off

axis tight ;
axis([ax_vol_rate_artery(1) ax_vol_rate_artery(2) ax_vol_rate_artery(3) ax_vol_rate_artery(4)]);
Plot_volume_rate_vein = getframe(gcf);
Plot_volume_rate_video_vein = zeros(size(Plot_volume_rate_vein.cdata,1),size(Plot_volume_rate_vein.cdata,2),3,size(v_RMS,3));
% yticks(sort([0 20 40 round(mean_volume_rate_vein(1))]))

for tt = 1:length(fullTime)
    figure(126); 
    curve1 = total_avg_blood_volume_rate_vein(1:tt)+0.5*total_std_blood_volume_rate_vein(1:tt);
    curve2 = total_avg_blood_volume_rate_vein(1:tt)-0.5*total_std_blood_volume_rate_vein(1:tt);
    tmp_fullTime = [fullTime(1:tt), fliplr(fullTime(1:tt))];
    inBetween = [curve1', fliplr(curve2')];
    fill( tmp_fullTime, inBetween, Color_std);
    hold on;
    plot(fullTime(1:tt),curve1,"Color", Color_std, 'LineWidth', 2);
    plot(fullTime(1:tt), curve2, "Color",Color_std, 'LineWidth', 2);
    plot(fullTime(1:tt),total_avg_blood_volume_rate_vein(1:tt),'-k','LineWidth',2);
    plot(fullTime(1:tt),mean_volume_rate_vein(1:tt),':k','LineWidth',2);
%     plot(fullTime(1:tt),total_avg_blood_volume_rate_vein(1:tt),'-k',...
%         fullTime(1:tt),mean_volume_rate_vein(1:tt),':k','LineWidth',2);
    ylabel('Blood volume rate (µL/min)')
    xlabel('Time (s)')
    title("Total blood volume rate in veins")
    fontsize(gca,14,"points") ;
    yticks(sort([0 20 40 round(mean_volume_rate_vein(1))]))
    set(gca,'PlotBoxAspectRatio',[2.5 1 1])
    axis([ax_vol_rate_artery(1) ax_vol_rate_artery(2) ax_vol_rate_artery(3) ax_vol_rate_artery(4)]);
    Plot_volume_rate_vein = getframe(gcf);
    Plot_volume_rate_video_vein(:,:,:,tt) = Plot_volume_rate_vein.cdata;

end

plot2txt(fullTime,total_avg_blood_volume_rate_vein, 'TotalBloodVolumeRateVeins', ToolBox)
plot2txt(fullTime,total_std_blood_volume_rate_vein, 'TotalBloodVolumeRateVeins', ToolBox)

w = VideoWriter(fullfile(ToolBox.PW_path_avi,strcat(ToolBox.main_foldername,'_plot_volume_rate_video_vein.avi')));
tmp = mat2gray( Plot_volume_rate_video_vein);
open(w)
for j = 1:size( Plot_volume_rate_video_vein,4)
    writeVideo(w,tmp(:,:,:,j)) ;  
end
close(w);

%% txt file output with measured pulse wave parameters

fileID = fopen(fullfile(ToolBox.PW_path_txt,strcat(ToolBox.main_foldername,'_pulseWaveOutputParameters.txt')),'a') ;
fprintf(fileID,[...
    'Value of total arterial blood volume rate (µL/min) :\n%d\n' ...
    'Value of total venous blood volume rate (µL/min) :\n%d\n'], ...
    mean_volume_rate_artery(1), ...
    mean_volume_rate_vein(1));
fclose(fileID) ;

for ii=1:size(avg_blood_volume_rate_artery,1)
    fileID = fopen(fullfile(ToolBox.PW_path_txt,strcat(ToolBox.main_foldername,'_pulseWaveOutputParameters.txt')),'a') ;
    fprintf(fileID,[...
        'Artery n°A%d : cross_section (mm^2) : \n %d \n' ...
        'Artery n°A%d : vessel diameter (µm) : \n %d \n' ...
        'Artery n°A%d : average velocity (mm/s) : \n %d \n' ...
        'Artery n°A%d : blood rate (µL/min) : \n %d \n \n'], ...
        ii, ...
        cross_section_area_artery(ii), ...
        ii, ...
        2*sqrt(cross_section_area_artery(ii)/pi)*1000, ... % calculation of the diameter knowing the disc area
        ii, ...
        avg_blood_velocity_artery(ii), ...
        ii, ...
        avg_blood_volume_rate_artery(ii));% mm^3/s -> µL/min 
    fclose(fileID) ;
end

for ii=1:size(avg_blood_volume_rate_vein,1)
    fileID = fopen(fullfile(ToolBox.PW_path_txt,strcat(ToolBox.main_foldername,'_pulseWaveOutputParameters.txt')),'a') ;
    fprintf(fileID,[...
        'Vein n°V%d : cross_section (mm^2) : \n %d \n ' ...
        'Vein n°V%d : vessel diameter (µm) : \n %d \n ' ...
        'Vein n°V%d : average velocity (mm/s) : \n %d \n ' ...
        'Vein n°V%d : blood rate (µL/min) : \n %d \n \n'], ...
        ii, ...
        cross_section_area_vein(ii), ...
        ii, ...
        2*sqrt(cross_section_area_vein(ii)/pi)*1000, ... % calculation of the diameter knowing the disc area
        ii, ...
        avg_blood_velocity_vein(ii), ...
        ii, ...
        avg_blood_volume_rate_vein(ii));% mm^3/s -> µL/min 
    fclose(fileID) ;
end

%% Saving figures

%print('-f111222','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_Artery_Sections.png'))) ;
%print('-f111223','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_Vein_Sections.png'))) ;
%print('-f121','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_mean_blood_volume_rate_in_arteries.png'))) ;
print('-f125','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_Plot_blood_volume_rate_in_arteries.png'))) ;
%print('-f122','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_mean_blood_volume_rate_in_veins.png'))) ;
print('-f126','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_Plot_blood_volume_rate_in_veins.png'))) ;
print('-f57','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_Plot_velocity_in_artery_section.png'))) ;
print('-f58','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_Plot_velocity_in_vein_section.png'))) ;



%print('-f115','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_AVG_Velocity_in_arteriy_sections.png'))) ;
%print('-f120','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_Total_blood_volume_rate_in_arteries.png'))) ;

%imwrite(F_total_blood_flow.cdata, fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_AVG_Velocity_in_arteriy_sections.png')));
imwrite(F_Total_blood_volume_rate_artery.cdata, fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_Total_blood_volume_rate_in_arteries.png')));
imwrite(F_Total_blood_volume_rate_vein.cdata, fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_Total_blood_volume_rate_in_veins.png')));
imwrite(F_numerotation.cdata, fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_arteries_veins_numerotation.png')));



%close all 
end