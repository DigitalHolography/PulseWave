function [] = velocityHistogramm(v_RMS_all, maskArtery,maskVein ,ToolBox)
%% Init of histogram axis
[size_vRMS_1,size_vRMS_2,size_vRMS_3] = size(v_RMS_all);

v_histo_artery = round(v_RMS_all.*maskArtery);
v_min_artery = min(v_histo_artery,[],'all');
v_max_artery = max(v_histo_artery,[],'all');

v_histo_vein = round(v_RMS_all.*maskVein); 
v_min_vein = min(v_histo_vein,[],'all');
v_max_vein = max(v_histo_vein,[],'all');

v_max_all = max(v_max_artery,v_max_vein);
v_min_all = min(v_min_artery,v_min_vein);
v_max_all_display = round(0.8*v_max_all);
v_min_all_display = round(0.8*v_min_all);



yAx = [v_min_all v_max_all ] ;
%yAx_display = [-20  80] ; 
yAx_display = yAx;
%FIXME trouver un moyen de croper proprement sans décaler le zéro


%% Velocity Histogram in arteries
%FIXME prctile 10% Y = percentil(X,[5 95])

X = linspace(v_min_all,v_max_all,v_max_all-v_min_all+1);
n = size(X,2);
xAx = [0 n*ToolBox.stride/(1000*ToolBox.fs)];
histo_artery = zeros(size(X,2),size_vRMS_3);
%histo_video_artery = zeros(size(X,2),size(dataCubeM2M0,3),3,size(dataCubeM2M0,3));

f_distrib_artery = figure(157);
f_distrib_artery.Position(3:4) = [600 275];
index_min = find(X == v_min_all_display);
index_max = find(X == v_max_all_display);
imagesc(xAx,yAx_display,histo_artery(index_min:index_max,:))
set(gca,'YDir','normal')
set(gca,'PlotBoxAspectRatio',[2.5 1 1])
colormap("hot")
f = getframe(gcf);
[M,N,~] = size(f.cdata);
histo_video_artery = zeros(M,N,3,size_vRMS_3);


%FIXME avoir une ligne à zéro de trois pixel 
%FIXME getframe pour la couleur 
for t = 1:size_vRMS_3
    for x = 1:size_vRMS_1
        for y = 1:size_vRMS_2
            if ( v_histo_artery(x,y,t) ~= 0) 
            i = find(X == v_histo_artery(x,y,t) );
            histo_artery(i,t) = histo_artery(i,t) + 1;
            end


        end
    end
    %histo_video_artery(:,:,t) = flip(histo_artery,1);
    figure(157)
    imagesc(xAx,yAx_display,histo_artery(index_min:index_max,:))
    set(gca,'YDir','normal')
    title("Velocity distribution in arteries")
    ylabel('Velocity (mm.s^{-1})')
    xlabel('Time (s)')
    set(gca,'PlotBoxAspectRatio',[2.5 1 1])
    f = getframe(gcf);
    histo_video_artery(:,:,:,t) = imresize(f.cdata,[M N]);

end 
velocity_dist_arteries = frame2im(getframe(gca));

w = VideoWriter(fullfile(ToolBox.PW_path_avi,strcat(ToolBox.main_foldername,'_velocity_histogram_artery.avi')));
tmp = mat2gray(histo_video_artery);
open(w)
for j = 1:size(histo_video_artery,4)
    writeVideo(w,tmp(:,:,:,j)) ;  
end
close(w);


%% Velocity Histogram in veins

X = linspace(v_min_all,v_max_all,v_max_all-v_min_all+1);
n = size(X,2);
xAx = [0 n*ToolBox.stride/(1000*ToolBox.fs)];
histo_vein = zeros(size(X,2),size_vRMS_3);
%histo_video_vein = zeros(size(X,2),size(dataCubeM2M0,3),size(dataCubeM2M0,3));
f = getframe(gcf);
[M,N,~] = size(f.cdata);
histo_video_vein = zeros(M,N,3,size_vRMS_3);



f_distrib_vein = figure(158);
f_distrib_vein.Position(3:4) = [600 275];
imagesc(xAx,yAx_display,histo_vein(index_min:index_max,:))
set(gca,'YDir','reverse')
set(gca,'PlotBoxAspectRatio',[2.5 1 1])
colormap("bone")



for t = 1:size_vRMS_3
    for x = 1:size_vRMS_1
        for y = 1:size_vRMS_2
            if( v_histo_vein(x,y,t) ~= 0)
            i = find(X == v_histo_vein(x,y,t) );
            histo_vein(i,t) = histo_vein(i,t) + 1;
            end


        end
    end
    %histo_video_vein(:,:,t) = flip(histo_vein,1);
    figure(158)
    imagesc(xAx,yAx_display,histo_vein(index_min:index_max,:))
    set(gca,'YDir','normal')
    ylabel('Velocity (mm.s^{-1})')
    xlabel('Time (s)')
    title("Velocity distribution in veins")
    set(gca,'PlotBoxAspectRatio',[2.5 1 1])
    f = getframe(gcf);
    %histo_video_vein(:,:,:,t) = imresize(f.cdata,[M N]);
    histo_video_vein(:,:,:,t) = f.cdata;
end 

velocity_dist_veins = frame2im(getframe(gca));

w = VideoWriter(fullfile(ToolBox.PW_path_avi,strcat(ToolBox.main_foldername,'_velocity_histogram_veins.avi')));
tmp = mat2gray(histo_video_vein);
open(w)
for j = 1:size(histo_video_vein,4)
    writeVideo(w,tmp(:,:,:,j)) ;  
end
close(w);

print('-f157','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_velocity_distribution_arteries_fullcycle.png'))) ;
print('-f158','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_velocity_distribution_veins_fullcycle.png'))) ;

imwrite(velocity_dist_veins,fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_Velocity_distribution_in_veins.png')),'png') ;
imwrite(velocity_dist_arteries,fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_Velocity_distribution_in_arteries.png')),'png') ;
end