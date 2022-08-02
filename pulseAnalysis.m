function pulseAnalysis(one_pulse_video, one_cycle_dir, filename, sys_index_list)

[RI_map, ARI] = construct_resistivity_index(one_pulse_video,one_cycle_dir);

%     videoObject = VideoReader(one_pulse_video);
%     h=videoObject.height;
%     w=videoObject.width;
%     
%     NumberFrames=videoObject.NumFrames;
%     
%     I = zeros(h,w, NumberFrames);
% 
%     %Video retrieval frame by frame
%     for i = 1 : NumberFrames
%         Frame = read(videoObject, i);
%         I(:,:,i)=rgb2gray(Frame);
%     end
    
    nb_frames = size(one_pulse_video,3) ;
    blur_time_sys = ceil(nb_frames/50);
    blur_time_dia = ceil(nb_frames/40);

    I = one_pulse_video ; 
    % normalize through laser momentary intensity FIXME
    for pp = 1:size(I, 3)
        I(:,:,pp) = I(:,:,pp) ./ mean(I(:,:,pp), [1 2]);
    end
    
    % identify dominant arteries pulse
    stdPulse = std(I, 0, 3);
    stdPulse = imbinarize(im2gray(stdPulse), 'adaptive', 'ForegroundPolarity', 'bright', 'Sensitivity', 0.05);
%     imwrite(mat2gray(stdPulse), 'stdPulse.png') ; %FIXME : mauvais dossier
    I_arteries = I .* stdPulse;
    % calculate the vector of pulse in arteries
    pulse_arteries = squeeze(mean(I_arteries, [1 2]));
    
    filename_mat = filename(1:length(filename)-3) ;
    cache_exists = exist(fullfile(strcat(strcat(one_cycle_dir,'/../../'),filename_mat,'.mat'))) ;
    if cache_exists
        load(fullfile(strcat(strcat(one_cycle_dir,'/../../'),filename_mat,'.mat'))) ;
        batch_stride = cache.batch_stride ; 
        Fs = cache.Fs ; 
        fileID = fopen(fullfile(one_cycle_dir,'pulse.txt'),'w') ;
        average_cycle_length = 0 ;
        nb_of_averaged_cycles =0;
        if size(sys_index_list,2)==1
            average_cycle_length = nb_frames ;
            nb_of_averaged_cycles =1 
        else 
            for i = 2:size(sys_index_list,2) 
                average_cycle_length = average_cycle_length + (sys_index_list(i)-sys_index_list(i-1)) ; 
                nb_of_averaged_cycles = nb_of_averaged_cycles +1;
            end
            average_cycle_length = average_cycle_length / (length(sys_index_list)-1) ;
        end
        T = linspace(0,batch_stride/Fs*average_cycle_length,nb_frames)  ;
        fprintf(fileID,'%d %d\n',T,pulse_arteries);
        fclose(fileID) ;
    end
%     pulse_init = pulse - mean(pulse, "all");
%     C = I - mean(I, 3);
    C = I;
    pulse_arteries_3d = zeros(size(I_arteries));
    for mm = 1:size(I, 1)
        for pp = 1:size(I, 2)
            pulse_arteries_3d(mm,pp,:) = pulse_arteries;
        end
    end

    
    I0 = I - mean(I, 3);
    for kk = 1:size(I, 3)
        R = I0 .* circshift(pulse_arteries_3d, kk, 3);
        C(:,:,kk) = squeeze(mean(R, 3));
    end
%     C = imgaussfilt3(C, 2);
    [max_C_3, id_max] = max(C, [], 3);
    figure(999);
    imagesc(id_max);
    title('Artery-vein map');
    colormap default
    colorbar;
    axis square
    axis off
    
    figure(33)
    imagesc(max_C_3)
    title('image des coef max de tous les pixels');
    colorbar ;
    fontsize(gca,12,"points") ;
    set(gca, 'LineWidth', 2);
    axis off
    axis image


%     imwrite(mat2gray(single(id_max)),'segmentationAV.png') %FIXME : mauvais dossier
    print('-f999','-dpng',fullfile(one_cycle_dir,strcat(filename,'segmentationAV_fig.png'))) ; 
    

    % pulses 
    nb_frames = size(one_pulse_video,3) ;

    %% arteries mask
    arteries = zeros(size(id_max)); 
    arteries (or(id_max < 0.1*nb_frames, id_max > 0.9 * nb_frames)) = 1 ; 
%     arteries = createArteryMask(one_pulse_video) ; 
    
    
    figure(10)
    imagesc(arteries) ;
    colormap gray
    title('Arteries mask'); 
    fontsize(gca,12,"points") ;
    set(gca, 'LineWidth', 2);
    axis off
    axis image
%     imwrite(mat2gray(single(arteries)),strcat(filename,'_arteries_mask.png'),'png') ; %FIXME : mauvais dossier
    print('-f10','-dpng',fullfile(one_cycle_dir,strcat(filename,'_arteries_mask.png'))) ;
    

    %% veins mask

    veins = zeros(size(id_max)) ; 
    veins(and(id_max>0.1*nb_frames,id_max<0.3*nb_frames)) = 1 ;  

    figure(110)
    imagesc(veins) ;
    colormap gray
    title('Veins mask'); 
    fontsize(gca,12,"points") ;
    set(gca, 'LineWidth', 2);
    axis off
    axis image
    print('-f110','-dpng',fullfile(one_cycle_dir,strcat(filename,'_veins_mask.png'))) ;
   


    %% calculation pulse in arteries
    I_arteries = I .*arteries;
    pulse_arteries = squeeze(sum(I_arteries, [1 2]))/nnz(I_arteries);
    pulse_arteries = pulse_arteries - min(pulse_arteries) ;
    [~,idx] = min(pulse_arteries,[],1) ;
    pulse_arteries = circshift(pulse_arteries,-idx) ;
    
    %% calculation veins in arteries
    I_veins = I .*veins ; 
    pulse_veins = squeeze(sum(I_veins, [1 2]))/nnz(I_veins);
    pulse_veins = pulse_veins - min(pulse_veins) ; 
    pulse_veins = circshift(pulse_veins,-idx) ;

    max_plot = max(max(pulse_arteries(:)),max(pulse_veins(:))) ; 
    pulse_arteries = pulse_arteries ./ max_plot ; 
    pulse_veins = pulse_veins ./ max_plot ;


    %% plot pulses in veins and arteries


   sys_index_list_one_cycle = find_systole_index(I_arteries) - idx ;


    if cache_exists
        [~,idx_sys] = max(pulse_arteries) ;
 
        T_syst = batch_stride/Fs*average_cycle_length * (idx_sys-1) / nb_frames ; 

        
      % Si .raw, utiliser les 2 lignes de flip ci-dessous :
%               T_syst = T(nb_frames) - T_syst;
%         pulse_arteries = flip(pulse_arteries,1);
%         pulse_veins = flip(pulse_veins,1);
       
        figure(70)
        plot(T,pulse_arteries,'k-',T,pulse_veins,'k--',LineWidth=2) ;
        xline(T_syst,':',{'Systolic peak      '},LineWidth=2) ;
        legend('Arteries','Veins') ;
        fontsize(gca,12,"points") ;
        xlabel('Time (s)','FontSize',14) ;
        ylabel('Normalized intensity','FontSize',14) ;
        pbaspect([1.618 1 1]) ;
        set(gca, 'LineWidth', 2);
        title('normalized pulse wave in arteries and veins');
       
        axis tight
    

        pulse_arteries_blurred_sys = movavgvar(pulse_arteries(1:idx_sys), blur_time_sys);
        diff_pulse_sys = diff(pulse_arteries_blurred_sys) ;
        pulse_arteries_blurred_dia = movavgvar(pulse_arteries(idx_sys:end), blur_time_dia);
        diff_pulse_dia = diff(pulse_arteries_blurred_dia) ;
        idx_list_threshold_dia = find(pulse_arteries_blurred_dia(:) < max(pulse_arteries_blurred_dia(:))/exp(1));
        [max_diff_pulse,idx_T_diff_max] = max(diff_pulse_sys); %faire ca mieux
        acc = max_diff_pulse/(T(idx_T_diff_max)-T(idx_T_diff_max-1));
       

        figure(80)
        plot(T,pulse_arteries,'k.', ...
            T(1:idx_sys),pulse_arteries_blurred_sys(1:idx_sys),'k-', ...
            T(idx_sys:nb_frames),pulse_arteries_blurred_dia(1:(nb_frames-idx_sys+1)),'k-', ...
            LineWidth=2) ;
        xline(T(idx_T_diff_max-1),':',{},LineWidth=2)
        text(T(idx_T_diff_max-1),0.1,' (1)','FontSize',14);
        xline(T_syst,':',{},LineWidth=2) ;
        text(T_syst,0.1,'  (2)','FontSize',14);
        xline(T(idx_sys+idx_list_threshold_dia(1)),':',{},LineWidth=2);
        text(T(idx_sys+idx_list_threshold_dia(1)),0.1,'  (3)','FontSize',14);
%                 yline(1/exp(1),':',LineWidth=2); 
        legend('data','smoothed line') ;
        fontsize(gca,12,"points") ;
        xlabel('Time (s)','FontSize',14) ;
        ylabel('Normalized intensity','FontSize',14) ;
        pbaspect([1.618 1 1]) ;
        set(gca, 'LineWidth', 2);
        title('Normalized average pulse wave');
        axis tight




%         coord_T = [T(idx_T_diff_max-1) T(idx_sys+idx_list_threshold_dia(1))];
%         coord_y = [pulse_arteries_blurred_sys(idx_T_diff_max) pulse_arteries_blurred_dia(idx_list_threshold_dia(1))];
%         str = {' \leftarrow max. increase rate', '\leftarrow 1/e threshold' };
%         text(coord_T, coord_y,str,'FontSize',14);

        figure(90)
        plot(T(1:idx_sys-1), diff_pulse_sys(1:idx_sys-1),'k-', T(idx_sys-1:nb_frames-2), diff_pulse_dia(1:(nb_frames-idx_sys)), 'k-', LineWidth=2);
        x=0;
        yline(x,':',LineWidth=2) ;
        fontsize(gca,12,"points") ;
        xlabel('Time (s)','FontSize',14) ;
        ylabel('Pulse derivative','FontSize',14) ;
        pbaspect([1.618 1 1]) ;
        set(gca, 'LineWidth', 2);
        title('Derivative of pulse wave');
        axis tight


    else % no .mat present, hence no timeline

        [~,idx_sys] = max(pulse_arteries) ; 
        T_syst = idx_sys ; 
        T=1:nb_frames;

              % Si .raw, utiliser les 2 lignes de flip ci-dessous :
%         pulse_arteries = flip(pulse_arteries,1);
%         pulse_veins = flip(pulse_veins,1);
       
       figure(70)
        plot(T,pulse_arteries,'k-',T,pulse_veins,'k--',LineWidth=2) ;
        xline(T_syst,':',{'Systolic peak      '},LineWidth=2) ;
        legend('Arteries','Veins') ;
        fontsize(gca,12,"points") ;
        xlabel('Time (s)','FontSize',14) ;
        ylabel('Normalized intensity','FontSize',14) ;
        pbaspect([1.618 1 1]) ;
        set(gca, 'LineWidth', 2);
        title('Normalized pulse wave in arteries and veins');
        axis tight


        pulse_arteries_blurred_sys = movavgvar(pulse_arteries(1:idx_sys), blur_time_sys);
        diff_pulse_sys = diff(pulse_arteries_blurred_sys) ;
        pulse_arteries_blurred_dia = movavgvar(pulse_arteries(idx_sys:end), blur_time_dia);
        diff_pulse_dia = diff(pulse_arteries_blurred_dia) ;
        idx_list_threshold_dia = find(pulse_arteries_blurred_dia(:) < max(pulse_arteries_blurred_dia(:))/exp(1));
        [max_diff_pulse,idx_T_diff_max] = max(diff_pulse_sys); %faire ca mieux
        acc = max_diff_pulse/(T(idx_T_diff_max)-T(idx_T_diff_max-1));


        figure(80)
        plot(T,pulse_arteries,'k.', ...
            T(1:idx_sys),pulse_arteries_blurred_sys(1:idx_sys),'k-', ...
            T(idx_sys:nb_frames),pulse_arteries_blurred_dia(1:(nb_frames-idx_sys+1)),'k-', ...
            LineWidth=2) ;
        xline(T(idx_T_diff_max-1),':',{},LineWidth=2)
        text(T(idx_T_diff_max-1),0.1,' (1)','FontSize',14);
        xline(T_syst,':',{},LineWidth=2) ;
        text(T_syst,0.1,'  (2)','FontSize',14);
        xline(T(idx_sys+idx_list_threshold_dia(1)),':',{},LineWidth=2);
        text(T(idx_sys+idx_list_threshold_dia(1)),0.1,'  (3)','FontSize',14);
%                 yline(1/exp(1),':',LineWidth=2); 
        legend('data','smoothed line') ;
        fontsize(gca,12,"points") ;
        xlabel('frames','FontSize',14) ;
        ylabel('Normalized intensity','FontSize',14) ;
        pbaspect([1.618 1 1]) ;
        set(gca, 'LineWidth', 2);
        title('Normalized average pulse wave');
        axis tight

       
        [max_diff_pulse,idx_T_diff_max] = max(diff_pulse_sys); %faire ca mieux
        idx_list_threshold_dia = find(pulse_arteries_blurred_dia(:) < 1/exp(1));
        coord_T = [T(idx_T_diff_max-1) T(idx_sys+idx_list_threshold_dia(1))];
        coord_y = [pulse_arteries_blurred_sys(idx_T_diff_max) pulse_arteries_blurred_dia(idx_list_threshold_dia(1))];
        str = {' \leftarrow max. increase rate', '\leftarrow 1/e threshold' };
        text(coord_T, coord_y,str,'FontSize',14);

       figure(90)
        plot(T(1:idx_sys-1), diff_pulse_sys(1:idx_sys-1),'k-', T(idx_sys-1:nb_frames-2), diff_pulse_dia(1:(nb_frames-idx_sys)), 'k-', LineWidth=2);
        x=0;
        yline(x,':',LineWidth=2) ;
        fontsize(gca,12,"points") ;
        xlabel('Time (s)','FontSize',14) ;
        ylabel('Pulse derivative','FontSize',14) ;
        pbaspect([1.618 1 1]) ;
        set(gca, 'LineWidth', 2);
        title('Approximate derivate of pulse wave');
        axis tight

    end
    T_syst = T(idx_sys);
    disp('T syst');
    disp(T_syst);
%     systole_slope = 1 / T_syst ; 
    systole_area = sum(pulse_arteries(1:idx_sys)) ;
    diastole_area = sum(pulse_arteries(idx_sys:end));
%     systole_area = sum(pulse_arteries(1:T_syst)) ;
%     diastole_area = sum(pulse_arteries(T_syst:end));
    cc=systole_area;
    systole_area = systole_area/(diastole_area+systole_area);
    diastole_area = diastole_area/(diastole_area+cc);
 
    nb_of_detected_systoles = size(sys_index_list,2) ;
   
    [max_diff_pulse,idx_T_diff_max] = max(diff_pulse_sys); %faire ca mieux
    text(pulse_arteries_blurred_sys(idx_T_diff_max),idx_T_diff_max,'\leftarrow max. rate')

    [min_diff_pulse,idx_T_diff_min] = min(diff_pulse_dia);%faire ca mieux

% txt file output with measured pulse wave parameters  
    fileID = fopen(fullfile(one_cycle_dir,'pulseWaveParameters.txt'),'w') ;
        fprintf(fileID,[...
         'Value of pulse derivative at the maximum systolic increase :\n%d\n' ...
         'Maximal acceleration (m/s^(-2)) :\n%d\n' ...
         'Time of maximum systolic increase (s) :\n%d\n' ...
         'Time of systolic peak (s) :\n%d\n' ...
         'Time of the intersection between the diastolic descent and the threshold 1/e (s) :\n%d\n' ...
         'Number of detected systoles :\n%d\nNumber of averaged cycles :\n%d\n' ...
         'Area under the systolic rise curve :\n%d\n' ...
         'Area under the diastolic descent curve  :\n%d\n' ...
         'Average arterial resistivity index :\n%d\n'], ...
        max_diff_pulse, ...
        acc, ...
        T(idx_T_diff_max), ...
        T_syst, ...
        T(idx_sys+idx_list_threshold_dia(1)), ...
        nb_of_detected_systoles, ...
        nb_of_averaged_cycles, ...
        systole_area, ...
        diastole_area, ...
        ARI); 
    fclose(fileID) ;

    print('-f70','-dpng',fullfile(one_cycle_dir,strcat(filename,'_pulse_AV.png'))) ;
    print('-f80','-dpng',fullfile(one_cycle_dir,strcat(filename,'_pulse_wave.png'))) ;
    print('-f90','-dpng',fullfile(one_cycle_dir,strcat(filename,'_diff_pulse.png'))) ;
    print('-f70','-depsc',fullfile(one_cycle_dir,strcat(filename,'_pulse_AV.eps'))) ;
    print('-f80','-depsc',fullfile(one_cycle_dir,strcat(filename,'_pulse_wave.eps'))) ;
    print('-f90','-depsc',fullfile(one_cycle_dir,strcat(filename,'_diff_pulse.eps'))) ;

end