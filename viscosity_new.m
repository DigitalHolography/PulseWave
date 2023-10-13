function [] = viscosity_new(SubImage_cell , SubVideo_cell,  ToolBox)


nb_section = size(SubImage_cell, 2);
nb_frames = size(SubVideo_cell{1}, 3);
n_interp = 100;
%interpolation parameter
k = 2;

velocity_profiles = zeros(n_interp, nb_frames, nb_section);
velocity_profiles_std = zeros(n_interp, nb_frames, nb_section);

for ii = 1 : nb_section
    subImg = SubImage_cell{ii};
    subVideo =SubVideo_cell{ii};

     %% interpolate
    interp_size = 4*size(subImg, 1)-3;

    subVideo_interp = zeros(interp_size, interp_size, size(subVideo, 3));
    subImg_interp = interp2(subImg, k);
    for frame = 1 : size(subVideo, 3)
        subVideo_interp(:,:, frame) = interp2(subVideo(:,:, frame), k);
    end
    avg_profile = squeeze(sum(subImg_interp, 1)/(size(subImg_interp,1)));
    projVideo = squeeze(sum(subVideo_interp, 1)/size(subVideo_interp,1));
    projVideo_std = squeeze(std(subVideo_interp,0,1));
    list = find(avg_profile>(0.1*max(avg_profile,[],'all')));


    x = 1:size(projVideo,1);
    xinterp_wall2wall = linspace(list(1),list(end),n_interp);

    

    velocityProfileInterp = zeros(n_interp, size(projVideo,2));
    for tt = 1:size(projVideo,2)
        velocityProfileInterp(:,tt) = interp1(x, projVideo(:,tt), xinterp_wall2wall);
    end
    velocity_profiles(:, :, ii) = velocityProfileInterp;

    velocityProfileInterp_std = zeros(n_interp, size(projVideo,2));
    for tt = 1:size(projVideo,2)
        velocityProfileInterp_std(:,tt) = interp1(x, projVideo_std(:,tt), xinterp_wall2wall);
    end
    velocity_profiles_std(:, :, ii) = velocityProfileInterp_std;

end% ii (artery #)

average_velocity_profile = squeeze(mean(velocity_profiles, 3));
average_velocity_profile_std = squeeze(mean(velocity_profiles_std, 3));

v = VideoWriter(fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername,'_velocity_profile.avi')));% avi
vMP4 = VideoWriter(fullfile(ToolBox.PW_path_mp4, strcat(ToolBox.main_foldername,'_velocity_profile.mp4')),'MPEG-4');% mp4
open(v);
open(vMP4);
mimin = min(average_velocity_profile(:));
[mamax, idx_mamax] = max(average_velocity_profile(:));
[~,idx_syst] = ind2sub(size(average_velocity_profile),idx_mamax);
% x for normalize wall length for fiting
x = linspace(-1,1,length(average_velocity_profile(:,1)));

Vmax_list = zeros(nb_frames,1);
alpha_list = zeros(nb_frames,1);
beta_list = zeros(nb_frames,1);
eta_list = zeros(nb_frames,1);
viscosity_list = zeros(nb_frames,1);

average_velocity_profile_systole = average_velocity_profile(:,idx_syst);
average_velocity_profile_diastole = average_velocity_profile(:,end);

% 
% curve1 = total_avg_blood_volume_rate_vein+0.5*total_std_blood_volume_rate_vein;
% curve2 = total_avg_blood_volume_rate_vein-0.5*total_std_blood_volume_rate_vein;
% fullTime2 = [fullTime, fliplr(fullTime)];
% inBetween = [curve1', fliplr(curve2')];
% fill(fullTime2, inBetween, Color_std);
% hold on;
% plot(fullTime,curve1,"Color",Color_std, 'LineWidth', 2);
% plot(fullTime, curve2, "Color",Color_std, 'LineWidth', 2);
% plot(fullTime,total_avg_blood_volume_rate_vein,'-k','LineWidth',1);
% axis tight ;
% hold off

Color_std = [0.8 0.8 0.8];
fullTime = 1:n_interp;


for tt = 1 : nb_frames
    tmp_velocity_profile = squeeze(average_velocity_profile(:,tt));
    tmp_velocity_profile_plus_std = tmp_velocity_profile + 0.5*average_velocity_profile_std(:,tt);
    tmp_velocity_profile_minus_std = tmp_velocity_profile - 0.5*average_velocity_profile_std(:,tt);
    inBetween = [tmp_velocity_profile_plus_std', fliplr(tmp_velocity_profile_minus_std')];
    fullTime2 = [fullTime, fliplr(fullTime)];
    


    % Use the defined function as an input to fit the function of viscosity
    tmp_fittype = fittype('Vmax .* (1-(1-alpha).* (abs(0.7*x).^beta))',...
    'dependent',{'tmp_velocity_profile'},'independent',{'x'},...
    'coefficients',{'Vmax','alpha','beta'});
    % tmp_fittype = fittype('Vmax .* (1-(1-0.13).* (abs(x).^beta))',...
    % 'dependent',{'tmp_velocity_profile'},'independent',{'x'},...
    % 'coefficients',{'Vmax','beta'});
    [tmp_fit, R2_tmp_fit] = fit(x', tmp_velocity_profile, tmp_fittype, 'StartPoint', [40 0.7 2],'Lower', [0 -5 1.5], 'Upper', [80 3 6]);
    R2_tmp_fit = R2_tmp_fit.rsquare;
    
    fifig = figure(899);
    
    fill(fullTime2,inBetween, Color_std);
    hold on
    plot(tmp_velocity_profile,'-k', 'LineWidth',2) ;
    plot(tmp_velocity_profile_plus_std,"Color",Color_std, 'LineWidth',2)
    plot(tmp_velocity_profile_minus_std,"Color",Color_std, 'LineWidth',2)
    plot(tmp_fit(x), '-r', 'LineWidth',2);
    title('average wall-to-wall arterial velocity profile');
    legend(strcat('R² = ',string(round(R2_tmp_fit,2)),' Vmax = ', string(round(tmp_fit.Vmax,1)),' alpha = ', string(round(tmp_fit.alpha,1)),' beta = ', string(round(tmp_fit.beta,1))));
    fontsize(gca,12,"points") ;
    xticks(x);
    xticklabels({'-1','wall start', '-0.6','-0.4', '-0.2','0','0.2','0.4', '0.6', '0.8','wall end'});
    xlabel('section','FontSize',14) ;
    pbaspect([1.618 1 1]) ;
    set(gca, 'LineWidth', 2);
    axis tight;
    ylim([mimin mamax]);
    ylabel('quantitative velocity mm/s','FontSize',14) ;
    hold off
    writeVideo(v, getframe(fifig));
    writeVideo(vMP4, getframe(fifig));


    Vmax_list(tt) = tmp_fit.Vmax;
    alpha_list(tt) = tmp_fit.alpha;
    % alpha_list(tt) = 0.13;
    beta_list(tt) = tmp_fit.beta;
    eta_list(tt) = (tmp_fit.beta + 1)/(tmp_fit.beta + tmp_fit.alpha);
    viscosity_list(tt) = -(eta_list(tt)-1.459)/0.017;

end
close(v)
close(vMP4)
video = subVideo;

% Systole/Diastole velocity profile

x_section = linspace(-0.7,0.7,length(squeeze(average_velocity_profile_systole)));
fit_velocity_profile_systole = Vmax_list(idx_syst)*(1-(1-alpha_list(idx_syst)).*abs(x_section).^beta_list(idx_syst));
fit_velocity_profile_diastole = Vmax_list(end)*(1-(1-alpha_list(end)).*abs(x_section).^beta_list(end));


figure(668)
plot(x_section,average_velocity_profile_systole,'-k',...
x_section,average_velocity_profile_diastole,'-k',...
    x_section,fit_velocity_profile_systole,'-r',...
    x_section,fit_velocity_profile_diastole,'-r', 'LineWidth',2)
title('Systole and diastole arterial velocity profile');
fontsize(gca,12,"points") ;
% xticks(x);
xticklabels({'section'});
xlabel('section','FontSize',14) ;
pbaspect([1.618 1 1]) ;
set(gca, 'LineWidth', 2);
axis tight;
ylim([0.9*mimin 1.1*mamax]);
ylabel('velocity (mm/s)','FontSize',14) ;


figure(666)
plot(viscosity_list)
pbaspect([1.618 1 1]);
xlabel('Frame','FontSize',14);
ylabel('Viscosity (cP)','FontSize',14); 
set(gca, 'LineWidth', 2);
axis tight;

% png
print('-f668','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_velocity_cross_section.png'))) ;
print('-f666','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_viscosity_in_time.png'))) ;
% eps
print('-f668','-depsc',fullfile(ToolBox.PW_path_eps,strcat(ToolBox.main_foldername,'_velocity_cross_section.eps'))) ;
print('-f666','-depsc',fullfile(ToolBox.PW_path_eps,strcat(ToolBox.main_foldername,'_viscosity_in_time.eps'))) ;

list_fig_close = [666, 899];
for ii=1:length(list_fig_close)
    close(list_fig_close(ii));
end

end