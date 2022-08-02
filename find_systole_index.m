 function [sys_index_list, mask] = find_systole_index(video)
    arguments
        video 
    end
    %Video retrieval frame by frame
    

    for pp = 1:size(video, 3)
        video(:,:,pp) = video(:,:,pp) ./ mean(video(:,:,pp), [1 2]);
    end

%     mask = std(video, 0, 3);
%     mask = imbinarize(im2gray(mask), 'adaptive', 'ForegroundPolarity', 'bright', 'Sensitivity', 0.2);
%     pulse = squeeze(mean(video .* mask, [1 2]));
% 
%     pulse_init = pulse - mean(pulse, "all");
%     C = video;
%     for kk = 1:size(video, 3)
%         C(:,:,kk) = video(:,:,kk) - squeeze(mean(video, 3));
%     end
%     pulse_init_3d = zeros(size(video));
%     for mm = 1:size(video, 1)
%         for pp = 1:size(video, 2)
%             pulse_init_3d(mm,pp,:) = pulse_init;
%         end
%     end
%     C = C .* pulse_init_3d;
%     C = squeeze(mean(C, 3));
%     mask = C > max(C(:))*0.2;

    mask = createArteryMask(video) ;

%     figure(1)
%     imagesc(mask);
%     title('Artery mask (no time lag)');
%     axis off
%     axis equal
%     colormap gray



    pulse = squeeze(mean(video .* mask, [1 2]));

    pulse_init = pulse - mean(pulse, "all");
    y = pulse_init;
    y = y/max(pulse_init);
    plot(y, 'r--');
    hold off
    pulse_init = detrend(pulse_init);


%     
%     pulse_init = flip(pulse_init,1); %FIXE ME 


    figure(33)
    plot(pulse_init,'-k','LineWidth',2) ; 
    title('pulse init');
    fontsize(gca,12,"points") ;
    xlabel('Frames','FontSize',14) ;
    pbaspect([1.618 1 1]) ;
    set(gca, 'LineWidth', 2);
    axis tight

    diff_signal = diff(pulse_init);
%     figure(3)
%     plot(diff_signal);
%     title('pulse init derivate');
%     1;

%     spectrum_signal = fft(diff_signal);
%     w = 1:length(spectrum_signal);
%     figure(41);
%     plot(w, abs(spectrum_signal));
%     signal_center = floor(length(spectrum_signal)/2);
%     [~, omega1] = max(abs(spectrum_signal(1:signal_center)));
%     [~, omega2] = max(abs(spectrum_signal(signal_center:end)));
%     omega2 = omega2 + signal_center - 1;
%     tmp = zeros(1, length(spectrum_signal));
%     tmp(omega1) = spectrum_signal(omega1);
%     tmp(omega2) = spectrum_signal(omega2);
%     tmp = ifft(tmp);
%     figure(42)
%     plot(w, real(tmp), w, diff_signal);
%     pseudo_period = length(w) / omega1;
%     distance_between_peaks = 0.7 * pseudo_period;
% 
%     [~,index_list] = findpeaks(diff_signal,1:length(diff_signal),'MinPeakDistance',distance_between_peaks,'SortStr','descend');
%     [~,index_list_tmp] = findpeaks(tmp,1:length(tmp),'MinPeakDistance',distance_between_peaks,'SortStr','descend');
%     index_list = sort(index_list(1:length(index_list)));
%     fprintf("num Cycles: %d\n", length(index_list)-1);
%     sys_index_list = zeros(1, length(index_list_tmp));
%     for ii = 1:length(sys_index_list)
%         tab = abs(index_list-index_list_tmp(ii));
%         [~, new_index] = min(tab);
%         sys_index_list(ii) = index_list(new_index);
%     end
%     sys_index_list = sort(sys_index_list, 'ascend');
    [~, sys_index_list] = findpeaks(diff_signal, 1:length(diff_signal), 'MinPeakHeight', max(diff_signal) * 0.7);
    1;
end