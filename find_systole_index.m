function [sys_index_list, fullPulseWave] = find_systole_index(video,path,maskArtery)


PW_params = Parameters(path);

fullPulseWave = squeeze(sum(video .*maskArtery, [1 2])/nnz(maskArtery));
pulse_init = detrend(fullPulseWave);
pulse_init = pulse_init - mean(pulse_init, "all");

[~,fullArterialPulseRmOut] = discardPulseWaveOutliers(pulse_init,3);
fullArterialPulseClean = smoothdata(fullArterialPulseRmOut,'lowess');

diff_signal = diff(fullArterialPulseClean);
[~, sys_index_list] = findpeaks(diff_signal, 1:length(diff_signal), 'MinPeakHeight', max(diff_signal) * PW_params.systoleThreshold);

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

end