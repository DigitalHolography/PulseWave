function artery_mask = createArteryMask(video)
%FIXME add threshold parameter
    mask = std(video, 0, 3);
    mask = imbinarize(im2gray(mask), 'adaptive', 'ForegroundPolarity', 'bright', 'Sensitivity', 0.2);
    pulse = squeeze(mean(video .* mask, [1 2]));

    pulse_init = pulse - mean(pulse, "all");
    C = video;
    for kk = 1:size(video, 3)
        C(:,:,kk) = video(:,:,kk) - squeeze(mean(video, 3));
    end
    pulse_init_3d = zeros(size(video));
    for mm = 1:size(video, 1)
        for pp = 1:size(video, 2)
            pulse_init_3d(mm,pp,:) = pulse_init;
        end
    end
    figure(7)
    y = pulse_init;
    y = y/max(pulse_init);
    plot(y);
    hold on
    C = C .* pulse_init_3d;
    C = squeeze(mean(C, 3));
    artery_mask = C > max(C(:)) * 0.2;

end
