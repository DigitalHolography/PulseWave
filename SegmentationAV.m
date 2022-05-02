function SegmentationAV(one_pulse_video)

    videoObject = VideoReader(one_pulse_video);
    h=videoObject.height;
    w=videoObject.width;
    
    NumberFrames=videoObject.NumFrames;
    
    I = zeros(h,w, NumberFrames);

    %Video retrieval frame by frame
    for i = 1 : NumberFrames
        Frame = read(videoObject, i);
        I(:,:,i)=rgb2gray(Frame);
    end
    
    % normalize through laser momentary intensity FIXME
    for pp = 1:size(I, 3)
        I(:,:,pp) = I(:,:,pp) ./ mean(I(:,:,pp), [1 2]);
    end
    
    % identify dominant arteries pulse
    mask_arteries = std(I, 0, 3);
    mask_arteries = imbinarize(im2gray(mask_arteries), 'adaptive', 'ForegroundPolarity', 'bright', 'Sensitivity', 0.2);
    I_arteries = I .* mask_arteries;

    % calculate the vector of pulse in arteries
    pulse_arteries = squeeze(mean(I_arteries, [1 2]));

%     pulse_init = pulse - mean(pulse, "all");
%     C = I - mean(I, 3);
    C = I;
    pulse_arteries_3d = zeros(size(I_arteries));
    for mm = 1:size(I, 1)
        for pp = 1:size(I, 2)
            pulse_arteries_3d(mm,pp,:) = pulse_arteries;
        end
    end

    tic
    I0 = I - mean(I, 3);
    for kk = 1:size(I, 3)
        R = I0 .* circshift(pulse_arteries_3d, kk, 3);
        C(:,:,kk) = squeeze(mean(R, 3));
    end
%     C = imgaussfilt3(C, 2);
    [max_C_3, id_max] = max(C, [], 3);
    figure(45);
    imagesc(id_max);
    toc

end