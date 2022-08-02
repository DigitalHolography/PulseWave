function one_cycle_video = create_one_cycle(video, mask, sys_index_list, Ninterp)
    %   one_cycle() : identifies pulse cycles and average them to one video
    %   sys_index_list : list of systole indexes in video
    arguments
        video 
        mask
        sys_index_list {mustBeNumeric, mustBeNonnegative} = []
        Ninterp {mustBeNumeric, mustBePositive} = 64
    end
    %Video retrieval frame by frame
    
    M = length(sys_index_list)-1;
    one_cycle_video = zeros(size(video,1), size(video,2), Ninterp, M);
    for ii = 1:M
        for id_x = 1 : size(video,1)
            for id_y = 1 : size(video,2)
                interp_range = linspace(sys_index_list(ii),sys_index_list(ii+1)-1,Ninterp);
                one_cycle_video(id_x,id_y,:,ii) = interp1((sys_index_list(ii):sys_index_list(ii+1)-1),squeeze(video(id_x, id_y,sys_index_list(ii):sys_index_list(ii+1)-1)),interp_range);
            end
        end
    end

%     one_cycle_video = mat2gray(squeeze(mean(one_cycle_video, 4)));
% FIXME
    one_cycle_video = (squeeze(mean(one_cycle_video, 4)));

%     shift = floor(size(one_cycle_video, 3) / 32) 
% FIXME
%     mask = createArteryMask(one_cycle_video) ;
    signal = squeeze(sum(one_cycle_video .* mask,[1 2 ])./ nnz(mask));
   
    [min_val,shift] = min(signal)  ; %problème lorsqu'on modifie la valeur de shift à cette ligne
    disp([min_val,shift]) ; 
     one_cycle_video = circshift(one_cycle_video, -shift, 3);

    one_cycle_video = mat2gray(mean(one_cycle_video, 4));
end