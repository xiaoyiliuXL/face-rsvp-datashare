function HR_mat_perm = label_permutation_func(dat, nIter)

    det_face_task = dat;
    nSubj = size(det_face_task, 3);
    
    % create an empty output mat
    HR_mat_perm = nan(2, 20, nSubj, nIter);
    % 1 row: HR
    % 2 row: refresh rate
    
    %% subject loop
    for iSubj = 1:nSubj
        
        curr_HR_face = det_face_task(1, :,iSubj);
        nIsi = numel(curr_HR_face);
        
        %% iteration loop
        for iIter = 1:nIter
            
            HR_mat_perm(1, :, iSubj, iIter) = randsample(curr_HR_face, nIsi);
            HR_mat_perm(2, :, iSubj, iIter) = det_face_task(2, :,iSubj);
        end
    end
end