function HR_mat_perm = label_permutation_func(dat, nIter)

    
    nSubj = size(dat, 3);
    
    % create an empty output mat
    HR_mat_perm = nan(2, 20, nSubj, nIter);
    % 1 row: HR
    % 2 row: refresh rate
    
    %% subject loop
    for iSubj = 1:nSubj
        
        curr_HR_face = dat(1, :,iSubj);
        nIsi = numel(curr_HR_face);
        
        %% iteration loop
        for iIter = 1:nIter
            
            HR_mat_perm(1, :, iSubj, iIter) = randsample(curr_HR_face, nIsi);
            HR_mat_perm(2, :, iSubj, iIter) = dat(2, :,iSubj);
        end
    end
end