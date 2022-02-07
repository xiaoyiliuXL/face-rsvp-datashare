%% load
clc; clear all;

load('../data/raw_data.mat');
% A single matrix with all the data (480 trials each subject)

    % 1 row: target gender: 0-female; 1-male
    % 2 row: participant gender
    % 3 row: target race
    % 4 row: participant race
    % 5 row: isi
    % 6 row: face detection acc
    % 7 row: gender detection acc
    % 8 row: participant age
    % 9 row: refresh rate
%% Get average HR for each condition

nFem = sum(squeeze(mat_data(2,1,:)) == 0);
nMale = sum(~squeeze(mat_data(2,1,:)) == 0);
nSubj = nFem + nMale;

big_HR_mat = nan(4, 20, nSubj);
% big average HR matrix (20 ISIs)
    % 1 row: age
    % 2 row: face task acc
    % 3 row: gender acc
    % 4 row: refresh rate
    
%% participant loop

for iSubj = 1:nSubj
    
    % get current subject's data
    curr_mat = mat_data(:,:,iSubj);
    isi_vals = unique(curr_mat(5,:));
    curr_age = unique(curr_mat(8,:));
    
    curr_refresh_rate = unique(curr_mat(9,:));
    
    %% get average HR
    loopIsi = 0;
    
    % ISI loops
    for iIsi = isi_vals
        
        loopIsi = loopIsi + 1;
        
        % current isi mask
        isi_mask = curr_mat(5,:) == iIsi;
        
        % compute HR at each ISI condition
        curr_HR_face = nanmean(curr_mat(6,isi_mask));
        curr_HR_gender = nanmean(curr_mat(7,isi_mask));
        
        % put the average HR in current ISI condition in place
        big_HR_mat(2, loopIsi, iSubj) = curr_HR_face;
        big_HR_mat(3, loopIsi, iSubj) = curr_HR_gender;
        
    end
    
    % participant loop
    big_HR_mat(1, :, iSubj) = repelem(curr_age, 20);
    big_HR_mat(4, :, iSubj) = repelem(curr_refresh_rate, 20);
end

save('big_HR_mat.mat', 'big_HR_mat')
