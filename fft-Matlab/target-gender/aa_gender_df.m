%% Create dataframe for gender effects analysis %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;

%% load data files

load('../data/big_HR_mat.mat');
    % 1 row: target gender(0-female; 1-male)
    % 2 row: gender group (0-same; 1-other)
    % 3 row: age
    % 4 row: isi
    % 5 row: face task acc
    % 6 row: gender acc
    % 7 row: refresh rate

nSubj = size(big_HR_mat, 3);

%% female target face data
big_HR_t_fem_mat = big_HR_mat(:, [1:20], :);

HR_t_fem_mat = nan(4, 20, nSubj);
% 1 row: age
% 2 row: face task acc
% 3 row: gender task acc
% 4 row: refresh rate

for iSubj = 1:nSubj
    curr_mat = squeeze(big_HR_t_fem_mat(:, :, iSubj));
    HR_t_fem_mat(1, :, iSubj) = repelem(unique(curr_mat(3,:)), 20);
    HR_t_fem_mat(4, :, iSubj) = repelem(unique(curr_mat(7,:)), 20);
    
    isi_vals = unique(curr_mat(4,:));
    
    % compute mean HR for each isi
    loop_isi = 0;
    for iIsi = isi_vals
        loop_isi = loop_isi + 1;
        
        % isi mask
        isi_mask = curr_mat(4, :) == iIsi;
        
        % mean HR
        curr_HR_face = mean(curr_mat(5, isi_mask));
        curr_HR_gender = mean(curr_mat(6, isi_mask));
        
        % add to data frame
        HR_t_fem_mat(2, loop_isi, iSubj) = curr_HR_face;
        HR_t_fem_mat(3, loop_isi, iSubj) = curr_HR_gender;
        
    end
end

save('HR_t_fem_mat.mat', 'HR_t_fem_mat');

%% male target face data
big_HR_t_male_mat = big_HR_mat(:, [21:40], :);

HR_t_male_mat = nan(4, 20, nSubj);
% 1 row: age
% 2 row: face task acc
% 3 row: gender task acc
% 4 row: refresh rate

for iSubj = 1:nSubj
    curr_mat = squeeze(big_HR_t_male_mat(:, :, iSubj));
    HR_t_male_mat(1, :, iSubj) = repelem(unique(curr_mat(3,:)), 20);
    HR_t_male_mat(4, :, iSubj) = repelem(unique(curr_mat(7,:)), 20);
    
    isi_vals = unique(curr_mat(4,:));
    
    % compute mean HR for each isi
    loop_isi = 0;
    for iIsi = isi_vals
        loop_isi = loop_isi + 1;
        
        % isi mask
        isi_mask = curr_mat(4, :) == iIsi;
        
        % mean HR
        curr_HR_face = mean(curr_mat(5, isi_mask));
        curr_HR_gender = mean(curr_mat(6, isi_mask));
        
        % add to data frame
        HR_t_male_mat(2, loop_isi, iSubj) = curr_HR_face;
        HR_t_male_mat(3, loop_isi, iSubj) = curr_HR_gender;
        
    end
end

save('HR_t_male_mat.mat', 'HR_t_male_mat');

%% female subject 

load('../data/HR_mat_fem.mat');

nFem = size(HR_mat_fem, 3);

HR_s_fem_mat = nan(4, 20, nFem);
% 1 row: age
% 2 row: face task acc
% 3 row: gender task acc
% 4 row: refresh rate

for iSubj = 1:nFem
    curr_mat = squeeze(HR_mat_fem(:, :, iSubj));
    HR_s_fem_mat(1, :, iSubj) = repelem(unique(curr_mat(3,:)), 20);
    HR_s_fem_mat(4, :, iSubj) = repelem(unique(curr_mat(7,:)), 20);
    
    isi_vals = unique(curr_mat(4,:));
    
    % compute mean HR for each isi
    loop_isi = 0;
    for iIsi = isi_vals
        loop_isi = loop_isi + 1;
        
        % isi mask
        isi_mask = curr_mat(4, :) == iIsi;
        
        % mean HR
        curr_HR_face = mean(curr_mat(5, isi_mask));
        curr_HR_gender = mean(curr_mat(6, isi_mask));
        
        % add to data frame
        HR_s_fem_mat(2, loop_isi, iSubj) = curr_HR_face;
        HR_s_fem_mat(3, loop_isi, iSubj) = curr_HR_gender;
        
    end
end

save('HR_s_fem_mat.mat', 'HR_s_fem_mat');

%% male subject 

load('../data/HR_mat_male.mat');

nMale = size(HR_mat_male, 3);

HR_s_male_mat = nan(4, 20, nMale);
% 1 row: age
% 2 row: face task acc
% 3 row: gender task acc
% 4 row: refresh rate

for iSubj = 1:nMale
    curr_mat = squeeze(HR_mat_male(:, :, iSubj));
    HR_s_male_mat(1, :, iSubj) = repelem(unique(curr_mat(3,:)), 20);
    HR_s_male_mat(4, :, iSubj) = repelem(unique(curr_mat(7,:)), 20);
    
    isi_vals = unique(curr_mat(4,:));
    
    % compute mean HR for each isi
    loop_isi = 0;
    for iIsi = isi_vals
        loop_isi = loop_isi + 1;
        
        % isi mask
        isi_mask = curr_mat(4, :) == iIsi;
        
        % mean HR
        curr_HR_face = mean(curr_mat(5, isi_mask));
        curr_HR_gender = mean(curr_mat(6, isi_mask));
        
        % add to data frame
        HR_s_male_mat(2, loop_isi, iSubj) = curr_HR_face;
        HR_s_male_mat(3, loop_isi, iSubj) = curr_HR_gender;
        
    end
end

save('HR_s_male_mat.mat', 'HR_s_male_mat');

%% same gender data
load('../data/big_HR_mat.mat');
load('../data/HR_mat_male.mat');

nMale = size(HR_mat_male, 3);
nSubj = size(big_HR_mat, 3);
nFem = nSubj - nMale;

HR_same_mat = nan(4, 20, nSubj);
% 1 row: age
% 2 row: face task acc
% 3 row: gender task acc
% 4 row: refresh rate

HR_same_mat(1, :, [1:nFem]) = big_HR_mat(3,[1:20],[1:nFem]);
HR_same_mat(1, :, [nFem+1:nSubj]) = big_HR_mat(3,[21:40],[nFem+1:nSubj]);

HR_same_mat(4, :, [1:nFem]) = big_HR_mat(7,[1:20],[1:nFem]);
HR_same_mat(4, :, [nFem+1:nSubj]) = big_HR_mat(7,[21:40],[nFem+1:nSubj]);

HR_same_mat([2:3], :, [1:nFem]) = big_HR_mat([5:6],[1:20],[1:nFem]);
HR_same_mat([2:3], :, [nFem+1:nSubj]) = big_HR_mat([5:6],[21:40],[nFem+1:nSubj]);

save('HR_same_mat.mat', 'HR_same_mat');

%% different gender data
load('../data/big_HR_mat.mat');
load('../data/HR_mat_male.mat');

nMale = size(HR_mat_male, 3);
nSubj = size(big_HR_mat, 3);
nFem = nSubj - nMale;

HR_diff_mat = nan(4, 20, nSubj);
% 1 row: age
% 2 row: face task acc
% 3 row: gender task acc
% 4 row: refresh rate

HR_diff_mat(1, :, [1:nFem]) = big_HR_mat(3,[21:40],[1:nFem]);
HR_diff_mat(1, :, [nFem+1:nSubj]) = big_HR_mat(3,[1:20],[nFem+1:nSubj]);

HR_diff_mat(4, :, [1:nFem]) = big_HR_mat(7,[21:40],[1:nFem]);
HR_diff_mat(4, :, [nFem+1:nSubj]) = big_HR_mat(7,[1:20],[nFem+1:nSubj]);

HR_diff_mat([2:3], :, [1:nFem]) = big_HR_mat([5:6],[21:40],[1:nFem]);
HR_diff_mat([2:3], :, [nFem+1:nSubj]) = big_HR_mat([5:6],[1:20],[nFem+1:nSubj]);

save('HR_diff_mat.mat', 'HR_diff_mat');

%% same gender data -male
load('../data/big_HR_mat.mat');
load('../data/HR_mat_male.mat');

nMale = size(HR_mat_male, 3);
nSubj = size(big_HR_mat, 3);
nFem = nSubj - nMale;

HR_same_s_male_mat = nan(4, 20, nMale);
% 1 row: age
% 2 row: face task acc
% 3 row: gender task acc
% 4 row: refresh rate

HR_same_s_male_mat(1, :, :) = big_HR_mat(3,[21:40],[nFem+1:nSubj]);
HR_same_s_male_mat(4, :, :) = big_HR_mat(7,[21:40],[nFem+1:nSubj]);
HR_same_s_male_mat([2:3], :, :) = big_HR_mat([5:6],[21:40],[nFem+1:nSubj]);

save('HR_same_s_male_mat.mat', 'HR_same_s_male_mat');

%% diff gender data -male
load('../data/big_HR_mat.mat');
load('../data/HR_mat_male.mat');

nMale = size(HR_mat_male, 3);
nSubj = size(big_HR_mat, 3);
nFem = nSubj - nMale;

HR_diff_s_male_mat = nan(4, 20, nMale);
% 1 row: age
% 2 row: face task acc
% 3 row: gender task acc
% 4 row: refresh rate

HR_diff_s_male_mat(1, :, :) = big_HR_mat(3,[1:20],[nFem+1:nSubj]);
HR_diff_s_male_mat(4, :, :) = big_HR_mat(7,[1:20],[nFem+1:nSubj]);
HR_diff_s_male_mat([2:3], :, :) = big_HR_mat([5:6],[1:20],[nFem+1:nSubj]);

save('HR_diff_s_male_mat.mat', 'HR_diff_s_male_mat');

%% same gender data -female
load('../data/big_HR_mat.mat');
load('../data/HR_mat_male.mat');

nMale = size(HR_mat_male, 3);
nSubj = size(big_HR_mat, 3);
nFem = nSubj - nMale;

HR_same_s_fem_mat = nan(4, 20, nFem);
% 1 row: age
% 2 row: face task acc
% 3 row: gender task acc
% 4 row: refresh rate

HR_same_s_fem_mat(1, :, :) = big_HR_mat(3,[1:20],[1:nFem]);
HR_same_s_fem_mat(4, :, :) = big_HR_mat(7,[1:20],[1:nFem]);
HR_same_s_fem_mat([2:3], :, :) = big_HR_mat([5:6],[1:20],[1:nFem]);

save('HR_same_s_fem_mat.mat', 'HR_same_s_fem_mat');

%% diff gender data -male
load('../data/big_HR_mat.mat');
load('../data/HR_mat_male.mat');

nMale = size(HR_mat_male, 3);
nSubj = size(big_HR_mat, 3);
nFem = nSubj - nMale;

HR_diff_s_fem_mat = nan(4, 20, nFem);
% 1 row: age
% 2 row: face task acc
% 3 row: gender task acc
% 4 row: refresh rate

HR_diff_s_fem_mat(1, :, :) = big_HR_mat(3,[21:40],[1:nFem]);
HR_diff_s_fem_mat(4, :, :) = big_HR_mat(7,[21:40],[1:nFem]);
HR_diff_s_fem_mat([2:3], :, :) = big_HR_mat([5:6],[21:40],[1:nFem]);

save('HR_diff_s_fem_mat.mat', 'HR_diff_s_fem_mat');