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
    
%% female target
nFem = sum(squeeze(mat_data(2,1,:)) == 0);
nMale = sum(~squeeze(mat_data(2,1,:)) == 0);
nSubj = nFem + nMale;

HR_mat_t_fem = nan(4, 20, nSubj); 

% participant loop

for iSubj = 1:nSubj
    
    % get current subject's data
    curr_mat = mat_data(:,:,iSubj);
    isi_vals = unique(curr_mat(5,:));
    curr_age = unique(curr_mat(8,:));
    
    curr_refresh_rate = unique(curr_mat(9,:));
    
    t_fem_mask = curr_mat(1,:) == 0;
    
    %% get average HR
    loopIsi = 0;
    
    % ISI loops
    for iIsi = isi_vals
        
        loopIsi = loopIsi + 1;
        
        % current isi mask
        isi_mask = curr_mat(5,:) == iIsi;
        
        % compute HR at each ISI condition
        curr_HR_face = nanmean(curr_mat(6,isi_mask & t_fem_mask));
        curr_HR_gender = nanmean(curr_mat(7,isi_mask & t_fem_mask));
        
        % put the average HR in current ISI condition in place
        HR_mat_t_fem(2, loopIsi, iSubj) = curr_HR_face;
        HR_mat_t_fem(3, loopIsi, iSubj) = curr_HR_gender;
        
    end
    
    % participant loop
    HR_mat_t_fem(1, :, iSubj) = repelem(curr_age, 20);
    HR_mat_t_fem(4, :, iSubj) = repelem(curr_refresh_rate, 20);
end

save('HR_mat_t_fem.mat', 'HR_mat_t_fem')

%% male target
nFem = sum(squeeze(mat_data(2,1,:)) == 0);
nMale = sum(~squeeze(mat_data(2,1,:)) == 0);
nSubj = nFem + nMale;

HR_mat_t_male = nan(4, 20, nSubj); 

% participant loop

for iSubj = 1:nSubj
    
    % get current subject's data
    curr_mat = mat_data(:,:,iSubj);
    isi_vals = unique(curr_mat(5,:));
    curr_age = unique(curr_mat(8,:));
    
    curr_refresh_rate = unique(curr_mat(9,:));
    
    t_male_mask = curr_mat(1,:) == 1;
    
    %% get average HR
    loopIsi = 0;
    
    % ISI loops
    for iIsi = isi_vals
        
        loopIsi = loopIsi + 1;
        
        % current isi mask
        isi_mask = curr_mat(5,:) == iIsi;
        
        % compute HR at each ISI condition
        curr_HR_face = nanmean(curr_mat(6,isi_mask & t_male_mask));
        curr_HR_gender = nanmean(curr_mat(7,isi_mask & t_male_mask));
        
        % put the average HR in current ISI condition in place
        HR_mat_t_male(2, loopIsi, iSubj) = curr_HR_face;
        HR_mat_t_male(3, loopIsi, iSubj) = curr_HR_gender;
        
    end
    
    % participant loop
    HR_mat_t_male(1, :, iSubj) = repelem(curr_age, 20);
    HR_mat_t_male(4, :, iSubj) = repelem(curr_refresh_rate, 20);
end

save('HR_mat_t_male.mat', 'HR_mat_t_male')

%% male target male subject
nFem = sum(squeeze(mat_data(2,1,:)) == 0);
nMale = sum(~squeeze(mat_data(2,1,:)) == 0);
nSubj = nFem + nMale;

HR_same_s_male_mat = nan(4, 20, nMale); 

% participant loop
iMale = 0;

for iSubj = 1:nSubj
    
    % get current subject's data
    curr_mat = mat_data(:,:,iSubj);
    curr_gender = curr_mat(2,:);
    
    if curr_gender == 1
        
        iMale = iMale + 1;
        
        isi_vals = unique(curr_mat(5,:));
        curr_age = unique(curr_mat(8,:));

        curr_refresh_rate = unique(curr_mat(9,:));

        t_male_mask = curr_mat(1,:) == 1;


        %% get average HR
        loopIsi = 0;

        % ISI loops
        for iIsi = isi_vals

            loopIsi = loopIsi + 1;

            % current isi mask
            isi_mask = curr_mat(5,:) == iIsi;

            % compute HR at each ISI condition
            curr_HR_face = nanmean(curr_mat(6,isi_mask & t_male_mask));
            curr_HR_gender = nanmean(curr_mat(7,isi_mask & t_male_mask));

            % put the average HR in current ISI condition in place
            HR_same_s_male_mat(2, loopIsi, iMale) = curr_HR_face;
            HR_same_s_male_mat(3, loopIsi, iMale) = curr_HR_gender;

        end

        % participant loop
        HR_same_s_male_mat(1, :, iMale) = repelem(curr_age, 20);
        HR_same_s_male_mat(4, :, iMale) = repelem(curr_refresh_rate, 20);
    end
    
end

save('HR_same_s_male_mat.mat', 'HR_same_s_male_mat')

%% female target male subject
nFem = sum(squeeze(mat_data(2,1,:)) == 0);
nMale = sum(~squeeze(mat_data(2,1,:)) == 0);
nSubj = nFem + nMale;

HR_diff_s_male_mat = nan(4, 20, nMale); 

% participant loop
iMale = 0;

for iSubj = 1:nSubj
    
    % get current subject's data
    curr_mat = mat_data(:,:,iSubj);
    curr_gender = curr_mat(2,:);
    
    if curr_gender == 1
        
        iMale = iMale + 1;
        
        isi_vals = unique(curr_mat(5,:));
        curr_age = unique(curr_mat(8,:));

        curr_refresh_rate = unique(curr_mat(9,:));

        t_female_mask = curr_mat(1,:) == 0;


        %% get average HR
        loopIsi = 0;

        % ISI loops
        for iIsi = isi_vals

            loopIsi = loopIsi + 1;

            % current isi mask
            isi_mask = curr_mat(5,:) == iIsi;

            % compute HR at each ISI condition
            curr_HR_face = nanmean(curr_mat(6,isi_mask & t_female_mask));
            curr_HR_gender = nanmean(curr_mat(7,isi_mask & t_female_mask));

            % put the average HR in current ISI condition in place
            HR_diff_s_male_mat(2, loopIsi, iMale) = curr_HR_face;
            HR_diff_s_male_mat(3, loopIsi, iMale) = curr_HR_gender;

        end

        % participant loop
        HR_diff_s_male_mat(1, :, iMale) = repelem(curr_age, 20);
        HR_diff_s_male_mat(4, :, iMale) = repelem(curr_refresh_rate, 20);
    end
    
end

save('HR_diff_s_male_mat.mat', 'HR_diff_s_male_mat')

%% female target female subject
nFem = sum(squeeze(mat_data(2,1,:)) == 0);
nMale = sum(~squeeze(mat_data(2,1,:)) == 0);
nSubj = nFem + nMale;

HR_same_s_fem_mat = nan(4, 20, nFem); 

% participant loop
iFem = 0;

for iSubj = 1:nSubj
    
    % get current subject's data
    curr_mat = mat_data(:,:,iSubj);
    curr_gender = curr_mat(2,:);
    
    if curr_gender == 0
        
        iFem = iFem + 1;
        
        isi_vals = unique(curr_mat(5,:));
        curr_age = unique(curr_mat(8,:));

        curr_refresh_rate = unique(curr_mat(9,:));

        t_female_mask = curr_mat(1,:) == 0;


        %% get average HR
        loopIsi = 0;

        % ISI loops
        for iIsi = isi_vals

            loopIsi = loopIsi + 1;

            % current isi mask
            isi_mask = curr_mat(5,:) == iIsi;

            % compute HR at each ISI condition
            curr_HR_face = nanmean(curr_mat(6,isi_mask & t_female_mask));
            curr_HR_gender = nanmean(curr_mat(7,isi_mask & t_female_mask));

            % put the average HR in current ISI condition in place
            HR_same_s_fem_mat(2, loopIsi, iFem) = curr_HR_face;
            HR_same_s_fem_mat(3, loopIsi, iFem) = curr_HR_gender;

        end

        % participant loop
        HR_same_s_fem_mat(1, :, iFem) = repelem(curr_age, 20);
        HR_same_s_fem_mat(4, :, iFem) = repelem(curr_refresh_rate, 20);
    end
    
end

save('HR_same_s_fem_mat.mat', 'HR_same_s_fem_mat')



%% male target female subject
nFem = sum(squeeze(mat_data(2,1,:)) == 0);
nMale = sum(~squeeze(mat_data(2,1,:)) == 0);
nSubj = nFem + nMale;

HR_diff_s_fem_mat = nan(4, 20, nFem); 

% participant loop
iFem = 0;

for iSubj = 1:nSubj
    
    % get current subject's data
    curr_mat = mat_data(:,:,iSubj);
    curr_gender = curr_mat(2,:);
    
    if curr_gender == 0
        
        iFem = iFem + 1;
        
        isi_vals = unique(curr_mat(5,:));
        curr_age = unique(curr_mat(8,:));

        curr_refresh_rate = unique(curr_mat(9,:));

        t_male_mask = curr_mat(1,:) == 1;


        %% get average HR
        loopIsi = 0;

        % ISI loops
        for iIsi = isi_vals

            loopIsi = loopIsi + 1;

            % current isi mask
            isi_mask = curr_mat(5,:) == iIsi;

            % compute HR at each ISI condition
            curr_HR_face = nanmean(curr_mat(6,isi_mask & t_male_mask));
            curr_HR_gender = nanmean(curr_mat(7,isi_mask & t_male_mask));

            % put the average HR in current ISI condition in place
            HR_diff_s_fem_mat(2, loopIsi, iFem) = curr_HR_face;
            HR_diff_s_fem_mat(3, loopIsi, iFem) = curr_HR_gender;

        end

        % participant loop
        HR_diff_s_fem_mat(1, :, iFem) = repelem(curr_age, 20);
        HR_diff_s_fem_mat(4, :, iFem) = repelem(curr_refresh_rate, 20);
    end
    
end

save('HR_diff_s_fem_mat.mat', 'HR_diff_s_fem_mat')

