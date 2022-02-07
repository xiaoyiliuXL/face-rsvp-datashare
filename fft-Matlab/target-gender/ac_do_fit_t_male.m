%% FFT analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all;
% clc;
% 
% addpath /Users/xl4251/Documents/MATLAB/CircStat2012a;

% FFT settings
n_perm = 10000;
foi = [0, 15]; % frequency of interest from 0 to 15

% graph settings
rng(0);
colors = [0 204 204;
          0 104 204;
         255 51 153;
         255 51 53]/255;

%% load data files

load('HR_t_male_mat.mat');
% 1 row: age
% 2 row: face task acc
% 3 row: gender task acc
% 4 row: refresh rate

nSubj = size(HR_t_male_mat, 3);

face_task_t_male = squeeze(HR_t_male_mat([2 4], :, :));

%% detrend at individual level
det_face_task_t_male = nan(size(face_task_t_male)); 
trend_face_task_t_male = nan(size(face_task_t_male, [2 3]));
% summary_face = []; summary_face.adjR2 = nan(nSubj, 1); summary_face.coeffs = nan(nSubj, 4);

for iSubj = 1: nSubj
    
    curr_face_dat = squeeze(face_task_t_male(:,:,iSubj));
    curr_fs = unique(curr_face_dat(2,:)) / 2; % each image correspons to 2 frames
    curr_time_bin = [0:19] ./ curr_fs;
    
    [curr_det_face, fit_curve] = exp_detrend(curr_time_bin, curr_face_dat(1,:));
    
    det_face_task_t_male(1,:,iSubj) = curr_det_face'; 
    det_face_task_t_male(2,:,iSubj) = repelem(unique(curr_face_dat(2,:)), 20); 
    trend_face_task_t_male(:,iSubj) = fit_curve;
%     summary_face.adjR2(iSubj) = curr_summary_face.adjR2;
%     summary_gender.adjR2(iSubj) = curr_summary_gender.adjR2;
end

% check group detrended data
% figure;
% plot(curr_time_bin, mean(squeeze(det_face_task_t_male(1,:,:)), 2)); hold on
% plot(curr_time_bin, mean(trend_face_task_t_male, 2)); hold on
% plot(curr_time_bin, mean(squeeze(face_task_t_male(1,:,:)), 2))

save('./det_face_task_t_male.mat', 'det_face_task_t_male')

%% Create required parameters for fft
isi_frames = 0:19;

params = [];
params.detrend_flag = -1; %don't do detrend (already done)
params.window = [];
params.power = 0;
params.subj_dim = 1;
params.verbose = -1;
params.lp_filter = 0;

%% Do fft
spctr_face_t_male = []; spctr_gender_t_male = [];

for iSubj = 1:nSubj
    
    curr_face_dat = squeeze(det_face_task_t_male(:,:,iSubj)); 
    curr_refresh_rate = unique(curr_face_dat(2,:));
 
%     if curr_refresh_rate == 75
%         params.f_sample = (curr_refresh_rate+3) / 2;
%     else
%         params.f_sample = curr_refresh_rate / 2;
%     end
 
    params.f_sample = 30;
    
    params.time_bins = isi_frames ./ params.f_sample;
    
    % zero pad parameter
    data_length_needed = params.f_sample/1.5;
    params.zero_pad = (data_length_needed - 20)/2;
    
    % do fft
    curr_spctr_face = fft_func(curr_face_dat(1,:), params);
    
    % truncate data if long
    if size(curr_spctr_face.spctr_out, 1) > 11
        curr_spctr_face.cmplx_out = curr_spctr_face.cmplx_out(1:11);
        curr_spctr_face.spctr_out = curr_spctr_face.spctr_out(1:11);
        curr_spctr_face.freqs = curr_spctr_face.freqs(1:11);
    end
   
    % combine results across participants
    if iSubj == 1
        spctr_face_t_male.cmplx_out = nan(size(curr_spctr_face.cmplx_out, 1), nSubj);
        spctr_face_t_male.spctr_out = nan(size(curr_spctr_face.spctr_out, 1), nSubj);
    end
    
    spctr_face_t_male.cmplx_out(:, iSubj) = curr_spctr_face.cmplx_out;
    spctr_face_t_male.spctr_out(:, iSubj) = curr_spctr_face.spctr_out;
    
end
spctr_face_t_male.freqs = curr_spctr_face.freqs;
mask_freq_interest = (unique(spctr_face_t_male.freqs) >= min(foi)) & (unique(spctr_face_t_male.freqs) <= max(foi));

%% see what's going on
% ampiltude
amp_face_t_male = mean(spctr_face_t_male.spctr_out, 2);

% phase locked amplitude
PL_face_t_male = abs(sum(spctr_face_t_male.cmplx_out, 2))/nSubj;

%% Permutation
[amp_mat_perm_t_male, PL_mat_perm_t_male] = deal(nan([size(spctr_face_t_male.spctr_out),...
    n_perm]));

HR_mat_perm_t_male = label_permutation_func(det_face_task_t_male, n_perm);

%% Do fft for each permutation
for iPerm = 1:n_perm
    
    for iSubj = 1:nSubj
        
        curr_face_perm = HR_mat_perm_t_male(1,:,iSubj,iPerm);
        curr_refresh_rate_perm = HR_mat_perm_t_male(2,1,iSubj,iPerm);
 
%         if curr_refresh_rate_perm == 75
%             params.f_sample = (curr_refresh_rate_perm+3) / 2;
%         else
%             params.f_sample = curr_refresh_rate_perm / 2;
%         end
        
        params.f_sample = 30;
        
        params.time_bins = isi_frames ./ params.f_sample;
        
        % zero pad parameter
        data_length_needed = params.f_sample/1.5;
        params.zero_pad = (data_length_needed - 20)/2;
        
        % do fft
        curr_spctr = fft_func(curr_face_perm, params);
        
        % truncate data if long
        if size(curr_spctr.spctr_out, 1) > 11
            curr_spctr.cmplx_out = curr_spctr.cmplx_out(1:11);
            curr_spctr.spctr_out = curr_spctr.spctr_out(1:11);
            curr_spctr.freqs = curr_spctr.freqs(1:11);
        end
        
        % big mat
        amp_mat_perm_t_male(:,iSubj,iPerm) = curr_spctr.spctr_out;
        PL_mat_perm_t_male(:,iSubj,iPerm) = curr_spctr.cmplx_out;
        
    end
    
    if mod(iPerm, 100) ==0
        fprintf('\n %d permutations', iPerm)
    end
    
end

%% 95% permutation
% get group average
perm_amp_face_t_male = squeeze(mean(amp_mat_perm_t_male(:,:,:), 2));
perm_PL_face_t_male = squeeze(abs(sum(PL_mat_perm_t_male(:,:,:), 2)) / nSubj);

% get 95%
amp_face_95_t_male = prctile(perm_amp_face_t_male, 95, 2);
PL_face_95_t_male = prctile(perm_PL_face_t_male, 95, 2);

%% Plot phase locked sum
set(groot, 'defaultAxesFontSize',14)

figure; 

subplot(1, 4, 1:4); hold on
plot(unique(spctr_face_t_male.freqs), PL_face_t_male, 'Color', colors(2, :), 'LineWidth', 7)
plot(unique(spctr_face_t_male.freqs), PL_face_95_t_male, '--','Color', colors(2, :), 'LineWidth', 4)
plot(unique(spctr_face_t_fem.freqs), PL_face_t_fem, 'Color', colors(1, :), 'LineWidth', 7)
plot(unique(spctr_face_t_fem.freqs), PL_face_95_t_fem, '--','Color', colors(1, :), 'LineWidth', 4)

legend('95 percentile permutations')
title('phase locked sum')
xlim(minmax(spctr_face_t_male.freqs))
xticks(spctr_face_t_male.freqs)
ylabel('amplitude (a.u.)')

figure;
ax = subplot(2, 1, 1); 
pol_ax = polaraxes('Units', ax.Units, 'Position',ax.Position);
delete(ax)
[peakfreqs_face_t_male, cntr_freq_face_t_male, pval_face_t_male, Zcirc_face_t_male] = peakfreqs(spctr_face_t_male, pol_ax, ...
    mask_freq_interest, colors(2, :));
title(sprintf('Face detection task\n F=%0.2f Hz, Rayleigh p=%0.3f', cntr_freq_face_t_male, pval_face_t_male), 'FontSize', 14)

ax = subplot(2, 1, 2); 
pol_ax = polaraxes('Units', ax.Units, 'Position',ax.Position);
delete(ax)
[peakfreqs_face_t_fem, cntr_freq_face_t_fem, pval_face_t_fem, Zcirc_face_t_fem] = peakfreqs(spctr_face_t_fem, pol_ax, ...
    mask_freq_interest, colors(1, :));
title(sprintf('Face detection task\n F=%0.2f Hz, Rayleigh p=%0.3f', cntr_freq_face_t_fem, pval_face_t_fem), 'FontSize', 14)

%% p value
pvals_face_t_male = compute_pvals_and_mc(PL_face_t_male, ...
    perm_PL_face_t_male, ...
    spctr_face_t_male.freqs)

