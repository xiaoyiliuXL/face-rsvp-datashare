%% FFT analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;

addpath /Users/xl4251/Documents/MATLAB/CircStat2012a;

% FFT settings
n_perm = 10000;
foi = [0, 15]; % frequency of interest from 0 to 15

% graph settings
rng(0);
colors = [0 204 204;
         255 51 153]/255;
     
%% load data files

load('HR_same_s_fem_mat.mat');
% 1 row: age
% 2 row: face task acc
% 3 row: gender task acc
% 4 row: refresh rate

nFem = size(HR_same_s_fem_mat, 3);

face_task_same_s_fem = squeeze(HR_same_s_fem_mat([2 4], :, :));

%% detrend at individual level
det_face_task_same_s_fem = nan(size(face_task_same_s_fem)); 
trend_face_task_same_s_fem = nan(size(face_task_same_s_fem, [2 3]));
% summary_face = []; summary_face.adjR2 = nan(nSubj, 1); summary_face.coeffs = nan(nSubj, 4);

for iSubj = 1: nFem
    
    curr_face_dat = squeeze(face_task_same_s_fem(:,:,iSubj));
    curr_fs = unique(curr_face_dat(2,:)) / 2; % each image correspons to 2 frames
    curr_time_bin = [0:19] ./ curr_fs;
    
    [curr_det_face, fit_curve] = exp_detrend(curr_time_bin, curr_face_dat(1,:));
    
    det_face_task_same_s_fem(1,:,iSubj) = curr_det_face'; 
    det_face_task_same_s_fem(2,:,iSubj) = repelem(unique(curr_face_dat(2,:)), 20); 
    trend_face_task_same_s_fem(:,iSubj) = fit_curve;
%     summary_face.adjR2(iSubj) = curr_summary_face.adjR2;
%     summary_gender.adjR2(iSubj) = curr_summary_gender.adjR2;
end

% check group detrended data
% plot(curr_time_bin, mean(squeeze(det_face_task_s_fem(1,:,:)), 2)); hold on
% plot(curr_time_bin, mean(trend_face_task_s_fem, 2)); hold on
% plot(curr_time_bin, mean(squeeze(face_task_s_fem(1,:,:)), 2))

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
spctr_face_same_s_fem = [];

for iSubj = 1:nFem
    
    curr_face_dat = squeeze(det_face_task_same_s_fem(:,:,iSubj)); 
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
        spctr_face_same_s_fem.cmplx_out = nan(size(curr_spctr_face.cmplx_out, 1), nFem);
        spctr_face_same_s_fem.spctr_out = nan(size(curr_spctr_face.spctr_out, 1), nFem);
    end
    
    spctr_face_same_s_fem.cmplx_out(:, iSubj) = curr_spctr_face.cmplx_out;
    spctr_face_same_s_fem.spctr_out(:, iSubj) = curr_spctr_face.spctr_out;
    
end
spctr_face_same_s_fem.freqs = curr_spctr_face.freqs;
mask_freq_interest = (unique(spctr_face_same_s_fem.freqs) >= min(foi)) & (unique(spctr_face_same_s_fem.freqs) <= max(foi));

%% see what's going on
% ampiltude
amp_face_same_s_fem = mean(spctr_face_same_s_fem.spctr_out, 2);

% phase locked amplitude
PL_face_same_s_fem = abs(sum(spctr_face_same_s_fem.cmplx_out, 2))/nFem;

%% Permutation
[amp_mat_perm_same_s_fem, PL_mat_perm_same_s_fem] = deal(nan([size(spctr_face_same_s_fem.spctr_out),...
    n_perm]));

HR_mat_perm_same_s_fem = label_permutation_func(det_face_task_same_s_fem, n_perm);

%% Do fft for each permutation
for iPerm = 1:n_perm
    
    for iSubj = 1:nFem
        
        curr_face_perm = HR_mat_perm_same_s_fem(1,:,iSubj,iPerm);
        curr_refresh_rate_perm = HR_mat_perm_same_s_fem(2,1,iSubj,iPerm);
 
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
        amp_mat_perm_same_s_fem(:,iSubj,iPerm) = curr_spctr.spctr_out;
        PL_mat_perm_same_s_fem(:,iSubj,iPerm) = curr_spctr.cmplx_out;
        
    end
    
    if mod(iPerm, 100) ==0
        fprintf('\n %d permutations', iPerm)
    end
    
end

%% 95% permutation
% get group average
perm_amp_face_same_s_fem = squeeze(mean(amp_mat_perm_same_s_fem(:,:,:), 2));
perm_PL_face_same_s_fem = squeeze(abs(sum(PL_mat_perm_same_s_fem(:,:,:), 2)) / nFem);

% get 95%
amp_face_95_same_s_fem = prctile(perm_amp_face_same_s_fem, 95, 2);
PL_face_95_same_s_fem = prctile(perm_PL_face_same_s_fem, 95, 2);

%% Plot phase locked sum
set(groot, 'defaultAxesFontSize',14)

figure; 

subplot(2, 3, 1:2); hold on
plot(unique(spctr_face_same_s_fem.freqs), PL_face_same_s_fem, 'Color', colors(1, :), 'LineWidth', 7)
plot(unique(spctr_face_same_s_fem.freqs), PL_face_95_same_s_fem, 'r--', 'LineWidth', 4)
legend('face detection task', '95 percentile permutations')
title('phase locked sum')
% xlim(minmax(unique(spctr_face.freqs)))
ylabel('amplitude (a.u.)')

ax = subplot(2, 3, 3); 
pol_ax = polaraxes('Units', ax.Units, 'Position',ax.Position);
delete(ax)
[peakfreqs_face, cntr_freq_face, pval_face, Zcirc_face] = peakfreqs(spctr_face_same_s_fem, pol_ax, ...
    mask_freq_interest, colors(1, :));
title(sprintf('Face detection task\n F=%0.2f Hz, Rayleigh p=%0.3f', cntr_freq_face, pval_face), 'FontSize', 14)

%% p value
pvals_face_same_s_fem = compute_pvals_and_mc(PL_face_same_s_fem, ...
    perm_PL_face_same_s_fem, ...
    spctr_face_same_s_fem.freqs)

