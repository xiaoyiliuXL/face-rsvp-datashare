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
         255 51 153]/255;


%% load data files

% load('../data/HR_mat.mat');
% 1 row: age
% 2 row: face task acc
% 3 row: gender task acc
% 4 row: refresh rate

nSubj = size(big_HR_mat, 3);

gender_task = squeeze(big_HR_mat([3 4], :, :));

%% detrend at individual level
det_gender_task = nan(size(gender_task)); 
trend_gender_task = nan(size(gender_task, [2 3]));
% summary_face = []; summary_face.adjR2 = nan(nSubj, 1); summary_face.coeffs = nan(nSubj, 4);

for iSubj = 1: nSubj
    
    curr_gender_dat = squeeze(gender_task(:,:,iSubj));
    curr_fs = unique(curr_gender_dat(2,:)) / 2; % each image correspons to 2 frames
    curr_time_bin = [0:19] ./ curr_fs;
    
    [curr_det_gender, fit_curve] = exp_detrend(curr_time_bin, curr_gender_dat(1,:));
    
    det_gender_task(1,:,iSubj) = curr_det_gender'; 
    det_gender_task(2,:,iSubj) = repelem(unique(curr_gender_dat(2,:)), 20); 
    trend_gender_task(:,iSubj) = fit_curve;
%     summary_face.adjR2(iSubj) = curr_summary_face.adjR2;
%     summary_gender.adjR2(iSubj) = curr_summary_gender.adjR2;
end

% % check group detrended data
% plot(curr_time_bin, mean(squeeze(det_gender_task(1,:,:)), 2)); hold on
% plot(curr_time_bin, mean(trend_gender_task, 2)); hold on
% plot(curr_time_bin, mean(squeeze(gender_task(1,:,:)), 2))

save('./det_gender.mat', 'det_gender_task')
%load('./det_gender.mat');

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
spctr_gender = []; spctr_gender = [];

for iSubj = 1:nSubj
    
    curr_gender_dat = squeeze(det_gender_task(:,:,iSubj)); 
    curr_refresh_rate = unique(curr_gender_dat(2,:));
 
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
    curr_spctr_gender = fft_func(curr_gender_dat(1,:), params);
    
    % truncate data if long
    if size(curr_spctr_gender.spctr_out, 1) > 11
        curr_spctr_gender.cmplx_out = curr_spctr_gender.cmplx_out(1:11);
        curr_spctr_gender.spctr_out = curr_spctr_gender.spctr_out(1:11);
        curr_spctr_gender.freqs = curr_spctr_gender.freqs(1:11);
    end
   
    % combine results across participants
    if iSubj == 1
        spctr_gender.cmplx_out = nan(size(curr_spctr_gender.cmplx_out, 1), nSubj);
        spctr_gender.spctr_out = nan(size(curr_spctr_gender.spctr_out, 1), nSubj);
    end
    
    spctr_gender.cmplx_out(:, iSubj) = curr_spctr_gender.cmplx_out;
    spctr_gender.spctr_out(:, iSubj) = curr_spctr_gender.spctr_out;
    
end
spctr_gender.freqs = curr_spctr_gender.freqs;
mask_freq_interest = (unique(spctr_gender.freqs) >= min(foi)) & (unique(spctr_gender.freqs) <= max(foi));

%% see what's going on
% ampiltude
amp_gender = mean(spctr_gender.spctr_out, 2);

% phase locked amplitude
PL_gender = abs(sum(spctr_gender.cmplx_out, 2))/nSubj;

%% Permutation
[amp_mat_perm, PL_mat_perm] = deal(nan([size(spctr_gender.spctr_out),...
    n_perm]));

HR_mat_perm = label_permutation_func(det_gender_task, n_perm);

%% Do fft for each permutation
for iPerm = 1:n_perm
    
    for iSubj = 1:nSubj
        
        curr_gender_perm = HR_mat_perm(1,:,iSubj,iPerm);
        curr_refresh_rate_perm = HR_mat_perm(2,1,iSubj,iPerm);
 
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
        curr_spctr = fft_func(curr_gender_perm, params);
        
        % truncate data if long
        if size(curr_spctr.spctr_out, 1) > 11
            curr_spctr.cmplx_out = curr_spctr.cmplx_out(1:11);
            curr_spctr.spctr_out = curr_spctr.spctr_out(1:11);
            curr_spctr.freqs = curr_spctr.freqs(1:11);
        end
        
        % big mat
        amp_mat_perm(:,iSubj,iPerm) = curr_spctr.spctr_out;
        PL_mat_perm(:,iSubj,iPerm) = curr_spctr.cmplx_out;
        
    end
    
    if mod(iPerm, 100) ==0
        fprintf('\n %d permutations', iPerm)
    end
    
end

%% 95% permutation
% get group average
perm_amp_gender = squeeze(mean(amp_mat_perm(:,:,:), 2));
perm_PL_gender = squeeze(abs(sum(PL_mat_perm(:,:,:), 2)) / nSubj);

% get 95%
amp_gender_95 = prctile(perm_amp_gender, 95, 2);
PL_gender_95 = prctile(perm_PL_gender, 95, 2);

%% Plot phase locked sum
figure;
set(groot, 'defaultAxesFontSize',14)

subplot(1, 4, 1:4); hold on
plot(spctr_face.freqs, PL_face, 'Color', colors(1, :), 'LineWidth', 7)
plot(spctr_face.freqs, PL_face_95, '--','Color', colors(1, :), 'LineWidth', 4)
plot(spctr_gender.freqs, PL_gender, 'Color', colors(2, :), 'LineWidth', 7)
plot(spctr_gender.freqs, PL_gender_95, '--','Color', colors(2, :), 'LineWidth', 4)

legend('95 percentile permutations')
title('phase locked sum')
xlim(minmax(spctr_face.freqs))
xticks(spctr_face.freqs)
ylabel('amplitude (a.u.)')

figure;
ax = subplot(2, 1, 1); 
pol_ax = polaraxes('Units', ax.Units, 'Position',ax.Position);
delete(ax)
[peakfreqs_face, cntr_freq_face, pval_face, Zcirc_face] = peakfreqs(spctr_face, pol_ax, ...
    mask_freq_interest, colors(1, :));
rticks([0.05 0.1 0.15])
rticklabels({'','',''})
title(sprintf('Face detection task\n F=%0.2f Hz, Rayleigh p=%0.3f', cntr_freq_face, pval_face), 'FontSize', 14)

ax = subplot(2, 1, 2); 
pol_ax = polaraxes('Units', ax.Units, 'Position',ax.Position);
delete(ax)
[peakfreqs_gender, cntr_freq_gender, pval_gender, Zcirc_gender] = peakfreqs(spctr_gender, pol_ax, ...
    mask_freq_interest, colors(2, :));
title(sprintf('Gender discrimination task\n F=%0.2f Hz, Rayleigh p=%0.3f', cntr_freq_gender, pval_gender), 'FontSize', 14)
