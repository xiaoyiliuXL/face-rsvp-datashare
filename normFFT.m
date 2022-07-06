clear 
clc
close all

load('detrended_face.mat')
load('isi_vals.mat')


%% compute FFT
params = [];
params.detrend_flag = -1;
params.window = [];
params.power = 0;
params.zero_pad = 0;
params.subj_dim = 2;
params.f_sample = 30;
params.time_bins = isi_vals'./1000;
params.verbose = -1;
params.lp_filter = 0;

det_dat = squeeze(detrended_face(1, :, :));
FFTout = cmpt_beh_spectra(det_dat, params);

%% test pure sinusoid amplitude
puresin = .25 * sin(7.5*2*pi*params.time_bins);
FFTooutpuresin = cmpt_beh_spectra(puresin, params);

figure; 
subplot(2, 3, 1:2);
plot(isi_vals, puresin, 'k', 'LineWidth', 3)
xlabel('ISI (ms)')
ylabel('accuracy')
title({'simulated pure sinusoid', '7.5 Hz, .25 amplitude'}) 
xlim(minmax(isi_vals))

subplot(2, 3, 4:5);
plot(FFTooutpuresin.freqs, FFTooutpuresin.spctr_out*2, 'r', 'LineWidth', 3)
xlabel('frequency (Hz)')
ylabel('amplitude (\approx accuracy)')
title('FFT spectra normalized')

% fetch amplitude at 7.5 Hz
amp_theta = FFTout.spctr_out(FFTout.freqs==7.5, :);

subplot(2, 3, [3, 6])
histogram(amp_theta*2)
ylabel('N occurrences')
xlabel('amplitude (\approx accuracy)')
title({'effect size distribution at 7.5 Hz', sprintf('mean=%0.3f, median=%0.3f',...
        mean(amp_theta*2), median(amp_theta*2))})

%%
% 
% figure; hold on;
% plot(FFTout.freqs, abs(sum(FFTout.cmplx_out, 2)))
% 
% figure;
% plot(isi_vals, mean(det_dat, 2))
% 
% 
% 
