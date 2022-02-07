%% Time domain figures

clear all
close all
clc

%% PREPARE FOR EXPONENTIAL FIT
load('../data/raw_data.mat');
load('big_HR_mat.mat');

ISI_frames=0:19;
f_sample = 30;
xvals = (ISI_frames/f_sample)';
str_task = {'face detection', 'gender discrimination'};

nSubj = size(big_HR_mat, 3);

% prepare fit
ft = fittype( 'a*exp(-x/b)+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0]; 
opts.Upper = [1 Inf 1];
% opts.StartPoint = [0.2614, 1.5340, 0.5217]; % parameter based on a first estimation on grand average for face detection
                                            % the choice of this srtarting
                                            % points does not inficiate the
                                            % validity of results since the
                                            % detrend is applied before the
                                            % permutations.

% define anonymous function to estimate the trend based on parameters
anon_exp = @(x, a, b, c) a .* exp(-x ./ b) + c;

% plot color
colors = [0 204 204;
          0 104 204;
         255 51 153;
         255 51 53]/255;


%% GROUP DECAY FUNCTION
avg_HR_mat = mean(big_HR_mat,3)';
sderr_HR_mat = (squeeze(std(big_HR_mat,[],3))/sqrt(nSubj))';

%% FACE
face_avg_mat = avg_HR_mat(:,2);
face_sderr_mat = sderr_HR_mat(:,2);

[xData_face, yData_face] = prepareCurveData(xvals, face_avg_mat);
[face_fits, gof_face] = fit(xData_face, yData_face, ft, opts );
coeffs_face = coeffvalues(face_fits);
trend_face = anon_exp(xData_face, coeffs_face(1), coeffs_face(2), coeffs_face(3));

det_face_task_sum = face_avg_mat-trend_face;
save('det_face_task_sum.mat', 'det_face_task_sum');

%% GENDER
gender_avg_mat = avg_HR_mat(:,3);
gender_sderr_mat = sderr_HR_mat(:,3);

[xData_gender, yData_gender] = prepareCurveData(xvals, gender_avg_mat);
[gender_fits, gof_gender] = fit(xData_gender, yData_gender, ft, opts );
coeffs_gender = coeffvalues(gender_fits);
trend_gender = anon_exp(xData_face, coeffs_gender(1), coeffs_gender(2), coeffs_gender(3));

det_gender_task_sum = gender_avg_mat-trend_gender;
save('det_gender_task_sum.mat', 'det_gender_task_sum');

%% FIGURE
figure(1)

ax = subplot(1,6, [1:6]);
hold on;
lineProps_face.col = {colors(1,:)};
lineProps_face.transparent = .7;
lineProps_face.style = '.';
mseb(xData_face', yData_face', face_sderr_mat(:,1)', lineProps_face);

hold on;
scatter(xData_face, yData_face, 50, colors(1,:), 'filled');

hold on;
plot(xvals, face_fits(xvals), 'Color', colors(1,:), 'LineWidth',3)
yline(0.431, '--', 'Color', colors(1,:), 'LineWidth',2)


hold on;
lineProps_gender.col = {colors(2,:)};
lineProps_gender.transparent = .7;
lineProps_gender.style = '.';
mseb(xData_gender', yData_gender', gender_sderr_mat(:,1)', lineProps_gender);

hold on;
scatter(xData_gender, yData_gender, 50, colors(2,:), 'filled');

hold on;
plot(xvals, gender_fits(xvals), 'Color', colors(2,:), 'LineWidth',3)
yline(0.564, '--', 'Color', colors(2,:), 'LineWidth',2)
yline(0.5, '--','Color', 'k', 'LineWidth',2)
ylabel('Hit rate');
xlabel('time (s)')


%% DETRENDED FIGURE OF FEMALE AND MALE TARGETS
load('det_face_task_t_fem.mat');
load('det_face_task_t_male.mat');

nSubj = size(det_face_task_t_fem, 3);

det_face_t_fem_avg_mat = mean(det_face_task_t_fem(1,:,:),3)';
det_face_t_fem_sderr_mat = (squeeze(std(det_face_task_t_fem(1,:,:),[],3))/sqrt(nSubj))';

det_face_t_male_avg_mat = mean(det_face_task_t_male(1,:,:),3)';
det_face_t_male_sderr_mat = (squeeze(std(det_face_task_t_fem(1,:,:),[],3))/sqrt(nSubj))';

figure(2)
ax = subplot(1,4, [1:4]);
hold on;
lineProps_face.col = {colors(1,:)};
lineProps_face.transparent = .01;
lineProps_face.style = '.';
mseb(xvals', det_face_t_fem_avg_mat, det_face_t_fem_sderr_mat', lineProps_face, 0.7);

hold on;
plot(xvals, det_face_t_fem_avg_mat, 'Color', colors(1,:), 'LineWidth',7)

hold on;
lineProps_gender.col = {colors(2,:)};
lineProps_gender.transparent = .01;
lineProps_gender.style = '.';
mseb(xvals', det_face_t_male_avg_mat, det_face_t_male_sderr_mat', lineProps_gender, 0.7);

hold on;
plot(xvals, det_face_t_male_avg_mat, 'Color', colors(2,:), 'LineWidth',7)


%% phase locked sum (run call_spectra_function.m first)
set(groot, 'defaultAxesFontSize',14)

figure; 

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
