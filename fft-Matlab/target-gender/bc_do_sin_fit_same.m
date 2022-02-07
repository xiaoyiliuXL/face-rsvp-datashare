%% Sinusoidal fitting analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;

% load data
load('./det_face_task_same.mat');
nSubj = size(det_face_task_same, 3);

% settings
n_perm = 10000;
ISI_frames = 0:19;
Fs = 30;
time_bins = (ISI_frames/Fs)';

% graph settings
rng(0);
colors = [0 204 204;
         255 51 153]/255;
     
%% get group mean and SD
mean_HR_det_same = squeeze(mean(det_face_task_same(1,:,:), 3));
sderr_HR_det_same = squeeze(std(det_face_task_same(1,:,:), [], 3) / sqrt(nSubj));

%% do sin fit
ft = fittype('sin1');
x = time_bins';

adj_R2 = nan(1);

% do fit
[xData, yData] = prepareCurveData(x, mean_HR_det_same);
opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Robust',...
        'Bisquare');
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf]; % -Inf lowerBands(iLoad,2) -Inf];
opts.Upper = [+Inf +Inf +Inf]; % +Inf upperBands(iLoad,2) +Inf];
% opts.StartPoint = [0.006965, 47.14,  -1.329];


[fitresult, statR, out] = fit(xData, yData, ft, opts);

%% plot
figure

ax = subplot(1, 3, [1, 2]);
hold on
% sd band
lineProps.col = {colors(1,:)};
lineProps.transparent = .5;
lineProps.style = '.';

mseb(xData', yData', sderr_HR_det_same, lineProps);

hold on
% mean dots
scatter(xData,yData, 50, colors(1,:), 'filled');

hold on
% fit curve
fplot(@(x) fitresult(x), [0 0.633], 'Color', colors(1,:), 'LineWidth',7)
xlim(minmax(x))

% title
swap = coeffvalues(fitresult);
freqs_fitted = swap(2)/(2*pi);
str_title1 = ['Best frequency sinusoid fit: ' ...
        num2str(round(freqs_fitted(1),2)) ' Hz'];
str_title2 = ['AdjR^2: ' num2str(round(statR.adjrsquare,2))];
    
real_adj_R2 = statR.adjrsquare;

title({'face detection', str_title1, str_title2});
    
ylabel({'2nd order', 'detrended HR'});

%% prepare permutation
shuffled_HR_det_same = label_permutation_func(det_face_task_same, n_perm);

mean_HR_det_shuffled_same = squeeze(mean(shuffled_HR_det_same(1,:,:,:), 3));

%% do perm

for iPerm = 1:n_perm
    
    [xData, yData] = prepareCurveData(x, mean_HR_det_shuffled_same(:, iPerm));
    
    [~, statP, ~] = fit(xData, yData , ft, opts);
    
    permuted_adj_R2(iPerm) = statP.adjrsquare;
    
    if mod(iPerm, 100) == 0
        fprintf('%d permuted fittings \n', iPerm)
    end
    
end

%% plot perm
ax = subplot(1,3,3);

h = histogram(permuted_adj_R2,100, 'FaceColor', [.1 .1 .1],...
        'Normalization', 'probability');
    
hold on
lineUp = max(h.Values); 

% compute empirical cdf on the permutations 
[curr_cdf(:,1), curr_cdf(:,2)] = ecdf(permuted_adj_R2);
% reverse the vector and find the first value exceeding alpha (.05)
curr_threshold = curr_cdf(find(abs(curr_cdf(:,1)-1)<.05,1),2);

% plot 95° percentile (set preferred value on y axis to make the line
% coherent with graph
plot([curr_threshold curr_threshold], [0 lineUp], '--k', 'LineWidth',4)

hold on

% real value obtained
curr_real_R2 = real_adj_R2;

% plot real value obtained
plot([curr_real_R2, curr_real_R2], [0 lineUp],'Color', colors(1,:),...
    'LineWidth', 4);

% see where your p is
pVal = 1-curr_cdf(find(curr_cdf(:,2)>=curr_real_R2,1)-1,1);

% add label for p value
text(curr_real_R2, lineUp-lineUp/10, ['p = ' num2str(pVal)])


    set(ax(1),'XTickLabel','')
    set(ax(1),'XColor',get(gca,'Color'))
    set(ax(1),'box','off')

    xlabel('Adj R^2')
    set(ax(1),'box','off')



    title('permuted AdjR^2')

ylim([0 lineUp])
xlim([-.3 1])
ylabel('probability')

