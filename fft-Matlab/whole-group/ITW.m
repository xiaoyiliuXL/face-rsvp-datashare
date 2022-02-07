%% INTEGRATION WINDOW ANALYSIS

clear all
close all
clc

load('../data/raw_data.mat');
load('big_HR_mat.mat');

ISI_frames=0:19;
f_sample = 30;
xvals = (ISI_frames/f_sample)';
str_task = {'face detection', 'gender discrimination'};

nsubjs = size(big_HR_mat, 3);

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
         255 51 153]/255;

%% GROUP DECAY FUNCTION
avg_HR_mat = mean(big_HR_mat,3)';
sderr_HR_mat = (squeeze(std(big_HR_mat,[],3))/sqrt(nsubjs))';

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



itw_mat = nan(size(avg_HR_mat,2), 2);

cvg_level = mean(avg_HR_mat([14:20],:), 1);
cvg_threshold_10 = cvg_level * 1.1;
cvg_threshold_5 = cvg_level * 1.05;

figure
cell_task_fits = cell(2,3);
plotCount = [1, 4];

ntasks = size(big_HR_mat, 1)-1;
for itask = 1:ntasks
    
    curr_avg_mat = avg_HR_mat(:,itask+1);
    
    [xData, yData] = prepareCurveData(xvals, curr_avg_mat);
    
    [cell_task_fits{itask}, gof] = fit( xData, yData, ft, opts );

    coeffs = coeffvalues(cell_task_fits{itask});
    trend_set = anon_exp(xData, coeffs(1), coeffs(2), coeffs(3));
    
    itw_10 = interp1(trend_set, xData, cvg_threshold_10(itask+1));
    itw_5 = interp1(trend_set, xData, cvg_threshold_5(itask+1));
    
    itw_mat(itask,1) = itw_10;
    itw_mat(itask,2) = itw_5;

    currPlot = plotCount(itask);
    ax = subplot(1,6,[currPlot, currPlot+1]);
    hold on
    lineProps.col = {colors(itask,:)};
    lineProps.transparent = .5;
    lineProps.style = '.';

    mseb(xData', yData', sderr_HR_mat(:,itask+1)', lineProps);

    hold on;
    scatter(xData,yData, 20, colors(itask,:), 'filled')
    hold on;
    plot(cell_task_fits{itask})
    hold on;
%     xlim(minmax(x))

%     title(str_task{itask});

    if itask<2
        
        set(ax(1),'XTickLabel','')
        set(ax(1),'XColor',get(gca,'Color'))
        set(ax(1),'box','off')
        
    else 
        
        ylabel({'decay function'});
        xlabel('time (s)')
    
        
    end

end

% %% INDIVIDUAL DECAY FUNCTION
% % itw_mat = nan(2, 2, nsubjs);
% subj_noDecay = [10 11 12 20 32 40 53 89 92 108 156 162];
% cvg_mat = nan(2,nsubjs);
% figure_idx = -1;
% ntasks = 1; %% only face task for now
% for itask = 1:ntasks
% 
%     delta_value = unique(mat_data(5,:,:));
%     curr_mat = mat_data([5 itask+5 8],:,:);
% %     curr_mat = squeeze(HR_mat(itask+1, :, :));
% 
%     for isubj = 1:nsubjs
% 
%         this_subj = squeeze(curr_mat(:,:, isubj));
%         cvg_mat(1, isubj) = unique(this_subj(3,:));
%         
%         %% individual HR mat with standard error
%         this_HR_mat = nan(2,20);
%         loopDelta = 0;
%         for iDelta = delta_value'
%             loopDelta = loopDelta + 1;
%             
%             delta_mask = this_subj(1,:) == iDelta;
%             curr_HR = mean(this_subj(2,delta_mask));
%             curr_se = std(this_subj(2,delta_mask))/sqrt(sum(delta_mask));
%             
%             this_HR_mat(1,loopDelta) = curr_HR;
%             this_HR_mat(2,loopDelta) = curr_se;
%         end
%         
%         %% individual itw
%         this_cvg_level = mean(this_HR_mat(1, 14:20));
%         cvg_mat(2,isubj) = this_cvg_level;
%         this_cvg_thre10 = this_cvg_level * 1.1;
%         this_cvg_thre5 = this_cvg_level * 1.05;
%         
%         [xData, yData] = prepareCurveData(xvals, this_HR_mat(1,:));
%         
%         % Fit model to data.
%         [fitresult, gof] = fit( xData, yData, ft, opts );
%         coeffs = coeffvalues(fitresult);
%         trend_set = anon_exp(xData, coeffs(1), coeffs(2), coeffs(3));
%         
%         if ~any(isubj == subj_noDecay)
%             itw_10 = interp1(trend_set, xData, this_cvg_thre10);
%             itw_5 = interp1(trend_set, xData, this_cvg_thre5);
%         end
% 
% %         itw_10 = interp1(trend_set, xData, this_cvg_thre10);
% %         itw_5 = interp1(trend_set, xData, this_cvg_thre5);
% 
% %         itw_mat(itask,1,isubj) = itw_10;
% %         itw_mat(itask,2,isubj) = itw_5;
%         
% %         ax = subplot(11,30,[isubj*2-1, isubj]);
% %         
% %         %% plots
% %         if mod(isubj, 10) == 1
% %             figure_idx = figure_idx + 1;
% %             figure;
% %         end
% %         
% %         plot_idx = (isubj-figure_idx*10)*2;
% %         subplot(10,2,[plot_idx-1 plot_idx]);
% %         hold on;
% %         lineProps.col = {colors(itask,:)};
% %         lineProps.transparent = .5;
% %         lineProps.style = '.';
% % 
% %         mseb(xData', yData', this_HR_mat(2,:), lineProps);
% % 
% %         hold on;
% %         scatter(xData,yData, 20, colors(itask,:), 'filled')
% %         hold on;
% %         plot(fitresult)
% %     %     xlim(minmax(x))
% % 
% %         title(sprintf('Subject %i',isubj));
% %         ylabel({'decay function'});
% %         
% %         if mod(isubj,10)==0 | isubj == 165
% %             xlabel('time (s)');
% %         end
% % 
% % %         if isubj<165
% % % 
% % %             set(ax(1),'XTickLabel','')
% % %             set(ax(1),'XColor',get(gca,'Color'))
% % %             set(ax(1),'box','off')
% % % 
% % %         else 
% % % 
% % %             xlabel('time (s)')
% % %         end
% % 
% 
%     end
% end
