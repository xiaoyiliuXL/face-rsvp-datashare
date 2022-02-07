function [det_dat, fit_curve, summary] = exp_detrend(xvals, dat)

%% Function for fitting exponential function and removing from raw data
%  AT INDIVIDUAL LEVEL
    
% nSubj = size(dat,2);

%% prepare cubic fit
ft = fittype('a*exp(-x/b)+c', 'independent', 'x', 'dependent', 'y');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0];
opts.Upper = [1 Inf 1];

% define anonymous function to estimate the trend based on parameters
anon_exp = @(x, a, b, c) a .* exp(-x ./ b) + c;

% data frame for detrended data and summary
[det_dat, fit_curve] = deal(nan(size(dat)));
summary.adjR2 = nan(1, 1);;
summary.coeffs = nan(1,3);
% summary.adjR2 = nan(nsubjs, 1);
% summary.coeffs = nan(nsubjs, 3);
warning('off', 'curvefit:fit:noStartPoint')


%% do detrend

[xData, yData] = prepareCurveData(xvals, dat); % xvals: isi

% fit model to data
[fitresult, gof] = fit(xData, yData, ft, opts );

% extract coefficients
coeffs = coeffvalues(fitresult);

% evaluate trend
trend_set = anon_exp(xData, coeffs(1), coeffs(2), coeffs(3));

% subtract trend
det_dat = yData - trend_set;

% store trend for easier check
fit_curve = trend_set;
        
% store coefficients and GoF
summary.adjR2 = gof.adjrsquare;
summary.coeffs(1, :) = coeffs;


%% loop thru participants
% for iSubj = 1:nSubj
%     
%     curr_dat = dat(:, iSubj);
%     [xData, yData] = prepareCurveData(xvals, curr_dat); % xvals: isi
%     
%     % fit model to data
%     [fitresult, gof] = fit( xData, yData, ft, opts );
%     
%     % extract coefficients
%     coeffs = coeffvalues(fitresult);
%     
%     % evaluate trend
%     trend_set = anon_exp(xData, coeffs(1), coeffs(2), coeffs(3), coeffs(4));
%     
%     % subtract trend 
%     det_dat(:, isubj) = yData - trend_set;
%     
%     % store trend for easier check
%     fit_curve(:, isubj) = trend_set;
%     
%     % store coefficients and GoF
%     summary.adjR2(isubj) = gof.adjrsquare;
%     summary.coeffs(isubj, :) = coeffs;
%     
% end

end