function [outStruct] = fft_func(data, params)

%% fft function
% 
%   data   <- array of HR, if detection task, or accuracy if 2AFC
%   params <- structure containing following fields:
%
%   .detrend_flag: apply detrend on data, to be used consciously. Possible
%       values ranging from -1 (no detrend) up to whatever value you want
%       indicating the order of the detrend itself.
%
%   .window: "hanning', 'hamming', 'tukey' or '[]' (none)
%
%   .power: [1 0] take the squared value of fft (or not, simple
%       amplitude)
%
%   .zero_pad: scalar, number of 0s to be applied to both signal tails (if
%       0 does nothing)
%
%   .subj_dim: [1 2] if data is a matrix (implying subjects on one dimension
%       and evolution of hit rate on the other, SPECIFY WHICH DIMENSION IS THE
%       ONE OF THE SUBJECTS!!!!!!!!
%
%   .time_bins: vect of timing for each point collected
%
%   .f_sample: sampling frequency
%
%   .lp_filter: [0 1]
%
%   .verbose: degree of verbosity (-1 0 1): 0 to suppress text output, -1
%   to suppress warnings as well. Don't suppress them if it's the first time
%   you use the function, just to have a sanity check on your data.
%   

        exp_fields = {'detrend_flag', 'hanning_flag','power', 'zero_pad',...
            'subj_dim', 'time_bins', 'f_sample', 'lp_filter', 'verbose'};

        %% test input of function
        f_num = numel(fieldnames(params));
        if f_num<9
            field_miss = exp_fields{isfield(params, exp_fields)};
            error(['the field(s): ' field_miss ' from params is missing'])
        end

        %% put data in the right way if necessary
        if params.subj_dim == 1
            n_subj = size(data,1);
            n_time = size(data,2);
            if params.verbose ==1
                fprintf('\n\n##########################################################\n\n')
                fprintf('n subj: %d \nn timepoints: %d\n', n_subj, n_time)
                fprintf('\n\n##########################################################\n\n')
            end
            data = data'; % this will return in an error if someone tries a ND matrix
        elseif params.subj_dim == 2
            n_subj = size(data,2);
            n_time = size(data,1);
            if params.verbose ==1
                fprintf('n subj: %d \nn timepoints: %d\n', n_subj, n_time)
            end

        else
            error('weird subject dimension. Not computing more than 2D matrix')
        end

%         %% Get rid of NAN
%         lgcl_nan = isnan(data);
%         if any(any(lgcl_nan))
%             check_consistency = sum(lgcl_nan,2);
%             exc_zero = check_consistency(check_consistency>0);
%             if all(exc_zero==n_subj)
%         %         this means that the amount of NaN is stable across participants.
%         %         Hence it is possible to cancel from data, with a warning
%                 indx_nan = find(lgcl_nan(:,1));
%                 data(indx_nan,:) = [];
%                 params.time_bins(indx_nan) = [];
%                 if params.verbose >= 0
% 
%                     warning('found consistent NaN across participants. Cutting them out')
%                     fprintf('\n\n##########################################################\n\n')
% 
%                 end
%             else
%                 error('found inconsistent NaN across particpants')
%             end 
%         end

        %% Check whether time bins match with the signal
        l_signal = size(data,1);
        l_signal_tBin = numel(params.time_bins);

        if l_signal~=l_signal_tBin
            error('number of time bins and the length of the sampled signal mismatch')
        end

        time_window = params.time_bins(end)-params.time_bins(1);
        curr_Fs = (numel(params.time_bins)-1)/time_window;

        if curr_Fs~=params.f_sample
            if params.verbose >= 0
                warning('detected discrepancy between declared and effective frequency sampling')
                fprintf('\ndeclared fSample: %d \ndetected fSample: %f2 \n', ...
                    params.f_sample, curr_Fs)
                fprintf('\n\n##########################################################\n')
            end
        end

        %% compute the required order of detrend on data

        if params.detrend_flag>=0

            data = local_do_detrend(data);
            outStruct.det_data = data;

            if params.verbose==1
                fprintf('applied %d order detrend', params.detrend_flag)
                fprintf('\n\n##########################################################\n\n')
            end

        else
            if params.verbose==1
                disp('no detrend applied')
            end

        end

        %% apply lowpass filter if specified

        if params.lp_filter 

            warning('lp filter still in beta version!')
            params.lp_filter_props = [];
            params.lp_filter_props.Fpass = params.f_sample/2 - 2;% % arbitrary indeed
            params.lp_filter_props.Fstop = params.f_sample/2-1;
            params.lp_filter_props.Apass = 1; % default values in matlab filter builder
            params.lp_filter_props.Astop = 60;

            data = local_do_LPfilter(data, params);

        end

        %% apply hanning window if specified

        if isempty(params.window)==0 

            data = local_do_window(data);

            if params.verbose==1
                disp([params.window ' window applied on data'])
            end

        else

            if params.verbose==1
                disp('no window applied')
            end

        end



%% spectra computation
padded_length = l_signal+params.zero_pad*2;
cmplx_signal = fft(data, padded_length)/padded_length;
out_fft = abs(cmplx_signal); 

if mod(padded_length,2)==0
    end_signal = padded_length/2+1;
else
    end_signal = ceil(padded_length/2);
end

% compute power if needed
if params.power
    
    out_fft = 2* out_fft.^2;
    
    if params.verbose==1
        disp('computed the power (as 2* out_fft.^2 )')
        disp('http://www.fieldtriptoolbox.org/tutorial/fourier')
    end
    
end

cut_cmplx = cmplx_signal(1:end_signal, :);
cut_amplitude = out_fft(1:end_signal, :);

%% output
outStruct.cmplx_out = cut_cmplx;
outStruct.spctr_out = cut_amplitude; % phase info removed
outStruct.cmplx_out_long = cmplx_signal;
outStruct.freqs = linspace(0,params.f_sample/2, end_signal);
outStruct.params = params;

%% HELPER FUNCTIONS


%% detrend
% several order

    function detrended_data = local_do_detrend(data)
        % helper to apply multiple order detrend on data
        % assume data structure as before: each column a subject, evolution
        % of time series along first dimension (rows)

        detrended_data = nan(size(data));

        for iPoly = 1:n_subj

            vect_data = data(:, iPoly);
            current_poly = polyfit(params.time_bins, vect_data, params.detrend_flag);
            trend = polyval(current_poly, params.time_bins);
            detrended_data(:, iPoly) = vect_data-trend;

        end

    end

%% LP filter: 
%not properly an anti alias, but ideally subserve the same function of
%smoothing signal too close to nyquist limit

    function filtered_data = local_do_LPfilter(data, params)
        
        filtered_data = nan(size(data));
        
        persistent filt_LP
        
            if isempty(filt_LP)
                
                h_filter = ...
                    fdesign.lowpass('fp,fst,ap,ast', params.lp_filter_props.Fpass,... 
                    params.lp_filter_props.Fstop, params.lp_filter_props.Apass,...
                    params.lp_filter_props.Astop, params.f_sample);
                
                filt_LP = design(h_filter, 'equiripple', ...
                            'MinOrder', 'any', ...
                            'StopbandShape', 'flat');
                
                set(filt_LP, 'PersistentMemory', true)
                
            end
            

        
        for iFilt = 1:n_subj
            
            filtered_data(:,iFilt) = filter(filt_LP, data(:,iFilt));
            
        end
        
    end


%% hanning windowing of data 

    function windowed_data = local_do_window(data)
        % helper to apply windowing on data
        % hanning, hamming, tukey allowed so far
        % assume data structure as before: each column a subject, evolution
        % of time series along first dimension (rows)
        
        switch params.window
            
            case 'hanning'
        
                this_window = hanning(l_signal);
                multi_wins = repmat(this_window, 1, n_subj);
                windowed_data = data.*multi_wins;
                
            case 'tukey'
                
                this_window = tukeywin(l_signal);
                multi_wins = repmat(this_window, 1, n_subj);
                windowed_data = data.*multi_wins;
                
            case 'hamming'
                
                this_window = hamming(l_signal);
                multi_wins = repmat(this_window, 1, n_subj);
                windowed_data = data.*multi_wins;
                
        end
                
                
    end
    

end