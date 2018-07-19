


function [peak_r_inx, rr_int_v, peak_value_v] = cardio_peak_detect2(signal, type, Fs, start_valid_data, varargin)
%
% Cardiac signal peak detection
%
% Input:
%
% signal:               Input signal.
% type:                 - 1 = ECG
%                       - 2 = PPG
%                       - 3 = PPG, min peaks
% Fs:                	Sampling frequency (Hz).
% start_valid_data:		Peaks found before this time will be discarded (secs).
%                       Typical value 1.5 secs.
%
% Output:
%
% peak_r_inx:           Vector of peaks indexes in number of samples. One sample corresponds to 1/Fs secs.
%                       Time is synchronized with the input signal.
% rr_int_v:             Vector of intervals between peaks in number of samples. One sample corresponds to 1/Fs secs.
% peak_value_v:         Input signal sampled at the vector peaks.
%




% clear;
% close all;


if type == 2
    type = 3;
    maxPeaks = true;
else
    maxPeaks = false;
end



switch type
    case 1 % ECG

        diff_order      = round(13*Fs/360);

    case 2 % PPG
        diff_order      = round(45*Fs/360);
    case 3 % PPG
        diff_order      = round(90*Fs/360); % ref
        signal          = -signal;
end



samples_before_valid	= round(start_valid_data*Fs);	% Samples before valid data.


% Double differentiation (band-pass filtering enhancing the feature to extracted) and squaring
diff_filt           = [1 zeros(1, diff_order) -1];
signal_diff         = filter(diff_filt, 1, signal);
signal_diff2        = filter(diff_filt, 1, signal_diff);
signal_diff2_sqr	= signal_diff2.^2;
% -------------------------------------------------------------------


% Peak detection parameters
threshold_ratio         = 0.27; % Detection threshold adaption to peak amplitude.
rr_int_max_ratio_init   = 1.1;  % If "current rr"/"average rr" > rr_int_max_ratio_init => restart peak search
                                % from the previous r-peak with a lower detection threshold.

% rr_int_max_ratio_inc	= 1.04; % Original
rr_int_max_ratio_inc	= 1.1;
threshold_dec           = 0.775;    % Peak detection threshold decrease factor per searchback iteration.
threshold_min_ratio     = 0.0541;	% The detection threshold doesn't go below: avg("peaks amplitude")*threshold_min_ratio.
rr_int_min_ratio        = 0.36226;	% Used to update the minimum interval where two consecutive peaks shouldn't be valid.
% rr_int_min              = round(100*Fs/360); % Initial value; updated during the peak search.
% Debug: start with longer min inti beat interval
rr_int_min              = round(150*Fs/360); % Initial value; updated during the peak search..
avg_win1                = 5;        % threshold_min average window.
avg_win2                = 20;       % rr_int_avg average window.
nb_iter                 = 20;       % Peak searchback maximum number of iterations.
% rr_int_inf              = round(90*Fs/360);
% rr_int_sup              = round(120*Fs/360);
if type == 2

%     rr_int_inf              = round(70*Fs/360);
%     rr_int_sup              = round(100*Fs/360);
    rr_int_inf              = round(80*Fs/360);
    rr_int_sup              = round(100*Fs/360);
elseif type == 3
    rr_int_inf              = round(70*Fs/360); % ref
    rr_int_sup              = round(100*Fs/360); % ref
else
%     rr_int_inf              = round(110*Fs/360); % type 3
%     rr_int_sup              = round(130*Fs/360);
    rr_int_inf              = round(100*Fs/360);
    rr_int_sup              = round(120*Fs/360);
end

% --------------------------------------------------------------------------------------------------


% Internal & output variables/vectors
threshold_min           = 0; % Updated during peak detection.
rr_int_max_ratio        = rr_int_max_ratio_init;
rr_int_cnt              = 0;
threshold               = threshold_min;
last_peak_amplitude     = 0;
peak_found              = false;
iter_cnt                = 0;
rr_int_avg              = -1;
peak_r_inx              = [];
peak_r_mag              = [];
rr_int_v                = [];
rr_int_avg_v            = [];
% --------------------------------------------------------------------------------------------------


threshold_v = []; % debug

if ~isempty(varargin)
    posPeaks = 1;
else
    posPeaks = 0;
end


%% R peaks detection -------------------------------------------------------------------------------
% k = diff_order + 2 + round(8*Fs/125);
k = diff_order + 2 + round(7*Fs/125);
while k <= length(signal_diff2_sqr)
    

%     if k >= 526500
%         bidon=0
%     end
    
    threshold_v(k) = threshold; % debug
    
    switch type
        case 1
            if posPeaks
                slope_test = (sign(signal_diff(k - 1)) < 0) && (sign(signal_diff(k - 1 - diff_order - 1)) > 0);
            else
                slope_test = (sign(signal_diff(k - 1)) * sign(signal_diff(k - 1 - diff_order -1))) < 0;
            end

        case 2
            % Positive and negative slopes respectively to the left and right of the peak
%             slope_test = (sign(signal_diff(k - 1) - round(1*Fs/125)) <= 0) && (sign(signal_diff(k - 1 - diff_order - round(8*Fs/125))) >= 0); 

%             slope_test = (sign(signal_diff(k - 1) - round(1*Fs/125)) <= 0) && (sign(signal_diff(k - 1 - diff_order - round(7*Fs/125))) >= 0);% ref
            
            slope_test = (sign(signal_diff(k - 1)) < 0) && (sign(signal_diff(k - 1 - diff_order - 1)) > 0);
        case 3
%             slope_test = (sign(signal_diff(k - 1) - round(Fs / 125)) <= 0) && (sign(signal_diff(k - 1 - diff_order + round(4*Fs/125))) > 0);
            slope_test = (sign(signal_diff(k - 1)) < 0) && (sign(signal_diff(k - 1 - diff_order - 1)) > 0);
    end
    
    
    if      (signal_diff2_sqr(k - 1) > threshold) && (signal_diff2_sqr(k) <= signal_diff2_sqr(k - 1)) &&...
            (signal_diff2_sqr(k - 1) > signal_diff2_sqr(k - 2)) && slope_test
            % Peak shape found.
            % Discard some artifacts by checking the signal derivative before and after the peak.
        
%           if (rr_int_cnt > 0) && (rr_int_cnt < round(120*Fs/350))
%           if (rr_int_cnt > 0) && (rr_int_cnt < round(115*Fs/350))  
%         if (rr_int_cnt > 0) && (rr_int_cnt < round(110*Fs/350)) 
%         if (rr_int_cnt > 0) && (rr_int_cnt < round(105*Fs/350))
%         if (rr_int_cnt > 0) && (rr_int_cnt < round(100*Fs/350))

%         if (rr_int_cnt > 0) && (rr_int_cnt < round(90*Fs/350))
            
        if (rr_int_cnt > 0) && (rr_int_cnt < rr_int_min)
            % After the 1st peak.
            % No close consecutive peaks.
            
            if signal_diff2_sqr(k - 1) > last_peak_amplitude
                % Replace the last peak by the current one
                if ~isempty(rr_int_v) % No rr_interval fo the 1st peak.
                    rr_int_v(end) = rr_int_v(end) + rr_int_cnt;
                end
                peak_r_inx(end)         = k - 1;
                peak_r_mag(end)         = signal_diff2_sqr(k - 1);
                last_peak_amplitude     = signal_diff2_sqr(k - 1);
                threshold_min           = mean(peak_r_mag(max(length(peak_r_mag) - avg_win1, 1):...
                                                length(peak_r_mag)))*threshold_min_ratio;
                                        % Proportional to a peaks amplitude sliding window average.
                
                threshold               = max(threshold_min, signal_diff2_sqr(k - 1)*threshold_ratio);
                rr_int_cnt              = 1;
                rr_int_max_ratio        = rr_int_max_ratio_init; % Restore the initial value
                iter_cnt                = 0;
                peak_found              = true;
                
            end
            
        else
            % Isolated peak (not close to the previous one).
            
            peak_r_inx              = [peak_r_inx (k - 1)];
            if (length(peak_r_inx) > 1) && (peak_r_inx(end - 1) > samples_before_valid)
                % At least 2 peaks.
                % Measure starts after "samples_before_valid".
                rr_int_v	= [rr_int_v rr_int_cnt];
            end
            peak_r_mag              = [peak_r_mag signal_diff2_sqr(k - 1)];
            last_peak_amplitude     = signal_diff2_sqr(k - 1);
            threshold_min           = mean(peak_r_mag(max(length(peak_r_mag) - avg_win1, 1):length(peak_r_mag)))*threshold_min_ratio;
            threshold               = max(threshold_min, signal_diff2_sqr(k - 1)*threshold_ratio);
            rr_int_cnt              = 1; % Peak found with 1 sample delay.
            rr_int_max_ratio        = rr_int_max_ratio_init; % Restore the initial value.
            iter_cnt                = 0;
            peak_found              = true;
            
        end
        
    end
    
    
    % Between R peaks
    
    % Recover from eventual initial large artifact detection
%     if (length(peak_r_inx) <= 3) && (rr_int_cnt/Fs > 1.0)
%     if (length(peak_r_inx) <= 45) && (rr_int_cnt/Fs > 1.0)
    if (rr_int_cnt/Fs > 1.0)
%     if (length(peak_r_inx) <= 20) && (rr_int_cnt/Fs > 1.0)
        threshold	= threshold*threshold_dec;
        peak_r_mag	= 0;
    end
    
    
    if rr_int_cnt ~= 0                  % Only if the 1st peak has already been found.
        rr_int_cnt = rr_int_cnt + 1;	% Measure RR interval.
    end
    
    %% RR average interval update
    if (~isempty(rr_int_v) && (k > samples_before_valid) && peak_found && (rr_int_cnt > rr_int_min))
        % Update after "samples_before_valid", when a new peak has been found.
        
        peak_found = false;
        rr_int_avg = mean(rr_int_v(max(length(rr_int_v) - (avg_win2 - 1), 1):end));
        
        % Debug: clip with min value
%         rr_int_avg = max(rr_int_avg, round(60/110*Fs));
        
        rr_int_avg_v        = [rr_int_avg_v rr_int_avg];      
        rr_int_min          = max(rr_int_inf, round(rr_int_min_ratio*rr_int_avg));
        rr_int_min          = min(rr_int_min, rr_int_sup);
        
    end
    %% --------------------------
    
    
    
    %% Potential missed peak: no peaks were detected within 'rr_int_max' samples
    if (rr_int_avg > 0) && (rr_int_cnt > round(rr_int_avg*rr_int_max_ratio)) &&...
       (iter_cnt < nb_iter)
        % New search with a lower threshold in an extended window
        iter_cnt            = iter_cnt + 1;
        rr_int_max_ratio	= rr_int_max_ratio*rr_int_max_ratio_inc;
        
        % Go back right after the last valid found peak and restart the search with a lower threshold
        k               = peak_r_inx(end) + 1 + rr_int_min;
        threshold       = max(threshold*threshold_dec, threshold_min);
        rr_int_cnt      = 1 + rr_int_min; % Peak found with 1 sample delay
    end
    %% ---------------------------------------------------------------------------------------------
    
    k = k + 1;
    
end % Peak search
%% -------------------------------------------------------------------------------------------------



peak_r_inx(peak_r_inx < samples_before_valid) = []; % Remove invalid peaks.


% delay           = diff_order + 1;
% delay           = diff_order + 2; % type 2
if (type==1) || (type == 2)
    delay           = diff_order + 1;
elseif type == 3
    delay           = diff_order - 3;
%     delay           = diff_order + 2;
end

peak_r_inx      = peak_r_inx - delay; % Compensate the filters delay.


% Local peak search ----------
% diff_f           = [1 zeros(1, 5) -1];
% tmp         = filter(diff_f, 1, signal);
if (length(peak_r_inx) > 1) && (type ~= 1)
    avg_window = 10;
    d_peak_r_inx = diff(peak_r_inx);
    d_peak_r_inx = [d_peak_r_inx(1) d_peak_r_inx];
    local_avg_beat_int = zeros(1, length(peak_r_inx));
    for k = 1:length(peak_r_inx)
        
        local_avg_beat_int(k) = mean( d_peak_r_inx(max(k - round(avg_window/2), 1):...
                                min(k + round(avg_window/2) - 1, length(d_peak_r_inx))) );
     
        if type == 2
            shift = round(1/4 * local_avg_beat_int(k)); % ref
%             shift = round(1/2 * local_avg_beat_int(k));
        elseif type == 3
%             shift = round(1/8 * local_avg_beat_int(k));
            shift = round(1/4 * local_avg_beat_int(k));
        end

        if ( ((peak_r_inx(k) - shift) > 0) && ( (peak_r_inx(k) + shift) <=  length(signal) ) )
            int = ( peak_r_inx(k) - shift: peak_r_inx(k) + shift );
            inx = find(signal(int) == max(signal(int)), 1);
%             inx = find(tmp(int) == max(tmp(int)), 1);
            peak_r_inx(k) = peak_r_inx(k) + inx - shift - 1;
        end
        
    end
end
% ----------------------------

peak_r_inx2 = [];
if maxPeaks
    for k = 1:length(peak_r_inx) - 1
        aux = signal(peak_r_inx(k):peak_r_inx(k + 1));
%         peak_r_inx2(k) = peak_r_inx(k) + find(aux == min(aux), 1); % bug xxxxxxxxxxxxx
        peak_r_inx2(k) = peak_r_inx(k) + find(aux == min(aux), 1) - 1;
    end
    peak_r_inx = peak_r_inx2;
end



rr_int_v        = diff(peak_r_inx);
peak_value_v	= signal(peak_r_inx);




