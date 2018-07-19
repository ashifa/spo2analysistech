

function [valid_segments, invalid_segments, validity_mask, signal_diff] = ppg_segmentation2(ppg, Fs)

% Look for corrupted segments of the input PPG signal.
%
% Inputs:     ppg signal in mV.
%             Sampling frequency of the input signal (Hz).
%
% Outputs:    Matrix of valid segments start/end indexes: [start1 end1]
%                                                         [start2 end2]
%                                                         [   ...     ]
%
%             Matrix of invalid segments start/end indexes: [start1 end1]
%                                                           [start2 end2]
%                                                           [   ...     ]
%             Signal validity mask:
%             Synchronized with the input PPG signal, with values:  1: valid signal
%                                                                   0: corrupted signal
%

ppg = ppg * 1e-4;  

% High-pass filtering
diff_order          = round(13*Fs/360);
% diff_order          = round(35*Fs/360);
diff_filt           = [1 zeros(1, diff_order) -1];
signal_diff         = filter(diff_filt, 1, ppg);


mask_window	= ones(1, round(1.5*Fs));
window_l	= length(mask_window);
% threshold	= -4;
threshold	= -8; % 07/2016 ref
% threshold	= -12;
% threshold	= 14; % 07/2016
% threshold	= -inf;

inx = (signal_diff < threshold); % ref
% inx = (abs(signal_diff) >= abs(threshold)); % 07/07/2016 
mask = conv(double(inx), mask_window); % as.numeric()
validity_mask(mask > 1) = 0;
validity_mask(mask < 1) = 1;
validity_mask = validity_mask(round(window_l/2): length(signal_diff));


aux             = [1 validity_mask];
invalid_start	= find(diff(aux) < 0);
invalid_end     = find(diff(aux) > 0);
invalid_end = invalid_end - 1;
if length(invalid_start) > length(invalid_end), invalid_end = [invalid_end length(validity_mask)]; end;


aux         = [0 validity_mask];
valid_start	= find(diff(aux) > 0);
valid_end   = find(diff(aux) < 0);
valid_end = valid_end - 1;
if length(valid_start) > length(valid_end), valid_end = [valid_end length(validity_mask)]; end;


invalid_segments	= [invalid_start' invalid_end'];
valid_segments      = [valid_start' valid_end'];







