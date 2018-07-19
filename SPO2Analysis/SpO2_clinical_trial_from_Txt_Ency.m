
% DynoSense Corp. CONFIDENTIAL 
% Clincial for SPO2
% V2: corrected time-index of RR,SpO2 data.

% THiS PATH SETTING and FILE NAME
path_name = '../Data/';
file_name1 = 'cleanedFile.txt';
%file_name1 = 'received-2-clean.txt';

% MUST PREPROCESS OUTPUT DATA using MakeFileClean Program

M = csvread([path_name file_name1],0,0);
sf=size(M);
PPG1=M(1:sf(1),4);
PPG0=M(1:sf(1),3); 
time=M(1:sf(1),1);
PpgIrDcOffset=M(1:sf(1),(sf(2)-1));
PpgRedDcOffset=M(1:sf(1),(sf(2)-1)); 

ls = 0;
PPG0(1:ls)           = [];
PPG1(1:ls)           = [];
PpgIrDcOffset(1:ls)  = [];
PpgRedDcOffset(1:ls) = [];

PpgIrDcOffset(1)  = 1;
PpgRedDcOffset(1) = 1;
% -------------------------------

inxPpgIrDcOffset  = find(PpgIrDcOffset);
inxPpgRedDcOffset = find(PpgRedDcOffset);

PPG0 = PPG0(~isnan(PPG0));
PPG1 = PPG1(~isnan(PPG1));

%Debug Plotting Use this to Figure out if Start_value_data Ok?
figure
plot(inxPpgIrDcOffset - 1 ,PpgIrDcOffset(inxPpgIrDcOffset)), hold on
plot((0:length(PPG0) - 1), PPG0), hold on
plot(inxPpgRedDcOffset - 1 ,PpgRedDcOffset(inxPpgRedDcOffset))
plot((0:length(PPG1) - 1), PPG1)
figure
plot(time,PPG0);

DcIrMat = [inxPpgIrDcOffset(1) PpgIrDcOffset(inxPpgIrDcOffset(1))];
for k = 2:length(inxPpgIrDcOffset)
    if (PpgIrDcOffset(inxPpgIrDcOffset(k)) ~= PpgIrDcOffset(inxPpgIrDcOffset(k - 1)))
        DcIrMat = [DcIrMat; [inxPpgIrDcOffset(k) PpgIrDcOffset(inxPpgIrDcOffset(k))]];
    end
end
DcRedMat = [inxPpgRedDcOffset(1) PpgRedDcOffset(inxPpgRedDcOffset(1))];
for k = 2:length(inxPpgRedDcOffset)
    if (PpgRedDcOffset(inxPpgRedDcOffset(k)) ~= PpgRedDcOffset(inxPpgRedDcOffset(k - 1)))
        DcRedMat = [DcRedMat; [inxPpgRedDcOffset(k) PpgIrDcOffset(inxPpgRedDcOffset(k))]];
    end
end

if size(DcIrMat, 1) == 1
    DcIrMat = [DcIrMat; [length(PPG0) DcIrMat(1, 2)]];
end
if size(DcRedMat, 1) == 1
    DcRedMat = [DcRedMat; [length(PPG1) DcRedMat(1, 2)]];
end

PPG0_with_dc = zeros(1, length(PPG0));
for k = 1:(size(DcIrMat, 1) - 1)
    PPG0_with_dc(DcIrMat(k, 1):DcIrMat(k + 1, 1) - 1) = PPG0(DcIrMat(k, 1):DcIrMat(k + 1, 1) - 1) + DcIrMat(k, 2);
end
PPG0_with_dc(DcIrMat(k + 1, 1):end) = PPG0(DcIrMat(k + 1, 1):end) + DcIrMat(k + 1, 2);

PPG1_with_dc = zeros(1, length(PPG1));
for k = 1:(size(DcRedMat, 1) - 1)
    PPG1_with_dc(DcRedMat(k, 1):DcRedMat(k + 1, 1) - 1) = PPG1(DcRedMat(k, 1):DcRedMat(k + 1, 1) - 1) + DcRedMat(k, 2);
end
PPG1_with_dc(DcRedMat(k + 1, 1):end) = PPG1(DcRedMat(k + 1, 1):end) + DcRedMat(k + 1, 2);

%plot((0:length(PPG0) - 1),PPG0_with_dc)
%hold on
%plot((0:length(PPG1) - 1),PPG1_with_dc)
%figure 

Fs = 250;
Fs = Fs*(1 - 0.0174);
ppg2=PPG1_with_dc;
ppg1=PPG0_with_dc; 

%figure;
%plot (ppg1); hold on;
%plot (ppg2);

inx1 = find(isnan(ppg1));
inx2 = find(isnan(ppg2));
inx1 = union(inx1, inx2);
ppg1(inx1) = [];
ppg2(inx1) = [];
clear inx1 inx2

[valid_segments1, invalid_segments1, validity_mask1, signal_diff1] = ppg_segmentation2(-ppg1', Fs);
[valid_segments2, invalid_segments2, validity_mask2, signal_diff2] = ppg_segmentation2(-ppg2', Fs);

lppg = min(length(validity_mask1), length(validity_mask2));
validity_mask = validity_mask1(1:lppg).*validity_mask2(1:lppg); % Intersection of ppg0 and ppg1 segments

validity_mask = [0 validity_mask 0];
validity_mask = discard_short_segments2(validity_mask, Fs, 3.5);
validity_mask = validity_mask(2:end - 1);

aux         = [0 validity_mask];
valid_start	= find(diff(aux) > 0);
valid_end   = find(diff(aux) < 0);
valid_end = valid_end - 1;
if length(valid_start) > length(valid_end), valid_end = [valid_end length(validity_mask)]; end;

valid_segments      = [valid_start' valid_end'];

load lpf_butterworth_4_fc5Hz_fs250Hz;
ppg1 = filtfilt(b, a, ppg1);
ppg2 = filtfilt(b, a, ppg2);

%figure;
%plot (ppg1); hold on;
%plot (ppg2);

home
snum = 0;
valid_segments = [1 valid_segments(end,end)]; % Override segmentation.

for k = 1:size(valid_segments, 1)
    

    
    lsegment = (valid_segments(k, 2) - valid_segments(k, 1))/Fs;
    
    ppg_segment1         = ppg1(valid_segments(k, 1):valid_segments(k, 2)); % Extract segment
    ppg_segment2         = ppg2(valid_segments(k, 1):valid_segments(k, 2)); % Extract segment
    
    % THIS PARATMETER MIHGT REQUIRE change 
    % It is function of data
    start_valid_data	= 70.0; % (secs)
    
    ppg_segment = ppg_segment1; % Override max energy selection
    
    [peak_r_inx1, rr_int_v1, peak_value_v1] = cardio_peak_detect2(ppg_segment, 2, Fs, start_valid_data);
    [peak_r_inx2, rr_int_v2, peak_value_v2] = cardio_peak_detect2(ppg_segment, 3, Fs, start_valid_data);
    
    peak_r_inx1_c = peak_r_inx1;
    peak_r_inx2_c = peak_r_inx2;
    
    cnth = 0;
    cntl = 0;
    RR = [];
    hbr_inst = [];
    inx_p = [];
    cor = zeros(1, length(peak_r_inx1_c) - 1);
    cor_inx = zeros(1, length(peak_r_inx1_c) - 1);
    for j = 1:length(peak_r_inx1_c) - 1
        inx = find(peak_r_inx2_c > peak_r_inx1_c(j), 1);
        if isempty(inx)
            break
        end
        if peak_r_inx2_c(inx) < peak_r_inx1_c(j + 1)
            inx_p(j) = peak_r_inx1_c(j);
            RR(j) = log( ppg_segment2(peak_r_inx1_c(j))/ppg_segment2(peak_r_inx2_c(inx)) )/...
                    log( ppg_segment1(peak_r_inx1_c(j))/ppg_segment1(peak_r_inx2_c(inx)) );
                
            hbr_inst(j) = Fs/(peak_r_inx1_c(j + 1) - peak_r_inx1_c(j))*60;
                           
                aux1 = ppg_segment2(peak_r_inx1_c(j): peak_r_inx1_c(j + 1));
                aux2 = ppg_segment1(peak_r_inx1_c(j): peak_r_inx1_c(j + 1));
                aux = corrcoef(aux1, aux2);
                cor(j) = aux(2);
                cor_inx(j) = (peak_r_inx1_c(j) - 1)/Fs;
                if cor(j) < 0.99
                    cntl = cntl + 1;
                else
                    cnth = cnth + 1;
                end
        end 
    end
    
    inx = find(inx_p == 0);
    inx_p(inx) = inx_p(max(inx - 1, 1)) + .1;

    inx = (cor_inx == 0);
    cor_inx(inx) = [];
    cor2 = cor;
    cor2(inx) = [];
    figure
    plot(cor_inx, cor2, '-*c'), grid on
    
    %Must Factor Optimized for Each Device
    p = [-8.0074   12.8082  -34.4576  113.8471]; 
    RR = real(RR);
    SpO2_inst = p(1)*RR.^3 + p(2)*RR.^2 + p(3)*RR + p(4);
    
    SpO2_avg = mean(SpO2_inst);




  shft=0;
    t = (0:length(ppg_segment1) - 1)/Fs + shft;
    %figure
    %h1 = plot(t, ppg_segment1); hold on, grid on
    %plot((peak_r_inx1 - 1)/Fs + shft,ppg_segment1(peak_r_inx1),'ko')
    %plot((peak_r_inx2 - 1)/Fs + shft,ppg_segment1(peak_r_inx2),'go')
    %h2 = plot(t, ppg_segment2,'r');
    %plot((peak_r_inx1 - 1)/Fs + shft,ppg_segment2(peak_r_inx1),'ko')
    %plot((peak_r_inx2 - 1)/Fs + shft,ppg_segment2(peak_r_inx2),'go')
    %legend([h1 h2] , 'PPG0', 'PPG1')
    %xlabel('time (sec)')
 
    %figure
    %hold on
   SpO2_inst_med = medfilt1(SpO2_inst, 60);
   SpO2_inst_med(SpO2_inst_med > 100) = 100;
     plot((inx_p - 1 + 0*valid_segments(k, 1))/Fs + shft, round(SpO2_inst_med), 'm.-')
    title('Instantaneous SpO2')
    xlabel('time (sec)')
    aSpo2 = round(SpO2_inst_med);
    aSpo2=aSpo2';
    bInx= (inx_p - 1 + 0*valid_segments(k, 1))/Fs + shft;
    bInx=bInx';
    inx_pR=inx_p';
    Spo2d=[inx_pR,bInx, aSpo2];
    csvwrite ('SPo2Data.csv', Spo2d);
    
    % Process RR
     RR_med = medfilt1(RR, 11);% 7    
    figure
    HBR_med = medfilt1(hbr_inst, 20);
    HBR_avg = filter(ones(1, 20)/20, 1, hbr_inst);
    plot((inx_p - 1 + 0*valid_segments(k, 1))/Fs + shft, round(HBR_avg), 'g-*')
    title('Instantaneous HR (med filt)')
    xlabel('time (sec)')
    aHbr = round(HBR_avg);
    aHbr=aHbr';
    bInx= (inx_p - 1 + 0*valid_segments(k, 1))/Fs + shft;
    bInx=bInx';
    HRd=[inx_pR,bInx, aHbr];
    csvwrite ('HRData.csv',HRd);
    
    
end

