
function processed_ecg_mask = discard_short_segments2(processed_ecg_mask, Fs, discard_time)


discrad_l = round(discard_time*Fs);

% aux = [0; processed_ecg_mask];
aux = [0 processed_ecg_mask];
d   = diff(aux);
dp  = find(d > 0);
dn  = find(d < 0);

for k = 1:length(dp)
    if ((dn(k) - 1 - dp(k)) <= discrad_l);
        processed_ecg_mask(dp(k):(dn(k) - 1)) = 0;
    end
end


  
