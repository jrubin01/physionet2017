function features_freq = get_features_frequency2(data,sr)
NFFT = 256;
f = (0:NFFT/2-1)/(NFFT/2)*(sr/2);
freq_range = [1,15;15,30;30,45;45,60;60,75;75,90;90,150;5,14;5,50];
s1 = data;
% s1 = s1.*hamming(length(s1));
Ft = fft(s1,NFFT);
power = abs(Ft(1:NFFT/2));

for bin=1:size(freq_range,1)
    idx = (f>=freq_range(bin,1)) & (f<freq_range(bin,2));
    P_S(1,bin) = median(power(idx));
end
features_freq = [P_S P_S(end-1)/P_S(end)];
end
