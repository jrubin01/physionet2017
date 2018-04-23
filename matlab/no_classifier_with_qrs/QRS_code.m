ecg = val(1,:);
fs = 250;

b1 = BaseLine1(ecg,fs*.3,'md');
ecg_fil = ecg-b1;

wt = modwt(ecg,5);
wtrec = zeros(size(wt));
wtrec(4:5,:) = wt(4:5,:);
y = imodwt(wtrec,'sym4');
y = abs(y).^2;
[~,R_index] = findpeaks(y,'MinPeakHeight',0.15,...
    'MinPeakDistance',40);

[~,R_index] = findpeaks(ecg_fil,'MinPeakHeight',50,...
    'MinPeakDistance',40);

% [qrs_amp_raw,R_index,~] = pan_tompkin(ecg_fil,fs,0);

figure;
plot(ecg)
hold on
plot(R_index,ecg(R_index),'ro')

rtnd = difffit(R_index')/fs*1000; % IBI
rtnd(end)=[];
hr = 60000./rtnd;

respiration = val(4,:);
respiration = -respiration;
% b1 = BaseLine1(respiration,fs*.3,'md');
% respiration = respiration-b1;
[peaksr,locsr] = findpeaks(respiration,'MinPeakHeight',1,'MinPeakDistance',900);
figure;
plot(respiration)
hold on
plot(locsr,peaksr,'ro')

rrd = difffit(locsr')/fs*1000; % IBI
rrd(end)=[];
rr = 60000./rrd;


  