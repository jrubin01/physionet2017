function [val,original_length,selected_length] = load_physionet_2017(full_path,fs,duration_second)
%% Load data
load(full_path)
% Selecting part of data
original_length = length(val);
% val = val(1,1:duration_second*fs); %Choose 9 seconds

%Removing ECG baseline
b1 = BaseLine1(val,fs*.3,'md');
val = val-b1;

%Correcting lead inversion if it is needed
if skewness(val)<-0.06
    val = -val;
end

%Normalizing ECG amplitude
% val = val/max(abs(val));
% val = (val-min(val))/(max(val)-min(val));
val = val/1000;
selected_length = length(val);
