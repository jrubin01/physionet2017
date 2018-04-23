function  adjRpeak = qrs_corrector_wrap(full_path,annot,fs)
% load(full_path)
val = load_physionet_2017(full_path,fs);
signal = val;

if isempty(annot)
    rPeak = [];
else
    rPeak = annot;
end
%     rPeak = ann_gqrs(:,1);
%     plot(data);hold on;plot(rPeak,data(rPeak),'*r')

%% adjust the beat label time to the peak
adjRpeak = rPeak; % allocate memory
winLenMsec = 75;
winLenSamp = round(winLenMsec*fs/1000);
for rPeakIdx = 1:length(rPeak)
    thisPeakSampNo = rPeak(rPeakIdx)+1;
    lowWinIdx = max(1,thisPeakSampNo-winLenSamp);
    hiWinIdx = min(length(signal),thisPeakSampNo+winLenSamp);
    [maxValue, maxIdx] = max(signal(lowWinIdx:hiWinIdx));
    [minValue, minIdx] = min(signal(lowWinIdx:hiWinIdx));
    midPntVal = mean(signal(lowWinIdx:hiWinIdx));
    if (midPntVal-minValue) < (maxValue-midPntVal)
        adjRpeak(rPeakIdx) = (lowWinIdx + maxIdx - 2);
    else
        adjRpeak(rPeakIdx) = (lowWinIdx + minIdx - 2);
    end
end