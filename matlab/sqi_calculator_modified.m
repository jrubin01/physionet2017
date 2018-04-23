function SQI = sqi_calculator_modified(datapath)
%Calculates SQI for ECG or PPG.
%taken from CO's code. Adapted by PC 21/06/2013. Adapted by Saman 2/23/2017
%Inputs:    data - 1D time series of data (ecg or ppg)
%           beatsSample - indices of beats (qrs spikes or ppg pulses) in sample
%           fs - sampling rate
%           IN - threshold for average correlation coefficient (default 0.66 ecg, 0.86 ppg)
%           plot_ex - for plot mode it should be 1.
%Outputs:   R2 - average correlation coefficient
%           SQI - 0 if bad and 1 if good
fs = 300;
%IN =0.66;
IN = 0.5;
plot_ex = 0;

%Load required data
load(datapath)
data = val'./1000; %convert amplitude to mV
clear val;
[ann_pt2,sign,en_thres] = qrs_detect2(data,0.25,0.6,fs);
beatsSample = ann_pt2';
%Initialize variables
SQI = 0;
R2 = 0;
avtempl = 0;
ts = 0;

%Find mean RR interval (sample) to define size of template
hrs = floor(mean(diff(beatsSample)));
% beatsTime = beatsSample./fs;
% HR = 60./diff(beatsTime); %Heart rate, [bpm]
% %Rule1: Check 40<=HR<=180
% rule1 = find (HR<40 | HR>180);
% %Rule2: Check all RR intervals <= 3 Second
% rule2 = find (diff(beatsTime)>3);
% %Rule3: ratio of max to min RR interval should be less than 2.2
% rule3 = find(max(diff(beatsTime))/min(diff(beatsTime))>2.2);
% rulesCheck = double([~isempty(rule1) ~isempty(rule2) ~isempty(rule3)]);
% %% Apply rules: Rule 1 (40<HR) || Rule 1 (HR<180) || Rule 2 (are there any gaps > 3s?) || Rule 3 (ratio of max to min RR interval should be less than 2.2)
% if ~isempty(rule1) || ~isempty(rule2) || ~isempty(rule3)
%     R2 = NaN;
%     SQI = 0;
% else
    %% template matching-first identify QRS complexes which have a complete template within the window of data.
    ts = [];
    j = find(beatsSample>hrs/2);
    l = find(beatsSample+floor(hrs/2)<length(data));
    if isempty(l)==1
        %It means that there was not enough QRS to create template
        return
    else
        %find QRS complexes and create time series cycle
        for k = j(1):l(end)
            t = data(beatsSample(k)-floor(hrs/2):beatsSample(k)+floor(hrs/2));
            tt = t/norm(t); tt = tt(:)';
            ts = [ts;tt];
        end
    end
    
    %find all templates in current window
    avtempl = mean(ts,1);
    
    %Now calculate correlation for every beat in this window
    r2 = nan(size(ts,1),1);
    for k = 1:size(ts,1)
        r2(k) = corr2(avtempl,ts(k,:));
    end
    
    %Calculate mean correlation coefficient
    R2 = mean(abs(r2));
    
    if R2<IN
        SQI = 0;
    else
        SQI = 1;
    end
    
    %% Plot template and inidividual beats
    if plot_ex
        paper_size = [6, 5];
        figure('Position', [200, 200, 100*paper_size(1), 100*paper_size(2)], 'Color',[1 1 1])
        lwidth1 = 3; lwidth2 = 2; ftsize = 14;
        time = 0:(length(avtempl)-1); time = time./fs;
        hold on,
        for beat_no = 1 : size(ts,1)
            plot(time, ts(beat_no,:), 'color', 0.7*[1 1 1], 'LineWidth', lwidth2)
        end
        plot(time, avtempl, 'r', 'LineWidth', lwidth1)
        %     set(gca, 'YTick', [])
        xlabel('Time [s]', 'FontSize', ftsize)
        xlim([0, time(end)])
        if IN==0.86
            ylab=ylabel('PPG', 'FontSize', ftsize, 'Rotation', 0);
        else
            ylab=ylabel('ECG', 'FontSize', ftsize, 'Rotation', 0);
        end
        set(ylab, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
        %     set(gca, 'FontSize', ftsize, 'XTick', [])
        %set ylim
        rang = range(ts(:));
        ylim([min(ts(:))-0.1*rang, max(ts(:))+0.1*rang]);
        %     if SQI && IN==0.86
        %         save_name = 'ppg_high_qual';
        %         annotation('textbox',[0.01, 0.9, 0.1,0.1],'String','(b)', 'FontSize', ftsize, 'LineStyle', 'None')
        %     elseif SQI && IN==0.66
        %         save_name = 'ecg_high_qual';
        %         annotation('textbox',[0.01, 0.9, 0.1,0.1],'String','(a)', 'FontSize', ftsize, 'LineStyle', 'None')
        %     elseif ~SQI && IN==0.66
        %         save_name = 'ecg_low_qual';
        %         annotation('textbox',[0.01, 0.9, 0.1,0.1],'String','(c)', 'FontSize', ftsize, 'LineStyle', 'None')
        %     elseif ~SQI && IN==0.86
        %         save_name = 'ppg_low_qual';
        %         annotation('textbox',[0.01, 0.9, 0.1,0.1],'String','(d)', 'FontSize', ftsize, 'LineStyle', 'None')
        %     end
    end
% end
end