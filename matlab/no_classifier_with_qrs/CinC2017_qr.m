clc
clear
close all
% addContainingDirAndSubDir
%%
full_path = cd;
selected_folder = uipickfiles('FilterSpec',full_path,'Type', { '*.mat',   'MAT-files'},'prompt','Select a file for loading');
if isequal(selected_folder,0) || isempty(selected_folder)
    msgbox('You didn''t select any folder for conversion (Run a program again)','Error','error')
    return
end

for ii= 1: length(selected_folder)
    full_path = selected_folder{ii};
    [~,recordName,~] = fileparts(full_path);
    fs = 300;
    %     load(selected_folder{ii})
    load(full_path)
    data = val'./1000; %convert amplitude to mV
    clear val;
    
        ann_pt2 = qrs_detector_wrap(full_path,fs,'pan2');

    %     fs = 300;
    %     %%
    %     [N,M] = size(data);
    %     temp = selected_folder{ii};
    %     recordName = temp(end-9:end-4);
    %     clear temp full_path;
    %% Write the data out to wfdb
    %TODO: this is required because gqrs/ ons run on matlab data directly
    %     wrsamp((transpose(1:N)-1),data*10000,recordName,fs,10000,'16+24');
        %% QRS Detection
    %     %=== call jqrs
    %     ann_jqrs = run_qrsdet_by_seg_ali(data,fs);
    
    %     %=== call gqrs
    %     % gqrs_thresh = 1;
    %     gqrs(recordName,[],[],1,[],'gqrs');
    %     % load in gqrs
    %     try
    %         ann_gqrs = rdann(recordName,'gqrs');
    %     catch
    %         % rdann sometimes crashes when loading empty qrs annotations
    %         ann_gqrs = [];
    %     end
    %
    %     %=== delete gqrs' output
    %     if SAVE_STUFF==0
    %         delete([recordName '.gqrs']);
    %     end
    %
    %     %=== call sqrs
    %     sqrs(recordName,[],[],1,[]);
    %     % load in sqrs
    %     try
    %         ann_sqrs = rdann(recordName,'qrs');
    %     catch
    %         % rdann sometimes crashes when loading empty qrs annotations
    %         ann_sqrs = [];
    %     end
    %
    %     %=== delete sqrs' output
    %     if SAVE_STUFF==0
    %         delete([recordName '.qrs']);
    %     end
    %
    %     %=== call wqrs
    %     wqrs(recordName,[],[],1,[],[],[],[]);
    %     % load in sqrs
    %     try
    %         ann_wqrs = rdann(recordName,'wqrs');
    %     catch
    %         % rdann sometimes crashes when loading empty qrs annotations
    %         ann_wqrs = [];
    %     end
    %
    %     %=== delete sqrs' output
    %     if SAVE_STUFF==0
    %         delete([recordName '.wqrs']);
    %     end
    %
    %     %===  call Pan-Tompkins method
    %     ann_pt = pantomQRSdetector(data,fs);
    %     %===  call Pan-Tompkins method (challenge code)
    %     [ann_pt2,sign,en_thres] = qrs_detect2(data,0.25,0.6,fs);
    %     ann_pt2 = ann_pt2';
    % %
    %     %=== Wavelet-based QRS
    %     % https://www.mathworks.com/examples/wavelet/mw/wavelet-ex77408607-r-wave-detection-in-the-ecg
    %     wt = modwt(data,5); %decompose the ECG waveform down to level 5 using the default 'sym4' wavelet.
    %     wtrec = zeros(size(wt));
    %     wtrec(4:5,:) = wt(4:5,:);   %reconstruct a frequency-localized version of the ECG waveform using only the wavelet coefficients at scales 4 and 5
    %     y = imodwt(wtrec,'sym4');
    %     y = abs(y).^2;
    %     [~,ann_wavelet] = findpeaks(y,'MinPeakHeight',0.15,'MinPeakDistance',50);
    %     ann_wavelet = ann_wavelet';
    %
    %     %=== Hilbert QRS detection
    %     %biosig should be installed
    %     H2 = qrsdetect(data,fs,1);
    %     % Extract QRS-info according to BIOSIG/T200/EVENTCODES.TXT
    %     idx = find(H2.EVENT.TYP == hex2dec('0501'));
    %     ann_hilbert = H2.EVENT.POS(idx);
    %
    %     %===  Filter Bank QRS Detection
    %     H2 = qrsdetect(data,fs,2);
    %     % Extract QRS-info according to BIOSIG/T200/EVENTCODES.TXT
    %     idx = find(H2.EVENT.TYP == hex2dec('0501'));
    %     ann_filterbank = H2.EVENT.POS(idx);
    %
    %     %===  Dynamic Plosion Index
    %     ann_dpi = dpi_qrs(data,fs,1800,5)';
    %
    %     clear en_thres H@ idx M N sign wt wtrec y
       %     %Saving the results of QRS detection
    %     fid = fopen([recordName '.dpi'], 'w');
    %     fprintf(fid, '%f\n', ann_dpi);
    %     fclose(fid);
    %     fid = fopen([recordName '.filterbank'], 'w');
    %     fprintf(fid, '%f\n', ann_filterbank);
    %     fclose(fid);
    %     fid = fopen([recordName '.hilber'], 'w');
    %     fprintf(fid, '%f\n', ann_hilbert);
    %     fclose(fid);
    %     fid = fopen([recordName '.jqrs'], 'w');
    %     fprintf(fid, '%f\n', ann_jqrs);
    %     fclose(fid);
    %     fid = fopen([recordName '.wavelet'], 'w');
    %     fprintf(fid, '%f\n', ann_wavelet);
    %     fclose(fid);
    [R2,sqi_correlation,rulesCheck] = sqi_calculator2(data',ann_pt2',fs,0.66,0);
    
    results.record{ii,1} = recordName;
    results.rsqi{ii,1} = R2;
    results.sqi{ii,1} = sqi_correlation;
    results.rulesCheck{ii,1} = rulesCheck;
end
load('D:\OneDrive - Philips\My Files\MATLAB\cinc2017\Jonathan\split0\REFERENCE.mat')
%%Evaluation
[C,order] = confusionmat(G_numeric,cell2mat(results.sqi))