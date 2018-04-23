function annot = qrs_detector_wrap(full_path,fs,qrs_method,start_sample,stop_sample)
%Wrapper to load ECG data for Physionet Challenge 2017 and perform QRS detection.
% 1- Load the data 
% 2- Removes baseline wander
% 3- Correct lead connection (if it is needed)
% 4- Perform QRS detection and returns annotations and saves the annotations in text file
%% Load data
val = load_physionet_2017(full_path,fs);

data = val'; %data should be in form of vector n*1
% data = val'./1000; %convert amplitude to mV [data should be in form of vector n*1
clear val b1;
%% Make sure the input argument are valid
result_path = ['.' filesep 'matlab' filesep 'results' filesep'];
if nargin<2
    fs = 300;
    qrs_method = 'jqrs';
    start_sample = 1;
    stop_sample = length(data);
elseif nargin < 3
    qrs_method = 'jqrs';
    %Select the whole data if start_sample and stop_sample are not provided
    start_sample = 1;
    stop_sample = length(data);
elseif nargin<4
    %Select the whole data if start_sample and stop_sample are not provided
    start_sample = 1;
    stop_sample = length(data);
elseif narging<5
    stop_sample = length(data);
end
if start_sample<1
    start_sample = 1;
end
if stop_sample>length(data)
    stop_sample = length(data);
end
% More work is required to capture all exceptions.

%%
data = data(start_sample:stop_sample);
[N,~] = size(data);
%     temp = selected_folder{ii};
%     recordName = temp(end-9:end-4);
[pathstr,recordName,~] = fileparts(full_path);
clear full_path;
%% QRS Detection
if isequal(qrs_method,'jqrs')
    %=== call jqrs
    annot = run_qrsdet_by_seg_ali(data,fs);
    %Saving the results of QRS detection
    fid = fopen([recordName '.jqrs'], 'w');
    fprintf(fid, '%f\n', annot);
    fclose(fid);
    %Move data to the results folder
    movefile('*.jqrs',[result_path 'qrs' filesep]);
elseif isequal(qrs_method,'gqrs')
    %% Write the data out to wfdb
    %this is required because gqrs run on .dat data
    if ~exist([recordName '.dat'], 'file') || ~exist([recordName '.hea'], 'file')
        wrsamp((transpose(1:N)-1),data*1000,recordName,fs,1000,'16+24');  %Save the data to .dat and creates .hea
    end
    %=== call gqrs
    % gqrs_thresh = 1;
    gqrs(recordName,[],[],1,[],'gqrs');     %save the data to .gqrs file
    % load in gqrs
    try
        annot = rdann(recordName,'gqrs');
        annot  = annot(:,1);
    catch
        % rdann sometimes crashes when loading empty qrs annotations
        annot = [];
    end
    %Move data to the results folder
    movefile('*.gqrs',[result_path 'qrs' filesep]);
elseif isequal(qrs_method,'sqrs')
    %% Write the data out to wfdb
    %this is required because sqrs run on .dat data
    if ~exist([recordName '.dat'], 'file') || ~exist([recordName '.hea'], 'file')
        wrsamp((transpose(1:N)-1),data*1000,recordName,fs,1000,'16+24');
    end
    %=== call sqrs
    sqrs(recordName,[],[],1,[]); %save the data to .qrs file
    % load in sqrs
    try
        annot = rdann(recordName,'qrs');
        annot  = annot(:,1);
    catch
        % rdann sometimes crashes when loading empty qrs annotations
        annot = [];
    end
    %Move data to the results folder
    movefile('*.qrs',[result_path 'qrs' filesep]);
elseif isequal(qrs_method,'wqrs')
    %% Write the data out to wfdb
    %this is required because wqrs run on .dat data
    if ~exist([recordName '.dat'], 'file') || ~exist([recordName '.hea'], 'file')
        wrsamp((transpose(1:N)-1),data*1000,recordName,fs,1000,'16+24');
    end
    %=== call wqrs
    wqrs(recordName,[],[],1,[],[],[],[]); %save the data to .wqrs file
    % load in sqrs
    try
        annot = rdann(recordName,'wqrs');
        annot  = annot(:,1);
    catch
        % rdann sometimes crashes when loading empty qrs annotations
        annot = [];
    end
    %Move data to the results folder
    movefile('*.wqrs',[result_path 'qrs' filesep]);
elseif isequal(qrs_method,'pan1')
    %===  call Pan-Tompkins method
    annot = pantomQRSdetector(data,fs);
    %Saving the results of QRS detection
    fid = fopen([recordName '.pan1'], 'w');
    fprintf(fid, '%f\n', annot);
    fclose(fid);
    %Move data to the results folder
    movefile('*.pan1',[result_path 'qrs' filesep]);
elseif isequal(qrs_method,'pan2')
    %===  call Pan-Tompkins method (challenge code)
    [annot,sign,en_thres] = qrs_detect2(data,0.25,0.6,fs);
    annot = annot';
    %Saving the results of QRS detection
    fid = fopen([recordName '.pan2'], 'w');
    fprintf(fid, '%f\n', annot);
    fclose(fid);
    %Move data to the results folder
    movefile('*.pan2',[result_path 'qrs' filesep]);
elseif isequal(qrs_method,'wavelet')
    %=== Wavelet-based QRS
    % https://www.mathworks.com/examples/wavelet/mw/wavelet-ex77408607-r-wave-detection-in-the-ecg
    wt = modwt(data,5); %decompose the ECG waveform down to level 5 using the default 'sym4' wavelet.
    wtrec = zeros(size(wt));
    wtrec(4:5,:) = wt(4:5,:);   %reconstruct a frequency-localized version of the ECG waveform using only the wavelet coefficients at scales 4 and 5
    y = imodwt(wtrec,'sym4');
    y = abs(y).^2;
    [~,annot] = findpeaks(y,'MinPeakHeight',0.04,'MinPeakDistance',50);
    annot = annot';
    fid = fopen([recordName '.wavelet'], 'w');
    %Saving the results of QRS detection
    fprintf(fid, '%f\n', annot);
    fclose(fid);
    %Move data to the results folder
    movefile('*.wavelet',[result_path 'qrs' filesep]);
elseif isequal(qrs_method,'hilbert')
    %=== Hilbert QRS detection
    %biosig should be installed
    H2 = qrsdetect(data,fs,1);
    % Extract QRS-info according to BIOSIG/T200/EVENTCODES.TXT
    idx = find(H2.EVENT.TYP == hex2dec('0501'));
    annot = H2.EVENT.POS(idx);
    %Saving the results of QRS detection
    fid = fopen([recordName '.hilber'], 'w');
    fprintf(fid, '%f\n', annot);
    fclose(fid);
    %Move data to the results folder
    movefile('*.hilber',[result_path 'qrs' filesep]);
elseif isequal(qrs_method,'filterbank')
    %===  Filter Bank QRS Detection
    H2 = qrsdetect(data,fs,2);
    % Extract QRS-info according to BIOSIG/T200/EVENTCODES.TXT
    idx = find(H2.EVENT.TYP == hex2dec('0501'));
    annot = H2.EVENT.POS(idx);
    %Saving the results of QRS detection
    fid = fopen([recordName '.filterbank'], 'w');
    fprintf(fid, '%f\n', annot);
    fclose(fid);
    %Move data to the results folder
    movefile('*.filterbank',[result_path 'qrs' filesep]);
elseif isequal(qrs_method,'dpi')
    %===  Dynamic Plosion Index
    annot = dpi_qrs(data,fs,1800,5)';
    %Saving the results of QRS detection
    fid = fopen([recordName '.dpi'], 'w');
    fprintf(fid, '%f\n', annot);
    fclose(fid);
    %Move data to the results folder
    movefile('*.dpi',[result_path 'qrs' filesep]);
elseif isequal(qrs_method,'F2AMM3')
    [R_inds,Q_inds,S_inds,QRS_On_II,QRS_Off_II,ecg_hat,Peak_activities]=F2AMM3_old(data,fs,1);
    annot = R_inds;
end
% movefile('*.hea',[result_path,'data\']);
% movefile('*.dat',[result_path,'data\']);