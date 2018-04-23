function y_pred = classifierNO(full_path)
fs = 300;
sqi_thresh = 0.66;
plot_mode = 0;
duration_second = 9;
%y_pred = classifierNO('A00002.mat')
%% Load data and do qrs detection
[val,~,~] = load_physionet_2017(full_path,fs,duration_second);
% [qrs_fusion,sqi_qrs_confirmation,val] = qrs_fusion(full_path,fs,plot_mode);
%gqrs
annot_gqrs = qrs_detector_wrap(full_path,fs,'gqrs');
annot_gqrs = qrs_corrector_wrap(full_path,annot_gqrs,fs);
annot_gqrs(annot_gqrs==0)=[];
qrs_method = annot_gqrs;
%Removes small and negative peaks
if ~isempty(qrs_method)
    if sum(val(qrs_method))>=0
        %             idx = find(val(gqrs)<=0.02);
        %             gqrs(idx) = [];
    else
        %             idx = find(val(gqrs)>0.02);
        %             gqrs(idx) = [];
        val = -val;
    end
end
%Pan-tompkins
annot_pan = qrs_detector_wrap(full_path,300,'pan2');
annot_pan(annot_pan==0)=[];
%% Signal quality measures
[R2,sqi_correlation,avtempl,ts] = sqi_calculator(val,qrs_method,fs,sqi_thresh,plot_mode);

%     [R2(2),~,~,~] = sqi_calculator(val,annot_gqrs,fs,sqi_thresh,plot_ex);
%     [R2(3),~,~,~] = sqi_calculator(val,annot_pan,fs,sqi_thresh,plot_ex);

[N_J,N_G] = bsqi(annot_pan/fs, annot_gqrs/fs,0.15);
bSQI = sum(N_J)/(size(annot_gqrs,1)+size(annot_pan,1)-sum(N_J));
features_sqi = [R2 bSQI];
%% RR features
if length(qrs_method)>2
    rr =  diff(qrs_method)/fs;
    hr = 60./rr;
    [PI, GI, SI] = HRA_Index(rr);
    features_rr = [length(qrs_method) min(rr) max(rr) mean(rr) median(rr) std(rr) sqrt(sum(rr.^2)/length(rr)) mean(hr) PI GI SI];
else
    rr= NaN;
    hr = NaN;
    features_rr = nan(1,11);
end
%% Spatial filing index features
if length(avtempl)>2
    data = (val-min(val))/(max(val)-min(val));
    %     data = (avtempl-min(avtempl))/(max(avtempl)-min(avtempl));
    delay = 4;
    blockNumber = 20;
    [orig,delayed] = series2ps_2d(data,delay,'no-normalized',plot_mode);
    [C,XEDGES,YEDGES] = histcounts2(orig,delayed,[blockNumber blockNumber],'XBinLimits',[0,1],'YBinLimits',[0,1],'Normalization','probability');
    Q = C.^2;
    %     QQ{1,ii = C;
    spatial = sum(sum(Q,2))/((blockNumber)^2);
    features_spatial = [C(:)' spatial];
else
    features_spatial = nan(1,401);
end
%% Poincare Features
if length(avtempl)>2
    sig = val/(max(abs(val)));
    %     sig = avtempl/(max(abs(avtempl)));
    delay =4;
    [orig,delayed] = series2ps_2d(sig,delay,'normalized',plot_mode);
    intersection_method = 2;
    m = 1;   % Slope of Poincare line
    b = 0;   % Intercept of Poincare line
    [pre_xsection,pre_ysection,pre_section_index,pre_section_type] = poincarecut_2d(orig,delayed,m,b,intersection_method);
    %--------------------------------------------------------------------------
    AA_filtering = 'no'; %To reduce effect of QRS complex & T wave
    if isequal(AA_filtering,'yes')
        % Region for AA selection
        % [xlim1 xlim2 ylim1 ylim2]
        filter_dim = [-0.05 0.05 -0.05 0.05];
        [prefilt_xsection,prefilt_ysection,prefilt_section_index,prefilt_section_type] = remove_aa(pre_xsection,pre_ysection,pre_section_index,pre_section_type,filter_dim);
    else
        filter_dim = [-1 1 -1 1];
        [prefilt_xsection,prefilt_ysection,prefilt_section_index,prefilt_section_type] = remove_aa(pre_xsection,pre_ysection,pre_section_index,pre_section_type,filter_dim);
    end
    
    [post_xsection,post_ysection,post_section_index,post_section_type,section_repeatation] = poincare_trim(prefilt_xsection,prefilt_ysection,prefilt_section_index,prefilt_section_type); %Remove same Points in Poincare Points
    % Step 6: Feature Extraction from Points of Section and its Figures
    x_ref = 0;  %Reference Points
    y_ref = 0;
    dist = ps2distance(post_xsection,post_ysection,x_ref,y_ref);
    diff_poincare_point = ps2diffdistance(pre_xsection,pre_ysection,x_ref,y_ref);
    % Plot Related Figures
    if isequal(plot_mode,'plot')
        poincare2timedomain   %Function to display points of section in time-domain
    end
    %We need "section probability" for plot and feature extraction parts
    section_probab = section_repeatation / sum(section_repeatation);
    % if sum(section_probab)>1 | sum(section_probab)<0
    %     disp('Error in section probability calculation')
    % end
    
    %Plot Figures Related to Poincare Section if it is selected
    if isequal(plot_mode,'plot')
        poincare_plot(orig,delayed,m,b,dist,diff_poincare_point,section_probab,post_xsection,post_ysection)
    end
    
    n_dist = 39; %12
    n_diffdist = 49;
    if isequal(plot_mode,'plot')
        figure('Name','Relative Frequency Histogram (for Points of Section)','NumberTitle','off')
        [N, X] = hist(dist,n_dist);
        bar(X, N./sum(N), 1);
        ylabel('Relative Frequency');
        xlabel('Points of Section Distances')
        axis tight
        %     title('Histogram of points of section')
        figure('Name','Relative Frequency Histogram (for Difference Between Points','NumberTitle','off')
        [N, X] = hist(diff_poincare_point,n_diffdist);
        bar(X, N./sum(N), 1);
        ylabel('Relative Frequency');
        xlabel('Distances between Consecutive Points of Section')
        axis tight
        %     title('Histogram of difference between points')
    end
    % Feature Extraction
    feature_numeric = [size(prefilt_section_index,2) size(post_section_index,2)];
    feature_dist = section_dist_feature(dist,n_dist);
    %             feature_diffdist = section_diff_feature(diff_poincare_point,n_diffdist);
    %             feature_multisection = section_multisection_feature(section_probab);
    %             features_poincare = [feature_numeric feature_dist feature_diffdist feature_multisection];
    features_poincare = [feature_numeric feature_dist];
else
    features_poincare = nan(1,13);
    %         features_poincare = nan(1,28);
    
end
%% Time domain features
%     %Method 1
%     if ~isempty(avtempl)
%         if skewness(avtempl)<0
%             avtempl = -avtempl;
%         end
%
%         [r_amp,r_sample,~,~] = findpeaks(avtempl,'MinPeakHeight',0.1,'npeaks',5,'SortStr','ascend');
%         r_amp = r_amp(end);
%         r_sample = r_sample(end);
%         DataInv = 1.01*max(avtempl) - avtempl;
%         [Minima,MinIdx] = findpeaks(DataInv);
%         Minima = avtempl(MinIdx);
%
%         idx_q = find(MinIdx<r_sample);
%         q_sample = MinIdx(idx_q(end));
%         q_amp = Minima(idx_q(end));
%         idx_t = find(MinIdx>r_sample);
%         s_sample = MinIdx(idx_t(1));
%         s_amp = Minima(idx_t(1));
%
%         if plot_ex
%             plot(avtempl)
%             hold on
%             plot(r_sample,r_amp,'*r')
%             plot(q_sample,q_amp,'*r')
%             plot(t_sample,t_amp,'*r')
%         end
%
% %         p_template = avtempl(1:q_sample);
% %         pwave_features = wave_extraction(p_template);
%
% %         qrs_template = avtempl(q_sample:s_sample);
% %         qrs_features = wave_extraction(qrs_template);
%         qrs_features = [];
%         qrs_features = [qrs_features (r_amp-q_amp) (r_amp-s_amp) (r_amp-q_amp)/(r_amp-s_amp) (s_sample-q_sample)];
%
% %         t_template = avtempl(s_sample:end);
% %         twave_features = wave_extraction(t_template);
%
%         p_series = [];
%         for i = 1:size(ts,1)
%             p_series = [p_series, ts(i,1:q_sample)];
%         end
%         %plot(p_series)
%         saen = SampEn(2,0.2*std(p_series), p_series,4);
%     else
%         saen = NaN;
%     end
%     features_time = [qrs_features saen];
%     features_time = [pwave_features qrs_features twave_features saen];
features_time = [];
%% Frequecy features
features_freq = get_features_frequency2(val,fs);

features = [features_sqi features_freq features_rr features_time features_spatial features_poincare];

%%
load('learned_parms_NO_gqrs')
features = (features-mu1)./sigma1;
y_pred  = boosted_predict(features,cvres_full.clf_full,'probability');
