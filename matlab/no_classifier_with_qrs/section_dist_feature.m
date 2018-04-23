%Features from Coordinate of Poincare Section
function [dist_features] = section_dist_feature(dist,n_dist)
mean_distance = mean(dist);
var_distance = var(dist);
cv_distance = abs(mean_distance/var_distance);
range_distance = range(dist);
interquartile_distance = iqr(dist);
percentile_distance = prctile(dist,[25 50 75]);
kurt_distance = kurtosis(dist);
skew_distance = skewness(dist);
h = hist(dist,n_dist);
%h = histogram_sam(dist,[min(dist),max(dist),n_dist]);
ent_distance = entropy_pdf(h);
%ent_distance = entropy_sam(dist);

%------------------------------
dist_features(1,1) = mean_distance;
dist_features(1,2) = var_distance;
dist_features(1,3)= cv_distance;
dist_features(1,4) = range_distance;
dist_features(1,5) = interquartile_distance;
dist_features(1,6:8) = percentile_distance;
dist_features(1,9) = kurt_distance;
dist_features(1,10) = skew_distance;
dist_features(1,11) = ent_distance;