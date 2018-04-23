function b = BaseLine1(x,L,approach)
%
% b = BaseLine1(x,L,approach),
% Baseline wander extraction from biomedical recordings, using a single 
% stage of median or moving average filtering.
%
% inputs:
% x: vector or matrix of noisy data (channels x samples)
% L: averaging window length (in samples)
% approach:
%   'md': median filtering
%   'mn': moving average
%
% output:
% b: vector or matrix of baseline wanders (channels x samples)

N = size(x,2);
b = zeros(size(x));
flen = floor(L/2);

if (strcmp(approach,'mn'))      % moving average filter
    for j = 1:N,
        index = max(j-flen,1):min(j+flen,N);
        b(:,j) = mean(x(:,index),2);
    end
elseif (strcmp(approach,'md'))  % median filter
    for j = 1:N,
        index = max(j-flen,1):min(j+flen,N);
        b(:,j) = median(x(:,index),2);
    end
end