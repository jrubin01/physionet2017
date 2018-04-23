function [orig,delayed] = series2ps_2d(signal,delay,type,plt)
%Convert Time Series to 2d Phase Space
if nargin<2
    delay = 1;
    type = 'non-normalized';
    plt = 'noplot';
elseif nargin<3
    type = 'non-normalized';
    plt = 'noplot';
elseif nargin<4
    plt = 'noplot';
end

if isequal(type,'normalized')
    orig = signal(1,1:size(signal,2)-delay)/max(abs(signal));
    delayed = signal(1,delay+1:size(signal,2))/max(abs(signal));
else
    orig = signal(1,1:size(signal,2)-delay);
    delayed = signal(1,delay+1:size(signal,2));
end

if isequal(plt,'plot')
    %Plot the allignment of original and delayed signal
    %figure;subplot(2,1,1);plot(orig);axis('tight');subplot(2,1,2);plot(delayed);axis('tight')
    %Plot Regular Phase Space
    figure('Name','Reconstructed Phase Space (RPS)','NumberTitle','off')
    %scatter(orig,delayed)
    plot(orig,delayed,'.r',orig,delayed)
    axis tight
    ylabel('x(t+\tau)')
    xlabel('x(t)')
    title(['Phase Space of Signal (',type,')'])
    hold on
end