function y = mean2(x)

y = sum(x(:),'double') / numel(x);
