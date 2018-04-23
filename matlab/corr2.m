function r = corr2(varargin)

[a,b] = ParseInputs(varargin{:});

a = a - mean2(a);
b = b - mean2(b);
r = sum(sum(a.*b))/sqrt(sum(sum(a.*a))*sum(sum(b.*b)));

%--------------------------------------------------------
function [A,B] = ParseInputs(varargin)

narginchk(2,2);

A = varargin{1};
B = varargin{2};

validateattributes(A, {'logical' 'numeric'}, {'real','2d'}, mfilename, 'A', 1);
validateattributes(B, {'logical' 'numeric'}, {'real','2d'}, mfilename, 'B', 2);

if any(size(A)~=size(B))
    error(message('images:corr2:notSameSize'))
end

if (~isa(A,'double'))
    A = double(A);
end

if (~isa(B,'double'))
    B = double(B);
end










