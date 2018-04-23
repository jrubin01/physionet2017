function [dist] = ps2distance(x_section,y_section,x_ref,y_ref)

%Distance
dist = sqrt((x_section-x_ref).^2 + (y_section-y_ref).^2);
ind = find(x_section<x_ref | y_section<y_ref);
dist(ind)= -dist(ind);

% if x>y
%     d = norm(x-y)
% else
%     d = -norm(x-y)
% end