function [diff_poincare_point] = ps2diffdistance(x_section,y_section,x_ref,y_ref)

%Distance
dist = sqrt((x_section-x_ref).^2 + (y_section-y_ref).^2);
ind = find(x_section<x_ref | y_section<y_ref);
dist(ind)= -dist(ind);

diff_poincare_point = abs(diff(dist));
% diff_poincare_point = diff(dist);