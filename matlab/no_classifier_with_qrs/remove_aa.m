function [x_section,y_section,section_index,section_type] = remove_aa(x_section,y_section,section_index,section_type,filter_dim)
% Just select a special points of section in special part of Phase Space

x_lim1 = filter_dim(1,1);
x_lim2 = filter_dim(1,2);
y_lim1 = filter_dim(1,3);
y_lim2 = filter_dim(1,4);

%plot(orig,delayed,'.')
idx = find(( x_section<=x_lim1 |  x_section>=x_lim2) | (y_section<=y_lim1 | y_section>=y_lim2));
x_section(idx)=[];
y_section(idx)=[];
section_index(idx) =[];
section_type(idx)=[];
%figure;plot(orig,delayed,'.')
