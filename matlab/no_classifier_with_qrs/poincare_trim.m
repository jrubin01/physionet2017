function [x_section_post,y_section_post,section_index_post,section_type_post,section_repeatation] = poincare_trim(x_section,y_section,section_index,section_type)
%Remove Same Points in Poincare Points
x_section_post = x_section;
y_section_post = y_section;
section_index_post = section_index;
section_type_post = section_type;
section_index = [section_index;ones(1,size(section_index,2))];
for i= 1: size(section_index,2)-1
    for j = i+1:size(section_index,2)
        if (x_section(1,i) == x_section(1,j)) & (y_section(1,i) == y_section(1,j))
            section_index(2,i) = section_index(2,i) + 1;
            x_section_post(1,j) = -10000;
            y_section_post(1,j) = -10000;
            section_index_post(1,j) = -10000;
            section_type_post(1,j) = -10000;
        end
    end
end
%Remove Same points
index = find(x_section_post == -10000);
x_section_post(index) = [];
index = find(y_section_post == -10000);
y_section_post(index) = [];
index = find(section_index_post == -10000);
section_index_post(index) = [];
index = find(section_type_post == -10000);
section_type_post(index) = [];
%Find Number of poincare section point repetition
section_repeatation = section_index(2,:);
section_repeatation(index) = [];
% %For Check
% %Remove Later
% sum(section_repeatation)


    
