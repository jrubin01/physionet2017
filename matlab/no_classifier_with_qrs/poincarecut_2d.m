function [x_section,y_section,section_index,section_type] = poincare_2d(orig,delayed,m,b,poincare_type)
%Find Poincare Section on 2D Phase Space
number_of_datapoints = size(orig,2);
i = 1:number_of_datapoints;
check(1,i) = delayed(1,i) - m.*orig(1,i) - b;
% Type 1 Poincare Section: Exact Poincare Section
if (poincare_type == 1)
    exact_section = find (check == 0);
    x_section = orig(1,exact_section);
    y_section = delayed(1,exact_section);
    section_index = exact_section;
    section_type = 10*ones(1,size(section_index,2)); %10=Exact 20=Interpolated
    %--------------------------------------------------------------------------
    %Type 2 Poincare Section: Single Sided FromLeft to Right (+-)
elseif (poincare_type == 2)
    section_pn = [];
    section_type = [];
    j = 1;
    for  i = 1:number_of_datapoints-1
        if (check(1,i)==0 & check(1,i+1)==0)       %Exact Section
            section_pn = [section_pn i];
            section_type = [section_type 10];
            x_section(1,j) = orig(1,i);
            y_section(1,j) = delayed(1,i);
            j = j + 1;
        elseif (check(1,i)==0 & check(1,i+1)<0)    %Exact Section
            section_pn = [section_pn i];
            section_type = [section_type 10];
            x_section(1,j) = orig(1,i);
            y_section(1,j) = delayed(1,i);
            j = j + 1;
        elseif (check(1,i)>0 & check(1,i+1)<=0)  %+- Section with Interpolation
            section_pn = [section_pn i];
            section_type = [section_type 20];
            if orig(1,i) == orig(1,i+1)
                x_section(1,j) = orig(1,i);
                y_section(1,j) = m .* x_section(1,j) + b;
            else
                assume_m = (delayed(1,i)-delayed(1,i+1))/(orig(1,i)-orig(1,i+1));
                assume_b = delayed(1,i)- assume_m.* orig(1,i);
                x_section(1,j) = (b-assume_b)/(assume_m-m);
                y_section(1,j) = m .* x_section(1,j) + b;
            end
            j = j + 1;
        end
    end
    section_index = section_pn;
    %--------------------------------------------------------------------------
    %Type 3 Poincare Section: Single Sided From Right to Left (-+)
elseif (poincare_type == 3)
    section_np = [];
    section_type = [];
    j = 1;
    for  i = 1:number_of_datapoints-1
        if (check(1,i)==0 & check(1,i+1)==0)       %Exact Section
            section_np = [section_np i];
            section_type = [section_type 10];
            x_section(1,j) = orig(1,i);
            y_section(1,j) = delayed(1,i);
            j = j + 1;
        elseif (check(1,i)==0 & check(1,i+1)>0)    %Exact Section
            section_np = [section_np i];
            section_type = [section_type 10];
            x_section(1,j) = orig(1,i);
            y_section(1,j) = delayed(1,i);
            j = j + 1;
        elseif (check(1,i)<0 & check(1,i+1)>=0)    %-+ Section with Interpolation
            section_np = [section_np i];
            section_type = [section_type 20];
            if orig(1,i) == orig(1,i+1)
                x_section(1,j) = orig(1,i);
                y_section(1,j) = m .* x_section(1,j) + b;
            else
                assume_m = (delayed(1,i)-delayed(1,i+1))/(orig(1,i)-orig(1,i+1));
                assume_b = delayed(1,i)- assume_m.* orig(1,i);
                
                x_section(1,j) = (b-assume_b)/(assume_m-m);
                y_section(1,j) = m .* x_section(1,j) + b;
            end
            j = j + 1;
        end
    end
    section_index = section_np;
    %--------------------------------------------------------------------------
    %Type 4 Poincare Section: Double Sided
elseif (poincare_type == 4)
    section_double = [];
    section_type = [];
    j = 1;
    for  i = 1:number_of_datapoints-1
        if (check(1,i)==0 & check(1,i+1)==0)        %Exact Section
            section_double = [section_double i];
            section_type = [section_type 10];
            x_section(1,j) = orig(1,i);
            y_section(1,j) = delayed(1,i);
            j = j + 1;
        elseif (check(1,i)==0 & check(1,i+1)<0)     %Exact Section
            section_double = [section_double i];
            section_type = [section_type 10];
            x_section(1,j) = orig(1,i);
            y_section(1,j) = delayed(1,i);
            j = j + 1;
        elseif (check(1,i)>0 & check(1,i+1)<=0)     %+- Section with Interpolation
            section_double = [section_double i];
            section_type = [section_type 20];
            if orig(1,i) == orig(1,i+1)
                x_section(1,j) = orig(1,i);
                y_section(1,j) = m .* x_section(1,j) + b;
            else
                assume_m = (delayed(1,i)-delayed(1,i+1))/(orig(1,i)-orig(1,i+1));
                assume_b = delayed(1,i)- assume_m.* orig(1,i);
                x_section(1,j) = (b-assume_b)/(assume_m-m);
                y_section(1,j) = m .* x_section(1,j) + b;
            end
            j = j + 1;
        elseif (check(1,i)==0 & check(1,i+1)>0)    %Exact Section
            section_double = [section_double i];
            section_type = [section_type 10];
            x_section(1,j) = orig(1,i);
            y_section(1,j) = delayed(1,i);
            j = j + 1;
        elseif (check(1,i)<0 & check(1,i+1)>=0) %-+ Section with Interpolation
            section_double = [section_double i];
            section_type = [section_type 20];
            if orig(1,i) == orig(1,i+1)
                x_section(1,j) = orig(1,i);
                y_section(1,j) = m .* x_section(1,j) + b;
            else
                assume_m = (delayed(1,i)-delayed(1,i+1))/(orig(1,i)-orig(1,i+1));
                assume_b = delayed(1,i)- assume_m.* orig(1,i);
                
                x_section(1,j) = (b-assume_b)/(assume_m-m);
                y_section(1,j) = m .* x_section(1,j) + b;
            end
            j = j + 1;
        end
    end
    section_index = section_double;
end

%%For Test
% clc
% clear all
% close all
% sig = [-5  -5 1 4 4 4 3 2 0 1 2 1 1 1 5 5];
% delay = 1;
% orig = sig(1,1:size(sig,2)-delay)/max(abs(sig));
% delayed = sig(1,delay+1:size(sig,2))/max(abs(sig));
% orig = [-0.2 0.1 0.4 0.4 0.8 0.9 1]
% delayed = [-0.2 0.1 0 0.6 0.6 0.9 1]
% % orig = [0.2 0.2];
% % delayed= [0.4 0];
% %
% % delayed = [0.2 0.2];
% % orig = [0 0.4];
% poincare_type = 1;
% %Plot Normalized Phase Space
% scatter(orig,delayed,'.')
% %plot(orig,delayed)
% axis tight
% grid on
% ylabel('x(t+\tau)')
% xlabel('x(t)')
% title('Phase Space of Signal')
% hold on