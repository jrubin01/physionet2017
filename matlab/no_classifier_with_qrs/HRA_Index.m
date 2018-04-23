function [PI, GI, SI] = HRA_Index(RR)
% Written by: Sadaf Moharerri 
% Edited by: Saman Parvaneh, PhD
% The goal of this function is to calculate the well-know heart rate assymetry (HRA) parametrs
% Input: RR (RR interval)
% Outputs: PI (Porta's HRA Index), GI (Guzik's HRA Index), and SI (Slope Index)

RRn = RR(1:(length(RR)-1));  %RRn  is corresponding to RR(n+1)
RRnm = RR(2:length(RR));     %RRnm  is corresponding to RR(n+1)
%% Calculating Porta's HRA Index
% Porta, Alberto, et al. "Temporal asymmetries of short-term heart period variability are linked to 
% autonomic regulation." American Journal of Physiology-Regulatory, Integrative and Comparative 
% Physiology 295.2 (2008): R550-R557.
DRR = RRnm - RRn;
u = find(DRR > 0);
b = find(DRR < 0);
PI = (length(b)/(length(b)+length(u))) * 100; % Porta's HRA Index
%% Calculating Guzik's HRA Index
% Piskorski, J., and P. Guzik. "Geometry of the Poincaré plot of RR intervals and its asymmetry in healthy adults." Physiological measurement 28.3 (2007): 287.
D = (RRnm - RRn)./sqrt(2);
Du = D(D > 0);
Db = D(D < 0);
% Guzik's HRA Index
GI = (sum(abs(Du)) / (sum(abs(Du)) + sum(abs(Db))))*100;
%% Calculating Slope Index
% Karmakar, C. K., A. H. Khandoker, and M. Palaniswami. "Phase asymmetry of heart rate variability signal." Physiological measurement 36.2 (2015): 303.
tetta = atand(RRnm./RRn);
Rtetta = 45 - tetta;
Rtetta_u = Rtetta(Rtetta < 0);
Rtetta_b = Rtetta(Rtetta > 0);
% Slope Index
SI = (sum(abs(Rtetta_u))/(sum(abs(Rtetta_u)) + sum(abs(Rtetta_b)))) * 100;