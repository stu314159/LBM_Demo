function [w,ex,ey,ez,bb_spd]=D3Q15_lattice_parameters()

w = ones(1,15);
w(1)=2/9;
w(2:7)=1/9;
w(8:end)=1/72;

% ex = [0 1 -1 0 0 0 0 1 -1 -1 1 1 -1 -1 1];
% ey = [0 0 0 1 -1 0 0 1 1 -1 -1 1 1 -1 -1];
% ez = [0 0 0 0 0 1 -1 1 1 1 1 -1 -1 -1 -1];

ex = [0 1 -1 0 0 0 0 1 -1 1 -1 1 -1 1 -1];
ey = [0 0 0 1 -1 0 0 1 1 -1 -1 1 1 -1 -1];
ez = [0 0 0 0 0 1 -1 1 1 1 1 -1 -1 -1 -1];

%bb_spd = [1 3 2 4 5 7 6 14 15 12 13 10 11 8  9];
%         1 2 3 4 5 6 7 8  9  10 11 12 13 14 15

all_spds = [ex' ey' ez'];
bb_spd = zeros(1,15);

for spd=1:15
  bb_spd(spd)= find((all_spds(:,1)==(-ex(spd))) & ...
      (all_spds(:,2)==(-ey(spd))) & (all_spds(:,3)==(-ez(spd)))); 
    
end

