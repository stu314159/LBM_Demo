function [w ex,ey,ez,bb_spd]=D3Q19_lattice_parameters()

% ex = [0 -1 0 0 -1 -1 -1 -1 0 0 1 0 0 1 1 1 1 0 0];
% ey = [0 0 -1 0 -1 1 0 0 -1 -1 0 1 0 1 -1 0 0 1 1];
% ez = [0 0 0 -1 0 0 -1 1 -1 1 0 0 1 0 0 1 -1 1 -1];

ex = [0 1 -1 0 0 0 0 1 -1 1 -1 1 -1 1 -1 0 0 0 0];
ey = [0 0 0 1 -1 0 0 1 1 -1 -1 0 0 0 0 1 -1 1 -1];
ez = [0 0 0 0 0 1 -1 0 0 0 0 1 1 -1 -1 1 1 -1 -1];


w = zeros(1,19);
for i = 1:19
    ind=abs(ex(i))+abs(ey(i))+abs(ez(i));
    switch ind
        case 0
            w(i)=1/3;
        case 1
            w(i) = 1/18;
            
        case 2
            w(i) = 1/36;
    end
end



all_spds = [ex' ey' ez'];
bb_spd = zeros(1,19);

for spd=1:19
    bb_spd(spd)= find((all_spds(:,1)==(-ex(spd))) & ...
        (all_spds(:,2)==(-ey(spd))) & (all_spds(:,3)==(-ez(spd))));
    
end
