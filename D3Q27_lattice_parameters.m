function [w ex,ey,ez,bb_spd]=D3Q27_lattice_parameters()

ex = [0 -1 0 0 -1 -1 -1 -1 0 0 -1 -1 -1 -1 1 0 0 1 1 1 1 0 0 1 1 1 1];
ey = [0 0 -1 0 -1 1 0 0 -1 -1 -1 -1 1 1 0 1 0 1 -1 0 0 1 1 1 1 -1 -1];
ez = [0 0 0 -1 0 0 -1 1 -1 1 -1 1 -1 1 0 0 1 0 0 1 -1 1 -1 1 -1 1 -1];

w = zeros(1,27);
for i = 1:27
    ind=abs(ex(i))+abs(ey(i))+abs(ez(i));
    switch ind
        case 0
            w(i)=8/27;
        case 1
            w(i) = 2/27;
            
        case 2
            w(i) = 1/54;
            
        case 3
            w(i) = 1/216;
    end
end


all_spds = [ex' ey' ez'];
bb_spd = zeros(1,19);

for spd=1:27
    bb_spd(spd)= find((all_spds(:,1)==(-ex(spd))) & ...
        (all_spds(:,2)==(-ey(spd))) & (all_spds(:,3)==(-ez(spd))));
    
end


