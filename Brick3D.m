function [gcoord,nodes,faces]=Brick3D(xm,xp,ym,yp,zm,zp,Nx,Ny,Nz)

x_space = linspace(xm,xp,Nx);
y_space = linspace(ym,yp,Ny);
z_space = linspace(zm,zp,Nz);

[X,Y,Z]=meshgrid(x_space,y_space,z_space);

ind = 1:(Nx*Ny*Nz); ind = reshape(ind,[Ny Nx Nz]);

gcoord = [X(:) Y(:) Z(:)];

% elements will be constructed using the same node-ordering in x-y plane as
% with 2D elements, with the 5-6-7-8 nodes occuring on the next higher z
% level

ind1 = ind(1:(end-1),1:(end-1),1:(end-1));
ind1 = ind1(:);
ind2 = ind(1:(end-1),2:end,1:(end-1)); ind2=ind2(:);
ind3 = ind(2:end,2:end,1:(end-1)); ind3=ind3(:);
ind4 = ind(2:end,1:(end-1),1:(end-1)); ind4=ind4(:);
ind5 = ind(1:(end-1),1:(end-1),2:end); ind5=ind5(:);
ind6 = ind(1:(end-1),2:end,2:end); ind6=ind6(:);
ind7 = ind(2:end,2:end,2:end); ind7=ind7(:);
ind8 = ind(2:end,1:(end-1),2:end); ind8=ind8(:);

nodes=[ind1,ind2,ind3,ind4,ind5,ind6,ind7,ind8];

faces.zy_m = find(gcoord(:,1)==xm);
faces.zy_p = find(gcoord(:,1)==xp);
faces.zx_m = find(gcoord(:,2)==ym);
faces.zx_p = find(gcoord(:,2)==yp);
faces.xy_m = find(gcoord(:,3)==zm);
faces.xy_p = find(gcoord(:,3)==zp);