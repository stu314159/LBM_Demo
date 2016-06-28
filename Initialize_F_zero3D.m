function [fIn,fOut,rho,ux,uy]=Initialize_F_zero3D(gcoord,ex,ey,ez,w,rho_init)

[~,numSpd]=size(ex);
[nnodes,~]=size(gcoord);

fIn = zeros(nnodes,numSpd);
rho = rho_init.*ones(nnodes,1);

ux = zeros(nnodes,1);
uy = zeros(nnodes,1);
uz = zeros(nnodes,1);

for i = 1:numSpd
   cu = 3*(ex(i)*ux + ey(i)*uy + ez(i)*uz);
   fIn(:,i) = w(i)*rho.*(1+cu+(1/2)*(cu.*cu)-(3/2)*(ux.^2 + uy.^2 + uz.^2));
    
end

fOut = fIn; % just to start 