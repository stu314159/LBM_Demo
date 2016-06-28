function stv = genStreamTgtVec3D(Nx,Ny,Nz,ex,ey,ez)

ind = 1:(Nx*Ny*Nz); ind = reshape(ind,[Ny Nx Nz]);

numSpd=length(ex);
nnodes = Nx*Ny*Nz;

stv = zeros(nnodes,numSpd);

for spd = 1:numSpd
   t = circshift(ind,[-ey(spd) -ex(spd) -ez(spd)]); t = t(:);
   stv(:,spd)=t;
        
end