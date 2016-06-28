function fIn_b = velocityBC_3D(fIn_b,w,ex,ey,ez,ux_p,uy_p,uz_p)

rho = sum(fIn_b,2);
ux = (fIn_b*ex')./rho;
uy = (fIn_b*ey')./rho;
uz = (fIn_b*ez')./rho;

dx = ux_p-ux;
dy = uy_p-uy;
dz = uz_p-uz;

numSpd = length(w);

for spd=2:numSpd
   cu=(3)*(ex(spd)*dx+ey(spd)*dy+ez(spd)*dz);
   fIn_b(:,spd)=fIn_b(:,spd)+w(spd)*(rho.*cu);
       
end

