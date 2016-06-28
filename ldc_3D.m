%ldc_3D.m
clear
clc
close('all');

initialization = 0;
% 0 = initialize fIn to zero speed
% 1 = initialize fIn to Poiseuille profile <--not used for this problem

dynamics = 1;
% 1 = LBGK
% 2 = TRT 
% 3 = MRT <--- D3Q15 and D3Q19 only

lattice_selection = 1;
% 1 = D3Q15
% 2 = D3Q19
% 3 = D3Q27

fluid = 1;
% 1 = glycerin
% 2 = glycol
% 3 = water

Num_ts = 10000;
ts_rep_freq = 100;
plot_freq = Num_ts-1;
Re = 5;

Lx_p = 1;
Ly_p = 1;
Lz_p = 1;

switch fluid
    case 1
        rho_p = 1260;
        nu_p = 1.49/rho_p;
        
    case 2
        rho_p = 965.3;
        nu_p = 0.06/rho_p;
        
    case 3
        rho_p = 1000;
        nu_p = 1e-3/rho_p;
        
end

Lo = Ly_p;
Uavg = nu_p*Re/Lo;

Uo = Uavg;
To = Lo/Uo;

Ld = 1; Td = 1; Ud = (To/Lo)*Uavg;
nu_d = 1/Re;

% convert to LBM units
dt = 5e-4;
Ny_divs = 32;
dx = 1/(Ny_divs-1);
u_lbm = (dt/dx)*Ud;
nu_lbm=(dt/(dx^2))*nu_d;
omega = get_BGK_Omega(nu_lbm);

u_conv_fact = (dt/dx)*(To/Lo);
t_conv_fact = (dt*To);
l_conv_fact = (dx*Lo);
f_conv_fact = (l_conv_fact^2)/(u_conv_fact^2);

rho_lbm = rho_p;

% generate LBM lattice
xm = 0; xp = Lx_p;
ym = 0; yp = Ly_p;
zm = 0; zp = Lz_p;

Ny = ceil((Ny_divs-1)*(Ly_p/Lo))+1;
Nx = ceil((Ny_divs-1)*(Lx_p/Lo))+1;
Nz = ceil((Ny_divs-1)*(Lz_p/Lo))+1;

[gcoord,~,faces]=Brick3D(xm,xp,ym,yp,zm,zp,Nx,Ny,Nz);
[nnodes,~]=size(gcoord);

ind = 1:(Nx*Ny*Nz); ind = reshape(ind,[Ny Nx Nz]);

x_space = linspace(xm,xp,Nx);
y_space = linspace(ym,yp,Ny);
z_space = linspace(zm,zp,Nz);
[X,Y,Z]=meshgrid(x_space,y_space,z_space);

quiv_sf = 3;
[x_quiv,y_quiv,z_quiv]=meshgrid(xm:(Lx_p/quiv_sf):xp,...
    ym:(Ly_p/quiv_sf):yp,zm:(Lz_p/quiv_sf):zp);

fprintf('dx = %g, dy = %g, dz = %g \n',x_space(2)-x_space(1),...
    y_space(2)-y_space(1),z_space(2)-z_space(1));


switch lattice_selection
    
    case 1
        [w,ex,ey,ez,bb_spd]=D3Q15_lattice_parameters();
        lattice = 'D3Q15';
    case 2
        [w,ex,ey,ez,bb_spd]=D3Q19_lattice_parameters();
        lattice = 'D3Q19';
    case 3
        [w,ex,ey,ez,bb_spd]=D3Q27_lattice_parameters();
        lattice = 'D3Q27';
        
end


streamTgtMat = genStreamTgtVec3D(Nx,Ny,Nz,ex,ey,ez);

numSpd = length(w);

LatticeSize = [Ny Nx Nz];
LatticeSpeeds = [ex; ey; ez];

snl = [faces.zx_m; faces.zx_p; faces.xy_m; faces.xy_p; faces.zy_p];
snl = unique(snl);

% lid-node-list is the zy plane where x is at a minimum
% may seem like an odd choice, but it obviously doesn't matter
lnl = faces.zy_m; 
lnl = setxor(lnl,intersect(lnl,snl)); % eliminate solid nodes from inl

ux_p_lid = zeros(length(lnl),1);
uy_p_lid = u_lbm*ones(length(lnl),1);
uz_p_lid = ux_p_lid;



% tag some nodes for visualization

z_pln = z_space(ceil(Nz/2));
vis_nodes = find(gcoord(:,3)==z_pln);

switch initialization
    case 0
        [fIn,fOut,rho,ux,uy]=Initialize_F_zero3D(gcoord,ex,ey,ez,w,rho_lbm);
    case 1
        % leave this out for now...I don't use it anyway...
    
end


switch dynamics
    
    case 1% BGK
        fEq = zeros(nnodes,numSpd);
    case 2 % TRT
        fEq = zeros(nnodes,numSpd);
        fNEq = zeros(nnodes,numSpd);
        fEven = zeros(nnodes,numSpd);
        fOdd = zeros(nnodes,numSpd);
        
    case 3 % MRT
        fEq = zeros(nnodes,numSpd);
        M = getMomentMatrix(lattice);
        S = getEwMatrixMRT(lattice,omega);
        omega_op = M\(S*M);
        
        
end

fprintf('Number of Lattice-points = %d.\n',nnodes);
fprintf('Number of time-steps = %d. \n',Num_ts);
%fprintf('Predicted execution time = %g.\n', predicted_ex_time);

fprintf('LBM viscosity = %g. \n',nu_lbm);
fprintf('LBM relaxation parameter (omega) = %g. \n',omega);
fprintf('LBM flow Mach number = %g. \n',u_lbm);

input_string = sprintf('Do you wish to continue? [Y/n] \n');

run_dec = input(input_string,'s');

if ((run_dec ~= 'n') && (run_dec ~= 'N'))
    
    fprintf('Ok! Cross your fingers!! \n');
    
    % commence time stepping
   % profile on
    tic;
    for ts = 1:Num_ts
        
        % say something comforting about my progress...
        if(mod(ts,ts_rep_freq)==0)
            fprintf('Executing time step number %d.\n',ts);
        end
        
        % compute density
        rho = sum(fIn,2);
        
        % compute velocities
        ux = (fIn*ex')./rho;
        uy = (fIn*ey')./rho;
        uz = (fIn*ez')./rho;
        
        % set macroscopic and Microscopic Dirichlet-type boundary
        % conditions
        
        % macroscopic BCs
        ux(lnl)=ux_p_lid;
        uy(lnl)=uy_p_lid;
        uz(lnl)=uz_p_lid;
        
        % microscopic BCs
        fIn(lnl,:)=velocityBC_3D(fIn(lnl,:),w,ex,ey,ez,...
            ux_p_lid,uy_p_lid,uz_p_lid);
        % Collide
        switch dynamics
            
            case 1 % LBGK
                for i = 1:numSpd
                    cu = 3*(ex(i)*ux+ey(i)*uy+ez(i)*uz);
                    fEq(:,i)=w(i)*rho.*(1+cu+(1/2)*(cu.*cu) - ...
                        (3/2)*(ux.^2 + uy.^2+uz.^2 ));
                    fOut(:,i)=fIn(:,i)-omega*(fIn(:,i)-fEq(:,i));
                end
                
            case 2 % TRT
                
                % find even part and odd part of off-equilibrium
                % for all speeds, then relax those parts...
                
                % compute fEq
                for i = 1:numSpd
                    cu = 3*(ex(i)*ux+ey(i)*uy+ez(i)*uz);
                    fEq(:,i)=w(i)*rho.*(1+cu+(1/2)*(cu.*cu) - ...
                        (3/2)*(ux.^2 + uy.^2+uz.^2 ));
                end
                
                % compute fNEq
                fNEq = fEq - fIn;
                
                % compute odd and even part of fNEq
                for i= 1:numSpd
                    fEven(:,i) = (1/2)*(fNEq(:,i)+fNEq(:,bb_spd(i)));
                    fOdd(:,i) = (1/2)*(fNEq(:,i)-fNEq(:,bb_spd(i)));
                    
                end
                
                % relax on these even and odd parts
                fOut = fIn + omega*fEven + fOdd;
                % omega for odd part is equal to 1.
                
            case 3 % MRT
                % compute fEq
                for i = 1:numSpd
                    cu = 3*(ex(i)*ux+ey(i)*uy+ez(i)*uz);
                    fEq(:,i)=w(i)*rho.*(1+cu+(1/2)*(cu.*cu) - ...
                        (3/2)*(ux.^2 + uy.^2+uz.^2 ));
                end
                % collide
                fOut = fIn - (fIn - fEq)*omega_op;
                
                
        end
        
         % bounce-back
        for i = 1:numSpd
            fOut(snl,i)=fIn(snl,bb_spd(i));
        end
        
        
        % stream
        for i = 1:numSpd
            fIn(streamTgtMat(:,i),i)=fOut(:,i);
        end
        
        % visualization        
        if(mod(ts,plot_freq)==0)
           figure(1)
            % velocity magnitude on mid-box slice
            ux_vp = ux(vis_nodes);
            uy_vp = uy(vis_nodes);
            uz_vp = uz(vis_nodes);
            u_vp = sqrt(ux_vp.^2+uy_vp.^2+uz_vp.^2)./u_conv_fact;
            u_vp =reshape(u_vp,[Ny Nx]);
            %imagesc(u_vp);
            contourf(u_vp,30);
            colorbar
            title('Velocity contour at mid-cavity slice');
            xlabel('x');
            ylabel('y');
            axis equal off
            
            figure(2)
            % density magnitude on mid-box slice
            rho_p = rho(vis_nodes);
            rho_p = reshape(rho_p,[Ny Nx]);
            contourf(rho_p,150);
            colorbar
            title('Density contour at mid-cavity slice');
            xlabel('x');
            ylabel('y');
            axis equal off
            
            % for figure 3, include pressure on bottom plane
            figure(3)
            rho_p = rho(faces.zy_p);
            rho_p = reshape(rho_p,[Nz Ny]);
            contourf(rho_p,50)
            colorbar
            title('Density Contour along bottom plane');
            xlabel('z');
            ylabel('y');
            axis equal off
            colorbar
            
            drawnow
            
            
        end
        
    end
    ex_time = toc;
    
    fprintf('Lattice-point updates per second = %g.\n',...
        Num_ts*nnodes/ex_time);
    %profile viewer
    
    
    
else
    fprintf('Run aborted.  Better luck next time!\n');
end
