%backward_step_3D.m


clear
clc
close('all');

initialization = 0;
% 0 = initialize fIn to zero speed
% 1 = initialize fIn to Poiseuille profile <--not used for this problem

dynamics = 2;
% 1 = LBGK
% 2 = TRT 
% 3 = MRT <--- D3Q15 and D3Q19 only

lattice_selection = 3;
% 1 = D3Q15
% 2 = D3Q19
% 3 = D3Q27

fluid = 1;
% 1 = glycerin
% 2 = glycol
% 3 = water

Num_ts = 10000;
ts_rep_freq = 100;
plot_freq = 50;
Re = 10;

Lx_p = 5;
Ly_p = 0.5;
Lz_p = 0.5;

step_height_frac = 1/3;
step_length_frac = 1/5;

step_height = Ly_p*step_height_frac;
step_x_pos = Lx_p*step_length_frac;

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

Lo = step_height;
Uavg = nu_p*Re/Lo;

% convert to dimensionless units
%Lo = Ly_p;

Uo = Uavg;
To = Lo/Uo;

Ld = 1; Td = 1; Ud = (To/Lo)*Uavg;
nu_d = 1/Re;

% convert to LBM units
dt = 5e-3;
Ny_divs = 6;
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



% get the solid nodes.  The flow will flow in the zy_m face and out the
% zy_p face.  All the other faces will be solid nodes.

snl = [faces.zx_m; faces.zx_p; faces.xy_m; faces.xy_p];
snl = unique(snl);

obst_node_list = find((gcoord(:,1)<(step_x_pos)) & (gcoord(:,2)<(step_height)));
snl = [snl;obst_node_list];

inl = find(gcoord(:,1)==0 & ~(gcoord(:,2)<(step_height) | gcoord(:,2)==(Ly_p)));
inl = setxor(inl,intersect(inl,snl)); % eliminate solid nodes from inl
onl = faces.zy_p;
onl = setxor(onl,intersect(onl,snl)); %eliminate solid nodes from onl



% compute parabolic inlet/outlet velocity boundary condition
Umax_in = (3/2)*u_lbm;
%by = Ly_p/2;
in_cl = step_height + (Ly_p - step_height)/2;
in_wdth=(Ly_p-step_height)/2;
ux_p_in = Umax_in*(1-((gcoord(inl,2)-in_cl)/in_wdth).^2);
uy_p_in = zeros(length(ux_p_in),1);
uz_p_in = uy_p_in;

% ux_p_in = u_lbm*ones(length(inl),1);
% uy_p_in = zeros(length(inl),1);

Umax_out = (3/2)*u_lbm*(2*in_wdth/Ly_p);
ux_p_out = Umax_out*(1-((gcoord(onl,2)-(Ly_p/2))/(Ly_p/2)).^2);
uy_p_out = zeros(length(onl),1);
uz_p_out = uy_p_out;


% tag some nodes for visualization
% xy-plane along the centerline to get an idea of the overall flow pattern
z_pln = z_space(ceil(Nz/2));
vis_nodes = find(gcoord(:,3)==z_pln);

% zy-plane at mid-channel to get a view of the flow profile
zy_pln = ind(:,ceil(Nx/2),:);
zy_pln = reshape(zy_pln,[Ny Nz]);

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
        ux(inl)=ux_p_in;
        uy(inl)=uy_p_in;
        uz(inl)=uz_p_in;
        
        ux(onl)=ux_p_out;
        uy(onl)=uy_p_out;
        uz(onl)=uz_p_out;
        
        % microscopic BCs
        fIn(inl,:)=velocityBC_3D(fIn(inl,:),w,ex,ey,ez,...
            ux_p_in,uy_p_in,uz_p_in);
        fIn(onl,:)=velocityBC_3D(fIn(onl,:),w,ex,ey,ez,...
            ux_p_out,uy_p_out,uz_p_out);
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
            ux_vp = ux(vis_nodes);
            uy_vp = uy(vis_nodes);
            uz_vp = uz(vis_nodes);
            u_vp = sqrt(ux_vp.^2+uy_vp.^2+uz_vp.^2)./u_conv_fact;
            u_vp =reshape(u_vp,[Ny Nx]);
            imagesc(u_vp);
            colorbar
            title('Midplane velocity profile')
            axis equal off
            drawnow
            
            figure(2)
            ux_cp = ux(zy_pln(:));
            uy_cp = uy(zy_pln(:));
            uz_cp = uz(zy_pln(:));
            u_cp = sqrt(ux_cp.^2+uy_cp.^2+uz_cp.^2)./u_conv_fact;
            u_cp = reshape(u_cp,[Ny Nz]);
            mesh(u_cp)
            view(3)
            title('Mid-channel velocity profile')
            drawnow
            
%             figure(3)
%            ux_quiv = reshape(ux,[Ny Nx Nz]);
%            uy_quiv = reshape(uy,[Ny Nx Nz]);
%            uz_quiv = reshape(uz,[Ny Nx Nz]);
%            quiver3(X,Y,Z,ux_quiv,uy_quiv,uz_quiv)
%            drawnow
%            pause(0.05);
            
            
        end
        
    end
    ex_time = toc;
    
    fprintf('Lattice-point updates per second = %g.\n',...
        Num_ts*nnodes/ex_time);
    
    
    
    
else
    fprintf('Run aborted.  Better luck next time!\n');
end