%ldc_2D.m

clear
clc
close('all');

dynamics = 1;
% 1 = LBGK
% 2 = TRT 
% 3 = MRT

Num_ts = 10000;
ts_rep_freq = 50;
plot_freq = 100;

Lx_p = 1;
Ly_p = 1;

rho_p = 965.3;
nu_p = 0.06/rho_p;

Re = 1000;
Uavg = nu_p*Re/Ly_p;

% convert to dimensionless units
Lo = Ly_p;
Uo = Uavg;
To = Lo/Uo;

Ld = 1; Td = 1; Ud = (To/Lo)*Uavg;
nu_d = 1/Re;

% convert to LBM units
dt = 5e-4;
Ny_divs = 151;
dx = 1/(Ny_divs-1);
u_lbm = (dt/dx)*Ud;
nu_lbm=(dt/(dx^2))*nu_d;
omega = get_BGK_Omega(nu_lbm);

u_conv_fact = (dt/dx)*(To/Lo);
t_conv_fact = (dt*To);

rho_lbm = rho_p;
rho_out = rho_lbm;

% generate LBM lattice
xm = 0; xp = Lx_p;
ym = 0; yp = Ly_p;

Ny = Ny_divs;
Nx = ceil((Ny_divs-1)*(Lx_p/Ly_p))+1;

[gcoord,~,~]=RecMesh(xm,xp,ym,yp,Nx,Ny);
[nnodes,~]=size(gcoord);

x_space = linspace(xm,xp,Nx);
y_space = linspace(ym,yp,Ny);

[X,Y]=meshgrid(x_space,y_space);

fprintf('dx = %g, dy = %g \n',x_space(2)-x_space(1),...
    y_space(2)-y_space(1));

[w,ex,ey,bb_spd]=D2Q9_lattice_parameters();
%stream_tgt = genTargetVecD2Q9r2(Nx,Ny);
LatticeSize = [Nx Ny];
LatticeSpeeds = [ex; ey];
stm = genStreamTgtMat(LatticeSize,LatticeSpeeds);

numSpd=9;
M = getMomentMatrix('D2Q9');

switch dynamics
    
    case 1
        S = omega.* eye(numSpd);
    case 2
        % TRT model as described in Kevin Tubbs' dissertation section 3.4
        S = zeros(numSpd);
        S(2,2)= omega;
        S(3,3)=omega;
        S(8,8)=omega;
        S(9,9)=omega;
        
        t_s = (1/2)+1/(12*((1/omega)-0.5));
       
        S(5,5)=1/t_s;
        S(7,7)=1/t_s;
    case 3
        % parameters taken from 
        % Chinese Physics Vol 15 No 8 Aug 2006
        % Simulating high Reynolds number flow in 2D lid driven cavity by
        % MRT etc...
        S = zeros(numSpd);
        S(2,2)=1.1;
        S(3,3)=1.0;
        S(5,5)=1.2;
        S(7,7)=1.2;
        S(8,8)=omega;
        S(9,9)=omega;
        
        
end

omega_op = M\(S*M);

% get solid node-list
snl=find((gcoord(:,1)==xp) | (gcoord(:,1)==xm) | (gcoord(:,2)==ym));

% get lid node-list
lnl=find((gcoord(:,2)==yp)); 
%exclude solid nodes from the lid node-list
lnl=setxor(lnl,intersect(snl,lnl));

% get bottom node-list
bnl=find(gcoord(:,2)==ym);

% get side-node lists...
side_nl=find((gcoord(:,1)==xp) | (gcoord(:,1)==xm));

% exclude side-nodes from the bottom node list
bnl = setxor(bnl,intersect(bnl,side_nl));

% prescribed velocity for the lid
ux_p = u_lbm*ones(length(lnl),1);
uy_p = zeros(length(lnl),1);

bottom_x = gcoord(bnl,1);

[fIn,fOut,rho,ux,uy]=Initialize_F_zero(gcoord,ex,ey,w,rho_lbm);

fEq = zeros(nnodes,numSpd);


fprintf('Number of Lattice-points = %d.\n',nnodes);
fprintf('Number of time-steps = %d. \n',Num_ts);
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
        
              
        % set macroscopic and Microscopic Dirichlet-type boundary
        % conditions
        % set macroscopic BC
        ux(snl)=0; uy(snl)=0;
        ux(lnl)=u_lbm; uy(lnl)=0;
       
        
         % set microscopic BC
        fIn(lnl,:)=velocityBC_D2Q9(fIn(lnl,:),w,ex,ey,ux_p,uy_p);
        
        % compute equilibrium
        for i = 1:numSpd
            cu = 3*(ex(i)*ux+ey(i)*uy);
            fEq(:,i)=w(i)*rho.*(1+cu+(1/2)*(cu.*cu) - ...
                (3/2)*(ux.^2 + uy.^2 ));
        end
        
        % Collide
        fOut= fIn - (fIn - fEq)*omega_op;
        
        
        % bounce-back
        for i = 1:numSpd
            fOut(snl,i)=fIn(snl,bb_spd(i));
        end
        
        
        % stream
        %fIn(stream_tgt)=fOut(:);
        for i = 1:numSpd
            fIn(stm(:,i),i)=fOut(:,i);
        end
        
        
        
        if(mod(ts,plot_freq)==0)
            % velocity
            %figure(1)
            figure(1)
            subplot(2,1,1)
            
            u_p = sqrt(ux.*ux+uy.*uy);
            u_p = reshape(u_p,Nx,Ny)./u_conv_fact;
            imagesc(u_p');
            %colorbar
            title('Velocity')
            axis equal off
            % density
            %figure(2)
            subplot(2,1,2)
            rho_p = reshape(rho,Nx,Ny);
            imagesc(rho_p');
            title('Density')
            axis equal off
            
           
            drawnow
            
            
        end
        
        
        
        
    end
    ex_time = toc;
    fprintf('Lattice-point updates per second = %g.\n',nnodes*Num_ts/ex_time);
    
    
else
    fprintf('Run aborted.  Better luck next time!\n');
end
