function [x, uout] = burgers(delta, nu, Mlist)
%   This m-file contains the code for the non-linear convection-diffusion
%   equation using a Finite Volume Method
%   - general boundary conditions: Dirichlet, Neumann, periodic
%   - options for explicit or implicit time marching (theta method)
%   - central discretization convective terms (2nd order)
%   - different linearizations of convective term
%   - remove linearization error at each time step,
%   - not energy conserving because in 1D there is no continuity equation
%
%   'horizontal', i.e. velocity points located on boundary
%   shifted boundary for Neumann to keep number of unknowns the same for
%   different boundary conditions

%   Benjamin Sanderse, March 2010
%   update June 2016



%%  input parameters


% flow properties and boundary conditions
Re      = 1/nu;       % Reynolds number
c       = 1;            % convection speed (can be set to zero for diffusion problems)

BC.left = 'dir';        % 'dir', 'neu', or 'per'. for periodic both sides must be 'per'
BC.right= 'dir';

uLe     = 1+delta;
uRi     = -1;
duLe    = 0; %Re/(exp(-Re)-1);
duRi    = 0; %Re/(1-exp(-Re));

exact_solution = 0; % 1 if exact solution available - used for errors and plotting

plot_solution = 0; % 1 if solution should be shown in a figure

% number of volumes in the x-direction
% when providing a list, all meshes will be evaluated and an error plot
% constructed
% Mlist   = [40 80 160 320 640 1280 2560 5120];
% Mlist   = 100000;

% mesh list iteration counter
n       = 1;

% domain and mesh
xLe     = -1;
xRi     = 1;
sx      = 1.0;                                      % stretch factor for non-uniform meshes
gridtype= 1;

% steady or unsteady problem
steady  = 1;

% time step and end time (for unsteady problems)
dt      = 0.01;
t       = 0.;                                       % start time
t_end   = 0.15;                                     % end time, not needed for steady simulations
    
theta   = 0.5;                                      % time integration method

% nonlinear solver convergence criterion (steady and unsteady)
acc     = 1e-8;
it_max  = 100;
it_newton = 3;

% loop over meshes
for Mi=Mlist
    
    L       = xRi-xLe;                                       % width of domain
    deltax  = L/Mi;                                     % uniform mesh size x-direction
    
    %% mesh
    if (gridtype==1)
        % uniform mesh
        [x,hx]  = nonuniform_grid(deltax,xLe,xRi,sx);
        
        % piecewise
        %     [x1,hx]  = nonuniform_grid(deltax,0,L/2,sx);
        %     x2 = x1+L/2;
        %     x = [x1 x2(2:end)];
        
        if (Re*deltax==2)
            warning('Re*h exactly 2!');
        end
        
    elseif (gridtype==3)
        % exponential mesh - refinement in boundary layer
        
        %     Shishkin point
        xs = xRi+log(1e-2)/Re;
        %     xs          = L*(1-(10/Re));  % corresponding to 2003 paper Veldman
        %     xs          = max(0.5*L,xs);
        ksi         = nonuniform_grid(deltax,xLe,xRi,1); % uniform grid
        sx          = ((1-xs)/xs)^2;  % stretch factor
        if (sx==1)
            warning('exponential grid is uniform!')
            [x,hx]  = nonuniform_grid(deltax,xLe,xRi,L,sx);
        else
            x       = L*(1-sx.^ksi)/(1-sx);
        end
    end
    
    % mesh size
    hx      = diff(x);
    M       = length(x)-1;
    Mlist(n)= M;
    
    % change the mesh for Neumann or periodic conditions
    if (strcmp(BC.left,'dir') && strcmp(BC.right,'neu'))
        hxnew = hx*(x(end)-x(1))/(x(end-1)-x(1));
        x     = [x(1) cumsum(hxnew)];
    end
    if (strcmp(BC.left,'neu') && strcmp(BC.right,'dir'))
        hxnew = hx*(x(1)-x(end))/(x(2)-x(end));
        xl    = x(1)-hxnew(1);
        x     = [xl xl+cumsum(hxnew)];
    end
    if (strcmp(BC.right,'neu') && strcmp(BC.left,'neu'))
        hxnew = hx*(x(end)-x(1))/(x(end-1)-x(2));
        xl    = x(1)-hxnew(1);
        x     = [xl xl+cumsum(hxnew)];
    end
    
    
    if (strcmp(BC.left,'per') && strcmp(BC.right,'per')) % periodic
        % change the domain so that x=L at M+1
        % then x(M+1) corresponds to x(1)
        hxnew = hx*(x(end)-x(1))/(x(end-1)-x(1));
        x     = [x(1) cumsum(hxnew)];
        hx    = diff(x);
        hx(end) = hx(1);
        x(end) = x(end-1)+hx(end);
        
        gx      = zeros(1,M+1);
        gx(1)   = (hx(1)+hx(end))/2;
        gx(2:M) = (hx(1:M-1)+hx(2:M))/2;
        gx(M+1) = (hx(M)+hx(end))/2;
        % sum(gx(2:M)) should be L_x
    else
        hx      = diff(x);
        gx      = zeros(1,M+1);
        gx(1)   = hx(1);
        gx(2:M) = (hx(1:M-1)+hx(2:M))/2;
        gx(M+1) = hx(end);
    end
    
    
    %% matrices for discretization
    
    % matrices with metrics
    diag1   = 1./hx';
    Om_ux   = spdiags(diag1,0,M,M);
    
    diag1   = gx(2:end-1)';
    Om_u    = spdiags(diag1,0,M-1,M-1);
    
    % differencing from 'pressure' points to velocity points
    diag1   = ones(M,1);
    Mu      = spdiags([-diag1 diag1],[0 1],M-1,M);
    
    % differencing from velocity to pressure points
    diag1   = ones(M,1);
    Mut     = spdiags([-diag1 diag1],[0 1],M,M+1);
    
    % averaging from velocity to pressure (like Mut)
    diag1   = 0.5*ones(M,1);
    A       = spdiags([diag1 diag1],[0 1],M,M+1);
    
    
    % boundary matrices:
    % u = Bin uin + Bb ub
    % Bbc u = ybc
    
    Bin     = spdiags(ones(M+1,1),-1,M+1,M-1);
    Bb      = spalloc(M+1,2,2);
    Bb(1,1) = 1;
    Bb(end,2)=1;
    
    % boundary conditions
    Bbc = spalloc(2,M+1,4);
    ybc = zeros(2,1);
    
    if (strcmp(BC.left,'dir') && strcmp(BC.right,'dir'))
        
        Bbc(1,1)   = 1;
        Bbc(2,end) = 1;
        
        ybc(1)     = uLe;
        ybc(2)     = uRi;
    end    
    if (strcmp(BC.left,'neu') && strcmp(BC.right,'dir'))
        
        Bbc(1,1)   = -1;
        Bbc(1,3)   = 1;
        Bbc(2,end) = 1;
        
        ybc(1)     = 2*hx(1)*duLe;
        ybc(2)     = uRi;
    end
    if (strcmp(BC.right,'neu') && strcmp(BC.left,'dir'))
        
        Bbc(1,1)   = 1;
        Bbc(2,end-2)= -1;
        Bbc(2,end) = 1;
        
        ybc(1)     = uLe;
        ybc(2)     = 2*hx(end)*duRi;
    end
    if (strcmp(BC.left,'neu') && strcmp(BC.right,'neu'))
        
        Bbc(1,1)   = -1;
        Bbc(1,3)   = 1;
        Bbc(2,end-2)= -1;
        Bbc(2,end) = 1;
        
        ybc(1)     = 2*hx(1)*duLe;
        ybc(2)     = 2*hx(end)*duRi;
    end
    if (strcmp(BC.right,'per') && strcmp(BC.left,'per'))
        Bbc(1,1)     = 1;
        Bbc(1,end-1) = -1;
        Bbc(2,2)     = 1;
        Bbc(2,end)   = -1;
    end
    
    % substituting u into Bbc u = ybc gives
    % u = (Bin - Bb*(Bbc*Bb)^(-1)*Bbc*Bin)*uin + Bb*(Bbc*Bb)^(-1)*ybc
    %   = B*uin + D*ybc
    D   = Bb*(Bbc*Bb\speye(2));     % = inv(Bbc*Bb)
    B   = Bin - D*Bbc*Bin;
    
    % identity matrix for implicit time
    I = speye(size(B,2));
    
    % interpolation matrix for convecting quantity
    In = A;
    

    %% initial condition (also necessary for steady nonlinear cases)
    % u_start       = zeros(M+1,1);
    
    % sine 1:
    % u_start       = sin(2*pi*x');
    % sine 2: (sufficiently smooth at x=L/4 and x-3*L/4)
    %         u_start       = (x>L/4 & x<3*L/4)'.*sin(pi*(x'-L/4)/(L/2)).^4;
    % sine with exact solution for Burgers equation:
    %     u_start   = (2/Re)*pi*sin(pi*x')./(2+cos(pi*x'));
    % hat:
    %         u_start       = (x>L/4 & x<3*L/4)'.*ones(M+1,1);
    % spike:
    %     u_start(M/2,1)= 1;
    % heaviside: make sure that M = odd, to prevent a point at x=L/2
    % u_start = (x>L/2)'.*ones(M+1,1);
    %     u_start = (x<L/2)'.*ones(M+1,1);
    %     if (mod(M,2)==0)
    %         u_start(M/2+1)= 0.5;
    %     end
    
    % linear:
    u_start = uLe + (uRi-uLe)/(xRi-xLe) * (x - xLe)';
    
    
    % solution at interior nodes
    u =  Bin'*u_start;
    
    
    %% solve steady state or unsteady problem
    
    eps     = dt/100;
       
    i       = 1; % iteration counter or time step counter
    
    % source term g, can be useful for method of manufactured solutions
    xin  = Bin'*x';
    g    = zeros(length(xin),1);
    
    % for example, take g such that exact solution u=(1-x)^p results
    % p    = 4;
    % g    = dt*( -c*p*(1-xin).^(p-1) - (1/Re)*p*(p-1)*(1-xin).^(p-2));
    
    if (steady==1)
        
        % iterative process to solve nonlinear system
        res = 1;
        Newton = 2;
        
        while (res > acc && i<=it_max)
           
            if i < it_newton
                % Picard:
                C    = 0.5*spdiags(In*(B*u+D*ybc),0,M,M); % convective velocity
                Q    = Om_u * ( -c*Mu*C*A + (1/Re)*Mu*Om_ux*Mut);
                R    = Q*B;
                yb   = Q*D*ybc;
                u    = -R\(yb+g);
                res   = 1+acc;
            else
                % Newton iteration (factor 2) 
                f    = (c*0.5*Mu*(A*(B*u+D*ybc)).^2 - (1/Re)*Mu*Om_ux*Mut*(B*u+D*ybc)) + Om_u*g;
                du   = - (c*0.5*Mu*Newton*spdiags(A*(B*u+D*ybc),0,M,M)*A - (1/Re)*Mu*Om_ux*Mut)*B \ f;
                u    = u + du;
                res  = max(abs(f));
            end
            i    = i+1;
            
        end

        if (i>=it_max)
            error('Newton solver did not converge, (delta, nu) = %.8e %.8e', delta, nu);
        end
        
    else
        
        % for time-dependent BC, specify something complying with BC in
        % loop (could be made a function)
        % ybc(1)   = sin(0);
        % ybc(end) = cos(0);

        % march in time
        while (t<=t_end-dt+eps)
            
            % Crank-Nicolson: (theta=0: FE, theta=1: BE)
            % time-dependent boundary conditions; check starting values
            %     ybc(1)  = sin(t+dt);
            %     ybc(end)= cos(t+dt);
            
            % linearization methods
            % Picard:
            %         C1    = 0.5*spdiags(In*(B*u+D*ybc),0,M,M); % convective velocity
            %         C2    = C1;
            % Extrapolated Picard:
            %         C1    = 0.5*spdiags(In*(B*u+D*ybc),0,M,M);
            %         C2    = 0.5*spdiags(In*(B*(2*u - u_old)+D*ybc),0,M,M);
            
            % Linear:
            %         C1    = speye(M);
            %         C2    = C1;
            
            u_old = u;
            du    = 1;
            
            % iterative process to solve nonlinear system
            while (du > acc)
                % Newton:
                C1   = c*Mu*spdiags(In*(B*u_old+D*ybc),0,M,M)*A;
                C2   = c*Mu*spdiags(In*(B*(u_old+u)+2*D*ybc),0,M,M)*A;         % implicit part
                C3   = c*Mu*spdiags(In*(B*u+D*ybc),0,M,M)*A;
                %             Diff = (1/Re)*Mu*Om_ux*Mut;
                Diff = spalloc(M-1,M+1,0);
                Q1   = (-(1/8)*C1 + 0.5*Diff); % explicit part
                Q2   = ( (1/4)*C2 - 0.5*Diff); % implicit part
                R1   = Q1*B;
                R2   = Q2*B;
                yb1  = Q1*D*ybc;
                yb2  = Q2*D*ybc;
                
                f      = (Om_u/dt + R1) *u_old + yb1 + yb2 + (1/8)*C3*(B*u+D*ybc);
                u_new  = (Om_u/dt + R2)\f;
                u      = u_new;
                
                % residual in finite volume setting
                u_temp  = B*(u_old + u_new)/2 + D*ybc;
                du      = max(abs( Om_u*(u_new - u_old)/dt + ...
                    ( 0.5*c*Mu*spdiags(In*u_temp,0,M,M)*A*u_temp - ...
                    Diff*u_temp) ));
                %             i    = i+1;
                %             figure(1)
                %             plot(xin,u_new,'bx-')
                %             pause;
            end
            
            %         plot total solution
            figure(1)
            plot(xin,u_new,'bx-')
            
            k(i)  = sum(u.^2.*gx(2:end-1)');
            q(i)  = sum(u.*gx(2:end-1)');
            
            i     = i+1;
            t     = t+dt; % time level of solution u just computed
            
        end
    end
    
    % compare with exact solution if available
    if (exact_solution==1)
        
        % boundary-layer:
        %  u_ex = u_exact_BL(xin,uLe,uRi,Re);
        %u_ex = u_exact_nonlinear_BL(xin,uLe,uRi,Re);
        u_ex = u_exact_nonlinear_BL2(xin,uLe,uRi,Re);
        
        % sine:
        % u_ex = sin(pi*x'/L)*exp(-pi^2*t/(Re*L^2));
        % transported sine:
        % u_ex = u_start;
        % manufactured solution:
        % u_ex = (1-xin).^p;
        % Heaviside:
        % u_ex = u_exact_linear_Heaviside(xin,t,c,Re,L);
        % u_ex = u_exact_nonlinear_Heaviside(xin,t,Re,L);
        % Burgers:
        % u_ex = u_exact_burgers(xin,t,Re);
        
    end
    
    % solution at interior points
    uh       = u;
    % solution including boundary points
    uout     = B*u+D*ybc;
    
    % compute max error and its position if exact solution is available
    if (exact_solution==1)
        e        = uh-u_ex;
        [errori(n) pos(n)] = max(abs(e));
        error2(n) = ((1/(M-1))*sum(abs(e.^2)))^(1/2);
    end
    
    if (plot_solution==1)
        
        figure(1)
        plot(xin,uh,'bx-')
        
        if (exact_solution)
            hold on
            plot(xin,u_ex,'r')
        end
    end
    
    n = n+1;
    
end

% plot error graphs in loglog plot
if (length(Mlist)>1 && exact_solution==1)
    figure(2)
    loglog(1./Mlist,errori,'bx-');
    hold on
    loglog(1./Mlist,error2,'rx-');
end
end

function [z,dz] = nonuniform_grid(deltaz,z_low,z_up,sz)
% generate a non-uniform grid, from z_low to z_up, starting deltaz and
% having stretch factor close to sz

% version December 2009, checks for uniform grids

if (sz==1)
    
    %     z  = z_low:deltaz:z_up;
    z = linspace(z_low,z_up,floor((z_up-z_low)/deltaz)+1);
    dz = diff(z);
    return;
end
if (sz<1)
    error('stretch factor < 1 not implemented');
end


i    = 1;
z(i) = z_low;

% first check the number of grid points by using the specified sz
while (z(i)+deltaz*(sz^i))<=z_up
    i    = i+1;
    z(i) = z(i-1) + deltaz*(sz^(i-2));
end

% adapt sz such that the upper boundary is exactly satisfied
% sum of powerseries should be z_up-z_low
S = z_up-z(1);           % sum wanted
n = length(z)-1;         % number of intervals

%Secant method for nonlinear iteration
s(1) = sz;
s(2) = sz^2;
i    = 1;
eps  = 1e-15;
while (abs(s(i+1)-s(i))>eps)
    %     s(i+2) = s(i+1) - (s(i+1)-s(i))*(S*(s(i+1)-1) - deltaz*s(i+1)*(s(i+1)^n-1)) / ...
    %                       ( (S*(s(i+1)-1) - deltaz*s(i+1)*(s(i+1)^n-1)) - (S*(s(i)-1) - deltaz*s(i)*(s(i)^n-1)));
    s(i+2) = s(i+1) - (s(i+1)-s(i))*( S*(s(i+1)-1) - deltaz*(s(i+1)^n-1) ) / ...
        ( (S*(s(i+1)-1) - deltaz*(s(i+1)^n-1)) - (S*(s(i)-1) - deltaz*(s(i)^n-1)));
    
    
    i=i+1;
end
sz = s(end);

if (abs(sz-1)<=eps) % sz=1
    %     z  = z_low:deltaz:z_up;
    z = linspace(z_low,z_up,floor((z_up-z_low)/deltaz)+1);
    dz = diff(z);
    return;
end


%check
% deltaz*(sz^n-1)/(sz-1)


%update z now based on the new sz
i    = 1;
z(i) = z_low;
while (i<n+1)
    z(i+1)  = z(i)+deltaz*(sz^(i-1));
    i       = i+1;
end

dz   = diff(z);
% for k=1:length(z)-1
%     dz(k)=z(k+1)-z(k);
% end

end
