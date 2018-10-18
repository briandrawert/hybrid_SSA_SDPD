function compact_fd_cylinder

%==========================================================================
% This program solves the system of coupled partial differential equations:
% 
%          dy1/dt = k * d2(y1)/dx2 - k12 * R1
%          dy2/dt = k * d2(y2)/dx2 - k12 * R2
%
% using compact finite differences in the space and a fourth order
% Runge-Kutta method in time.
%==========================================================================

%==========================================================================
% Define physical parameters
%==========================================================================

% space specs
Lx = 1;
dx = 0.1/4;
x = (0:dx:Lx)';
nx = length(x);

% time specs
Lt = 1;
dt = 0.0001;
t = (0:dt:Lt)';
nt = length(t);
freq_save = 2000; %saves results

% physical parameters
k   = 1e-1;
k12 = 1000;
neq = 2;

% compact fd coefficients
alpha = 2/11;
a_fd  = 12/11;
b_fd  = 3/11;


%==========================================================================
% Initial and boundary conditions
%==========================================================================
% ICs
y = zeros(nx,neq);

% BCs
y(1,1)  = 1;
y(1,2)  = 0;

y(nx,1) = 0;
y(nx,2) = 1;




%==========================================================================
% Form matrix A to solve 2nd order derivatives:
%==========================================================================
A = zeros(nx,nx);
A = compute_A_matrix(A,nx,alpha);



%==========================================================================
% Solve system of ODEs in time
%==========================================================================        
for tc = 1:nt
    
    %----------------------------------------------------------------------
    % Evaluate k1_RK
    %----------------------------------------------------------------------
    k1_RK = f(t,y,nx,dx,a_fd,b_fd,neq,k12,k,A);
    
    
    %----------------------------------------------------------------------
    % Evaluate k2_RK
    %----------------------------------------------------------------------
    k2_RK = f(t+dt/2,y+dt*k1_RK/2,nx,dx,a_fd,b_fd,neq,k12,k,A);
    
    
    %----------------------------------------------------------------------
    % Evaluate k3_RK
    %----------------------------------------------------------------------
    k3_RK = f(t+dt/2,y+dt*k2_RK/2,nx,dx,a_fd,b_fd,neq,k12,k,A);
   
    
    %----------------------------------------------------------------------
    % Evaluate k4_RK
    %----------------------------------------------------------------------
    k4_RK = f(t+dt,y+dt*k3_RK,nx,dx,a_fd,b_fd,neq,k12,k,A);
     
    %----------------------------------------------------------------------
    % Step to time t+dt
    %----------------------------------------------------------------------
    y = y + (dt/6) * (k1_RK + 2*k2_RK + 2*k3_RK + k4_RK);
     



    if (mod(tc,freq_save) == 0)
        figure(1);
        hold on;
        plot(x,y(:,1),'b');
        plot(x,y(:,2),'r');
        legend('C_A','C_B');
    end



end



% % Analytical solution (1d diffusion equation)
% 
% % Time control
% Lt        = 1;             %time span
% t         = Lt;            %time
% C0        = 1;             %concentration at right boundary
% n         = 200;           %number of Fourier modes
% C         = zeros(1,nx);
% Lx        = 1;
% D         = k;
% 
% for i=1:nx
%     sum = 0;
%     for k=1:n
%         beta = k*pi/Lx;
%         sum = sum - 2*C0/(beta*Lx) * sin(beta*x(i)) * exp(-D*t*beta^2);
%     end
%     C(i) = C0 * (-x(i)/Lx + 1) + sum;
% end
% 
% 
% 
% Compute L2 norm of solution at steady-state (pure convection case)
% L2_error = norm(C(:) - y(:,1))/nx;
% fprintf('L2 error = %.16f \n',L2_error);
% 

end



function RHS = f(t,y,nx,dx,a_fd,b_fd,neq,k12,k,A)

% pre-allocate variables
b       = zeros(nx,1);
RHS     = zeros(nx,neq);


%--------------------------------------------------------------------------
%loop over species:
%--------------------------------------------------------------------------
for j=1:neq
    %----------------------------------------------------------------------
    %compute b vector:
    %----------------------------------------------------------------------
    
    %boundary and near-boundary nodes
    b(2)    = (13*y(2,j)    - 27*y(3,j)      + 15*y(4,j)    - y(5,j))    /(dx^2);
    b(nx-1) = (13*y(nx-1,j) - 27*y(nx-2,j)   + 15*y(nx-3,j) - y(nx-4,j)) /(dx^2);

    %interior nodes
    for i=3:nx-2
        b(i) = (a_fd/(dx^2)) * (y(i-1,j)-2*y(i,j)+y(i+1,j)) + (b_fd/(4*(dx^2))) * (y(i-2,j)-2*y(i,j)+y(i+2,j));
    end
    
    d2y_dx2     = A(2:nx-1,2:nx-1)\b(2:nx-1);
    RHS(2:nx-1,j) = k*d2y_dx2 - k12 * y(2:nx-1,1) .* y(2:nx-1,2);
end


end



function A = compute_A_matrix(A,nx,alpha)

% Create main matrix structure
for i=2:nx-1
   for j=2:nx-1
       if (i==j) 
           A(i,j) = 1;
       elseif (i+1==j)
           A(i,j) = alpha;
       elseif (i-1==j)
           A(i,j) = alpha;
       end
   end
end

% Fix near-boundary nodes, i=2 and i=nx-1: (4th order, one-directional compact FDs)
A(2,2) = 1;
A(2,3) = 11;
A(nx-1,nx-1) = 1;
A(nx-1,nx-2) = 11;


end


