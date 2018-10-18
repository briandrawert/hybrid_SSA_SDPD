function heat1d_exact

% Analytical solution (1d diffusion equation, steady-state)

%==========================================================================
% Time control
Lt        = 10;             %time span
dt        = 1;             %time step
t         = 0:dt:Lt;       %time array
nt        = length(t) - 1; %number of time steps
freq_save = 2;             %freq store results 
%==========================================================================


%==========================================================================
% Physical parameters
C0 = 1;           %concentration at right boundary
D  = 1e-2;          %diffusivity [m^2/s]
%==========================================================================


%==========================================================================
% Geometrical parameters
Lx = 1;             %x-size [m]
dx = 1e-2;          %grid spacing
nx = Lx/dx + 1;     %number of grid points
x  = 0:dx:Lx;       %mesh, x-direction
n  = 100;           %number of Fourier modes
%==========================================================================


%==========================================================================
% Preallocate variables
C       = zeros(1,nx);
counter = 0;

%==========================================================================



%==========================================================================
% main loop (computes C(x,t) exact; i = x rank; j = t rank)
for j = 1:nt 
    for i=1:nx
        sum = 0;
        for k=1:n
            beta = k*pi/Lx;
            sum = sum - 2*C0/(beta*Lx) * sin(beta*x(i)) * exp(-D*t(j+1)*beta^2);
        end
        C(i) = C0 * (-x(i)/Lx + 1) + sum;
    end
    if (mod(j,freq_save) == 0)
        counter  = counter + 1;
        if (counter == 1)
            figure(1);
            hold on;
        end
        plot(x,C);
    end
end
%==========================================================================





end