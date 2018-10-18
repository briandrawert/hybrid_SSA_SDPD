function evaluate_input

% Density [kg/m^3]
rho = 1000;

eps = 1e-16;

% Lengths
Lxint = 1+eps;
Lyint = 1/34;

% Number of particles in x and y
Nxint = 35;
Nyint = 7;

% Number of wall particles per wall in x and y
Nxwall = 0;
Nywall = 0;

% Particle spacing in x, y and z directions
dx = Lxint/(Nxint-1);
dy = Lyint/(Nyint-1);
fprintf('dx            = %.16f [m] \n',dx);
fprintf('dy            = %.16f [m] \n',dy);

% Total length in x and y (including walls)
Lx = Lxint + Nxwall*2*dx;
Ly = Lyint + Nywall*2*dy;
fprintf('Lx            = %.16f [m] \n',Lx);
fprintf('Ly            = %.16f [m] \n',Ly);

% Total number of particles (interior)
Npx   = Nxint + 2*Nxwall;
Npy   = Nyint + 2*Nywall;
Np    = Npx*Npy;
Npint = Nxint*Nyint;
fprintf('Npint         = %d [particles] \n',Npint);

% Total volume (interior)
vint = Lxint*1;
fprintf('vint          = %d [m^3] \n',vint);

% Volume per particle (interior)
vi = vint / (Npint-1);
fprintf('vi            = %f [m^3] \n', vi);

% Total mass (interior)
mint = vint*rho;
fprintf('mint          = %f [kg] \n',mint);

% Mass per particle (interior)
mi = mint/(Npint-1);
fprintf('mi            = %.16f [kg] \n', mi);

%compute imposed C based on imposed Cd
Cd = 1000;
C = Cd/vi;
fprintf('C             = %f [molecules/m^3] \n', C );

end