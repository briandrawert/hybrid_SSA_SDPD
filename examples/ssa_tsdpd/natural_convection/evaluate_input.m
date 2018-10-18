function evaluate_input

% Density [kg/m^3]
rho = 1;


% Lengths
Lxint = 1;
Lyint = 1;


% Number of particles in x and y
Nxint = 51;
Nyint = 51;


% Number of wall particles per wall in x and y
Nxwall = 6;
Nywall = 6;


% Particle spacing in x, y and z directions
dx = Lxint/Nxint;
dy = Lyint/Nyint;
fprintf('dx            = %.16f [m] \n',dx);
fprintf('dy            = %.16f [m] \n',dy);


% Total length in x and y (including walls)
Lx = Lxint + Nxwall*2*dx;
Ly = Lyint + Nywall*2*dy;
fprintf('Lx            = %.16f [m] \n',Lx);
fprintf('Ly            = %.16f [m] \n',Ly);


% Total number of particles
Npx = Nxint + 2*Nxwall;
Npy = Nyint + 2*Nywall;
Np = Npx*Npy;
fprintf('Np            = %d [particles] \n',Np);


% Total volume
vtot = Lx*Ly*1;
fprintf('vtot          = %d [particles] \n',vtot);


% Volume per particle
vi = vtot / Np;
fprintf('vi           = %f [m^3] \n', vi);


% Total mass
mtot = vtot*rho;
fprintf('mtot         = %f [kg] \n',mtot);


% Mass per particle
mi = vtot/Np;
fprintf('mi           = %.16f [kg] \n', mi);


%compute imposed C based on imposed Cd
Cd = 1000;
C = Cd/vi;
fprintf('C               = %f [molecules/m^3] \n', C );

end