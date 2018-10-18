function evaluate_input

%number of dimensions
% note: if n_dimensions = 2, Lz = 1 by definition (pseudo-2D)
n_dimensions = 1;


%distance between particles [m]
dist = 1/16;
fprintf('dist            = %.16f [m] \n',dist);


%dimensions
Lx=[0,1+dist-1e-15];
Ly=[0,dist-1e-15];
Lz=[0,dist-1e-15];

%physical parameters [kg/m^3]
density = 1000;


%number of particles in x, y and z directions
nx = floor((Lx(2)-Lx(1))/dist) + 1;
ny = floor((Ly(2)-Ly(1))/dist) + 1;
nz = ceil((Lz(2)-Lz(1))/dist);

%total number of particles
if (n_dimensions == 3)
    np = nx*ny*nz;
elseif (n_dimensions == 2)
    np = nx*ny;
elseif (n_dimensions == 1)
    np = nx;
end
fprintf('np              = %d [particles] \n',np);

%compute total volume
if (n_dimensions == 3)
    volume_total = (Lx(2)-Lx(1))*(Ly(2)-Ly(1))*(Lz(2)-Lz(1));
    fprintf('volume_total    = %f [m^3] \n',volume_total);
elseif (n_dimensions == 2)
    volume_total = (Lx(2)-Lx(1))*(Ly(2)-Ly(1))*1.0;
    fprintf('volume_total    = %f [m^3] \n',volume_total);
elseif (n_dimensions == 1)
    volume_total = (Lx(2)-Lx(1))*1.0*1.0;
    fprintf('volume_total    = %f [m^3] \n',volume_total);
end

%compute total mass
mass_total = density * volume_total;
fprintf('mass_total      = %f [kg] \n',mass_total);

%compute mass per particle
mass_particle = mass_total / np;
fprintf('mass_particle   = %.16f [kg] \n', mass_particle);

%compute volume per particle
volume_particle = volume_total / np;
fprintf('volume_particle = %f [m^3] \n', volume_particle);

%compute imposed C based on imposed Cd
Cd = 1000;
C = Cd/volume_particle;
fprintf('C               = %f [molecules/m^3] \n', C );

end