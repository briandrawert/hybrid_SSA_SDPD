function data_pack_exact = analytical_heat_2d(nx,ny)

% Analytical solution (2d diffusion equation, steady-state)

% create mesh
a  = 1;          %0<= x,y <= a;
Lx = 1;          %x-size
Ly = 1;          %y-size
n  = 100;        %number of Fourier modes
dx = Lx/(nx-1);  %grid spacing in x-dir
dy = Ly/(ny-1);  %grid spacing in y-dir
x = 0:dx:Lx;     %x-vector
y = 0:dy:Ly;     %y-vector


field1  = 'field';    value1  = zeros(nx,ny);
field2  = 'h_line';   value2  = zeros(nx,1);
field3  = 'v_line';   value3  = zeros(1,ny);
field4  = 'X';        value4  = zeros(nx,ny);
field5  = 'Y';        value5  = zeros(nx,ny);
field6  = 'x';        value6  = zeros(nx,1);
field7  = 'y';        value7  = zeros(1,ny);

data_pack_exact = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,field7,value7);

% main loop
for i = 1:nx
    for j=1:ny
       sum = 0;
       for k=1:n
          sum = sum + (2*(1-(-1)^k)/(k*pi)) * sin(k*pi*x(i)/a) * sinh(k*pi*y(j)/a) / (sinh(k*pi));
       end
       data_pack_exact.field(i,j) = sum;
    end
end
 
data_pack_exact.h_line = data_pack_exact.field(:,ceil(ny/2));
data_pack_exact.v_line = data_pack_exact.field(ceil(nx/2),:);

[data_pack_exact.X,data_pack_exact.Y] = meshgrid(x,y);

data_pack_exact.x = y;
data_pack_exact.y = x;

end