function data_pack_exact = analytical_heat_1d(nx,D,C0,t_steps,dt,freq_results,x,Lx)

% Analytical solution (2d diffusion equation, steady-state)

% create mesh
Lx = 1;                       %x-size [m]
t = dt*freq_results*t_steps;  %t-vector
n  = 1000;                    %number of Fourier modes
%x = 0:dx:Lx;                  %x-vector
n_steps = length(t_steps);    %number of time steps


field1  = 'field';    value1  = zeros(nx,n_steps);
field2  = 'x';        value2  = zeros(nx,1);
field3  = 't';        value3  = zeros(n_steps,1);

data_pack_exact = struct(field1,value1,field2,value2,field3,value3);


% main loop (computes C(x,t) exact; i = x rank; j = t rank)
for j = 1:n_steps
    for i=1:nx
        sum = 0;
        for k=1:n
            beta = k*pi/Lx;
            sum = sum - 2*C0/(beta*Lx) * sin(beta*x(i)) * exp(-D*t(j)*beta^2);
        end
        data_pack_exact.field(i,j) = C0 * (-x(i)/Lx + 1) + sum;
    end
end


data_pack_exact.x = x;
data_pack_exact.t = t;

end