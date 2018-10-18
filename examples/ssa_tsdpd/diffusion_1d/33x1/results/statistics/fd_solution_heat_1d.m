function fd_solution_heat_1d

Lx = 1;
dx = 1e-2;
x = 0:dx:Lx;
nx = length(x);

D = 1e-2;
Tspan = 5;
dt = 1e-3;
nt = Tspan/dt + 1;


%IC
C0 = 1320000;
C = zeros(nx,nt);
C(nx,:)  = 0;
C(1,:)   = C0;

for j=1:nt-1
    for i=2:nx-1
        C(i,j+1) = C(i,j) + dt*D/(dx^2)*(C(i+1,j)-2*C(i,j)+C(i-1,j));
    end
end

figure(1);
plot(x,C(:,end)/C0);
hold on;
end
