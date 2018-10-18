function fd_solution_heat_2d

Lx = 1;
Ly = 1;
nx = 100;
ny = 100;
x = linspace(0,Lx,nx);
y = linspace(0,Ly,ny);
dx = Lx/(nx-1);
dy = Ly/(ny-1);


D = 1e-2;
Tspan = 1;
dt = 1e-3;
nt = Tspan/dt + 1;
save_at = 1000;


%IC
C0 = 2e6;
C = zeros(nx,ny);
C(1:nx,ny)  = C0;

[X,Y] = meshgrid(x,y);


for k = 1:nt
    
    for j=2:ny-1
        for i=2:nx-1
            C(i,j) = C(i,j) + dt*D/(dx^2)*(C(i+1,j)-2*C(i,j)+C(i-1,j)) + dt*D/(dy^2)*(C(i,j+1)-2*C(i,j)+C(i,j-1));
        end
    end
    
    if (k == save_at)
        plot(y,C(nx/2,:));
%         contourf(Y,X,C);
%         colorbar;
%         axis square;
    end
    
    
    
    
end

end
