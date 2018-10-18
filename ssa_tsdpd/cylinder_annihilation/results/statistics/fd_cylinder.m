function fd_cylinder

% space specs
Lx = 1;
dx = 0.001;
x = 0:dx:Lx';
nx = length(x);

% time specs
Lt = 1;
dt = 0.0001;
t = 0:dt:Lt';
nt = length(t);

% physical parameters
k = 0.1;
k_AB = 0.1;

% ICs
CA = zeros(nx,nt);
CB = zeros(nx,nt);


% BCs
CA(1,:)  = 1;
CA(nx,:) = 0;

CB(1,:)  = 0;
CB(nx,:) = 1;




for j = 1:nt
    
    for i=2:nx-1
        div_CA = (CA(i+1,j) - 2* CA(i,j) + CA(i-1,j))/(dx^2);
        div_CB = (CB(i+1,j) - 2* CB(i,j) + CB(i-1,j))/(dx^2);
        R_A    = CA(i,j)*CB(i,j);
        R_B    = R_A;
        
        CA(i,j+1) = dt*k*div_CA - k_AB*dt*R_A + CA(i,j);
        CB(i,j+1) = dt*k*div_CB - k_AB*dt*R_B + CB(i,j); 
    end    
    
end

figure(1)
hold on
% plot(x,CA(:,20),'r');
% plot(x,CB(:,20),'b');
% 
% 
% plot(x,CA(:,40),'r');
% plot(x,CB(:,40),'b');
% 

plot(x,CA(:,end),'r');
plot(x,CB(:,end),'b');








end