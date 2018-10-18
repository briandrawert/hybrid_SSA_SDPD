function data_pack_fd = fd_cylinder(nx,D,C0,t_solution,dt,freq_results,x)

% space specs
Lx = 1;
dx = 0.1/4;
x = (0:dx:Lx)';
nx = length(x);

% time specs
Lt = 1;
tvec = (0:dt:Lt)';
t = 0;
nt = length(t_solution);


field1  = 'C0';           value1  = zeros(nx,nt);
field2  = 'C1';           value2  = zeros(nx,nt);
field3  = 'x';            value3  = zeros(nx,1);

data_pack_fd = struct(field1,value1,field2,value2,field3,value3);

data_pack_fd.x = x;
counter = 0;


% physical parameters
k = 0.1;
k_AB = 0.1;

% ICs
CA    = zeros(nx,1);
CB    = zeros(nx,1);
CAnew = zeros(nx,1);
CBnew = zeros(nx,1);


% BCs
CA(1)  = C0;
CA(nx) = 0;
CB(1)  = 0;
CB(nx) = C0;


CAnew(1)  = C0;
CAnew(nx) = 0;
CBnew(1)  = 0;
CBnew(nx) = C0;



for tc = 1:length(tvec)
    
    for i=2:nx-1
        div_CA = (CA(i+1) - 2* CA(i) + CA(i-1))/(dx^2);
        div_CB = (CB(i+1) - 2* CB(i) + CB(i-1))/(dx^2);
        R_A    = CA(i)*CB(i);
        R_B    = R_A;
        
        CAnew(i) = dt*k*div_CA - k_AB*dt*R_A + CA(i);
        CBnew(i) = dt*k*div_CB - k_AB*dt*R_B + CB(i); 
    end    

    
    CA = CAnew;
    CB = CBnew;

    if (mod(tc,freq_results) == 0)
        counter = counter + 1;
        data_pack_fd.C0(:,counter) = CA;
        data_pack_fd.C1(:,counter) = CB;
    end
    
end




end