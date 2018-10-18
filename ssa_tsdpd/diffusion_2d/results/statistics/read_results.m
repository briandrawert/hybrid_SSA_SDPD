function [C_SSA,C_SDPD] = read_results(nx,ny,nt,data)

C_SSA  = zeros(nx,ny,nt);
C_SDPD = zeros(nx,ny,nt);

for k = 1:nt
    for j=1:ny
        for i=1:nx
            C_SDPD(i,j,k) = data(i+(j-1)*ny,1,k);
            C_SSA(i,j,k)  = data(i+(j-1)*ny,2,k);
        end
    end
end

end
