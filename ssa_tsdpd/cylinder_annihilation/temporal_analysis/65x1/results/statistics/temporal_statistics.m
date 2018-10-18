function [data_pack_SSA,data_pack_SDPD] = temporal_statistics(t_steps,n_realizations,nx,ny,data_pack)


nt = length(t_steps);

field1  = 'C0_mean';        value1     = zeros(nx,ny,nt);
field2  = 'C0_prime';       value2     = zeros(nx,ny,n_realizations,nt);
field3  = 'C0_variance';    value3     = zeros(nx,ny,nt);
field4  = 'C0_sd';          value4     = zeros(nx,ny,nt,1);

field5  = 'C1_mean';        value5     = zeros(nx,ny,nt);
field6  = 'C1_prime';       value6     = zeros(nx,ny,n_realizations,nt);
field7  = 'C1_variance';    value7     = zeros(nx,ny,nt);
field8  = 'C1_sd';          value8     = zeros(nx,ny,nt,1);

field9  = 'x';              value9     = zeros(nx,1);


data_pack_SSA = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,field7,value7,field8,value8,field9,value9);


field1  = 'C0';           value1  = zeros(nx,ny,nt);
field2  = 'C1';           value2  = zeros(nx,ny,nt);
field3  = 'x';            value3  = zeros(nx,1);

data_pack_SDPD = struct(field1,value1,field2,value2,field3,value3);



% =========================================================================
% SSA statistics
% =========================================================================
% data_pack(r,s,nx,ny,column_id)

% compute mean (1st order moment)
for s=1:nt
    for r = 1:n_realizations
        data_pack_SSA.C0_mean(:,:,s) = data_pack_SSA.C0_mean(:,:,s)+ data_pack(:,:,3,r,s)/n_realizations;
        data_pack_SSA.C1_mean(:,:,s) = data_pack_SSA.C1_mean(:,:,s)+ data_pack(:,:,4,r,s)/n_realizations;
    end
end

for i=1:nx
    data_pack_SSA.x(i)  = data_pack(i,1,5,1,1);
end

% compute fluctuations
for s=1:nt
    for r=1:n_realizations
        data_pack_SSA.C0_prime(:,:,r,s) = data_pack(:,:,3,r,s) - data_pack_SSA.C0_mean(:,:,s);
        data_pack_SSA.C1_prime(:,:,r,s) = data_pack(:,:,4,r,s) - data_pack_SSA.C1_mean(:,:,s);
    end
end

% compute variance (2nd order moment)
for s=1:nt
    for r=1:n_realizations
        for i=1:nx
            for j=1:ny
                data_pack_SSA.C0_variance(i,j,s) = data_pack_SSA.C0_variance(i,j,s) + (data_pack_SSA.C0_prime(i,j,r,s)^2)/(n_realizations);
                data_pack_SSA.C1_variance(i,j,s) = data_pack_SSA.C1_variance(i,j,s) + (data_pack_SSA.C1_prime(i,j,r,s)^2)/(n_realizations);
            end
        end
    end
end

% compute standard deviation
data_pack_SSA.C0_sd = sqrt(data_pack_SSA.C0_variance);
data_pack_SSA.C1_sd = sqrt(data_pack_SSA.C1_variance);




% =========================================================================
% SDPD data
% =========================================================================
data_pack_SDPD.x = data_pack_SSA.x;

r=1;
for s = 1:nt
    for i=1:nx
        for j=1:ny
            data_pack_SDPD.C0(i,j,s) = data_pack(i,j,1,r,s);
            data_pack_SDPD.C1(i,j,s) = data_pack(i,j,2,r,s);
        end
    end
end


end
