function [data_pack_SSA,data_pack_SDPD] = temporal_statistics(t_steps,n_realizations,nx,ny,data_pack)


nt = length(t_steps);

field1  = 'field_mean';        value1     = zeros(nx,ny,nt);
field2  = 'field_prime';       value2     = zeros(nx,ny,n_realizations,nt);
field3  = 'field_variance';    value3     = zeros(nx,ny,nt);
field4  = 'field_sd';          value4     = zeros(nx,ny,nt,1);
field5  = 'x';                 value5     = zeros(nx,1);


data_pack_SSA = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5);


field1  = 'field';        value1  = zeros(nx,ny,nt);
field2  = 'x';            value2  = zeros(nx,1);

data_pack_SDPD = struct(field1,value1,field2,value2);



% =========================================================================
% SSA statistics
% =========================================================================
% data_pack(r,s,nx,ny,column_id)

% compute mean (1st order moment)
for s=1:nt
    for r = 1:n_realizations
        data_pack_SSA.field_mean(:,:,s) = data_pack_SSA.field_mean(:,:,s)+ data_pack(:,:,2,r,s)/n_realizations;
    end
end

for i=1:nx
    data_pack_SSA.x(i)  = data_pack(i,1,3,1,1);
end

% compute fluctuations
for s=1:nt
    for r=1:n_realizations
        data_pack_SSA.field_prime(:,:,r,s) = data_pack(:,:,2,r,s) - data_pack_SSA.field_mean(:,:,s);
    end
end

% compute variance (2nd order moment)
for s=1:nt
    for r=1:n_realizations
        for i=1:nx
            for j=1:ny
                data_pack_SSA.field_variance(i,j,s) = data_pack_SSA.field_variance(i,j,s) + (data_pack_SSA.field_prime(i,j,r,s)^2)/(n_realizations);
            end
        end
    end
end

% compute standard deviation
data_pack_SSA.field_sd = sqrt(data_pack_SSA.field_variance);





% =========================================================================
% SDPD data
% =========================================================================
data_pack_SDPD.x = data_pack_SSA.x;

r=1;
for s = 1:nt
    for i=1:nx
        for j=1:ny
            data_pack_SDPD.field(i,j,s) = data_pack(i,j,1,r,s);
        end
    end
end


end
