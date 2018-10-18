


function [data_pack_statistics] = temporal_statistics(t1,t2,nx,ny,Y)

nt = t2 - t1 + 1;


field1  = 'field_mean';        value1  = zeros(nx,ny);
field2  = 'field_prime';       value2  = zeros(nx,ny,nt);
field3  = 'field_variance';    value3  = zeros(1,nt);
field4  = 'field_sd';          value4  = zeros(1,nt);

field5  = 'h_line_mean';       value5  = zeros(nx,1);
field6  = 'v_line_mean';       value6  = zeros(1,ny);
field7  = 'h_line_prime';      value7  = zeros(nx,1,nt);
field8  = 'v_line_prime';      value8  = zeros(1,ny,nt);
field9  = 'h_line_variance';   value9  = zeros(nx,1);
field10 = 'v_line_variance';   value10 = zeros(1,ny);
field11 = 'h_line_sd';         value11 = zeros(nx,1);
field12 = 'v_line_sd';         value12 = zeros(1,ny);




data_pack_statistics = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6, ...
      field7,value7,field8,value8,field9,value9,field10,value10,field11,value11,field12,value12);


% =========================================================================
% stats over the entire field
% =========================================================================

% compute mean (1st order moment)
for t = t1:t2
    data_pack_statistics.field_mean = data_pack_statistics.field_mean + Y(:,:,t)/nt;
end

% compute fluctuations
k=1;
for t = t1:t2
    data_pack_statistics.field_prime(:,:,k) = Y(:,:,t) - data_pack_statistics.field_mean;
    k=k+1;
end

% compute variance (2nd order moment)
k=1;
for t = t1:t2
    for i=1:nx
        for j=1:ny
            data_pack_statistics.field_variance(k) = data_pack_statistics.field_variance(k) + (data_pack_statistics.field_prime(i,j,k)^2)/(nx*ny);
        end
    end
    k=k+1;
end

% compute standard deviation
data_pack_statistics.field_sd = sqrt(data_pack_statistics.field_variance);



% =========================================================================
% stats over centerlines
% =========================================================================
%Extract centerlines
hor_centerline = Y(:,ceil(ny/2),:);
ver_centerline = Y(ceil(nx/2),:,:);

    % =====================================================================
    % Horizontal centerline
    % =====================================================================
    for t = t1:t2
        data_pack_statistics.h_line_mean = data_pack_statistics.h_line_mean + hor_centerline(:,:,t)/nt;
    end

    % compute fluctuations
    k=1;
    for t = t1:t2
        data_pack_statistics.h_line_prime(:,:,k) = hor_centerline(:,:,t) - data_pack_statistics.h_line_mean;
        k=k+1;
    end
    
    % compute variance
    for i=1:nx
        k=1;
        for t = t1:t2
            data_pack_statistics.h_line_variance(i) = data_pack_statistics.h_line_variance(i) + (data_pack_statistics.h_line_prime(i,:,k)^2)/(nt);
            k=k+1;
        end
    end
    
    %compute standard deviation
    data_pack_statistics.h_line_sd = sqrt(data_pack_statistics.h_line_variance);
    
    
    % =====================================================================
    % Vertical centerline
    % =====================================================================
    for t = t1:t2
        data_pack_statistics.v_line_mean = data_pack_statistics.v_line_mean + ver_centerline(:,:,t)/nt;
    end

    % compute fluctuations
    k=1;
    for t = t1:t2
        data_pack_statistics.v_line_prime(:,:,k) = ver_centerline(:,:,t) - data_pack_statistics.v_line_mean;
        k=k+1;
    end
    
    % compute variance
    
    for j=1:ny
        k=1;
        for t = t1:t2
            data_pack_statistics.v_line_variance(j) = data_pack_statistics.v_line_variance(j) + (data_pack_statistics.v_line_prime(:,j,k)^2)/(nt);
            k=k+1;
        end
    end
    
    %compute standard deviation
    data_pack_statistics.v_line_sd = sqrt(data_pack_statistics.v_line_variance);


end
