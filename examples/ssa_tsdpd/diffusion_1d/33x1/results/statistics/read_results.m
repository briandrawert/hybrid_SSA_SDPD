function [data_pack_SSA,data_pack_SDPD] = read_results(n_steps,n_realizations,nx,ny,data_pack)


%preallocate structures
field1  = 'field';
field2  = 'x';
value1{n_realizations}{n_steps}(nx,ny) = [];
value2{n_realizations}{n_steps}(nx,ny) = [];
data_pack_SDPD = struct(field1,value1,field2,value2);

field1  = 'field';
value1{n_realizations}{n_steps} = [];
data_pack_SSA  = struct(field1,value1);



for m = 1:n_realizations
    for n=1:n_steps
        for o=1:nx
            for p=1:ny
                data_pack_SDPD.field{m}{n}(o,p) = data_pack{m}{n}(o+(p-1)*ny,1);
                 data_pack_SSA.field{m}{n}(o,p) = data_pack{m}{n}(o+(p-1)*ny,2);

                data_pack_SDPD.x{m}{n}(o,p) = data_pack{m}{n}(o+(p-1)*ny,3);
                data_pack_SSA.x{m}{n}(o,p)  = data_pack{m}{n}(o+(p-1)*ny,3);

            end
        end
    end
end

end
