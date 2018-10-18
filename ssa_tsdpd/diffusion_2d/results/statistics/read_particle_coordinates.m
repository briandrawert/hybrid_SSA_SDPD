function [X,Y,x,y] = read_particle_coordinates(nx,ny,data_pack)

for i=1:nx
   x(i) = data_pack(i,3,end);
end

for j=1:ny
   k = (j-1)*ny + 1;
   y(j) = data_pack(k,4,end); 
end

[X,Y] = meshgrid(x,y);

end