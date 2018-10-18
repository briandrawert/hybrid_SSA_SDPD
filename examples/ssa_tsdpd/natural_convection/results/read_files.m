function data_pack = read_files(cols,rows,file_name)

nx = length(rows);
ny = length(cols);

data_pack = zeros(nx,ny);


aux=dlmread(file_name,',',1,0);

for j=1:length(cols)   
    data_pack(:,j) = aux(rows(1):rows(end),cols(j));
end


end