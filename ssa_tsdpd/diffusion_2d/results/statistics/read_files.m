function data_pack = read_files(prefix,mesh,extension,nt)

for i=1:nt
    file = [prefix mesh '.' num2str(i-1) extension];
    aux=dlmread(file,',',1,0);
    data_pack(:,:,i) = [aux(:,2),aux(:,3),aux(:,8),aux(:,9)];
end

end