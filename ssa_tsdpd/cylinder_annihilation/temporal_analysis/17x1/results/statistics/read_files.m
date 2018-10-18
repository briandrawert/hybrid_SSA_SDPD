function data_pack = read_files(n_realizations,t_steps,nx,ny,file_prefix)

n_steps   = length(t_steps);
data_pack = zeros(nx,ny,5,n_realizations,n_steps);
                  %x, y, cols, n_realizations, n_steps

for r = 1:n_realizations
    for s = 1:n_steps
        file = ['../../',num2str(r),'/',file_prefix,'.',num2str(t_steps(s)),'.csv'];
        aux=dlmread(file,',',1,0);
        %                       C_[0]     C_[1]     Cd_[0]    Cd_[1]   x
        data_pack(:,:,:,r,s) = [aux(:,2),aux(:,3),aux(:,4),aux(:,5),aux(:,12)];
    end
end

end