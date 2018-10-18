function postprocessing

% =========================================================================
% initialize workspace
% =========================================================================
close all;

% =========================================================================
% plot options
% =========================================================================
set(0,'DefaultTextInterpreter','latex');
linewidth  = 1.2;
markersize = 5;


% =========================================================================
% input parameters
% =========================================================================
nx = 17;                           % number of particles in x-direction
ny = 1;                            % number of particles in y-direction
C0 = 16000;                        % reference concentration
dt = 1e-4;                         % time step
freq_results = 2000;               % freq. os saved results
Lx = 1;                            % total length
D = 0.1;                           % diffusivity
t_steps = [1 2 3 4 5];             % read data in these specific time steps
t_solution = [2000 4000 6000 8000 10000];
n_realizations = 100;              % number of realizations of the SSA simulation
file_prefix = 'cylinder_annihilation_results';



% =========================================================================
% read files
% =========================================================================
data_pack = read_files(n_realizations,t_steps,nx,ny,file_prefix);




% =========================================================================
% compute temporal statistics of SSA realizations over interval [t1,t2]
% =========================================================================
[data_pack_SSA,data_pack_SDPD] = temporal_statistics(t_steps,n_realizations,nx,ny,data_pack);




% =========================================================================
% compute steady-state exact solution
% =========================================================================
% data_pack_fd = fd_cylinder(nx,D,C0,t_solution,dt,freq_results,data_pack_SDPD.x);
data_pack_fd = compact_fd_cylinder(nx,D,C0,t_solution,dt,freq_results,data_pack_SDPD.x);



% =========================================================================
% compute RMS of erros
% =========================================================================




% =========================================================================
% plot results
% =========================================================================


figure(1);
hold on;
for s=1:length(t_steps)
      plt1 = errorbar(data_pack_SSA.x,data_pack_SSA.C0_mean(:,:,s)/C0,data_pack_SSA.C0_sd(:,:,s)/C0,'-ok','LineWidth',linewidth,'MarkerSize',markersize);
      plt2 = plot(data_pack_SDPD.x,data_pack_SDPD.C0(:,:,s)/C0,'-bx','LineWidth',linewidth,'MarkerSize',markersize);
      plt3 = plot(data_pack_fd.x,data_pack_fd.C0(:,s)/C0,'--k','LineWidth',linewidth,'MarkerSize',markersize);      
      plt4 = errorbar(data_pack_SSA.x,data_pack_SSA.C1_mean(:,:,s)/C0,data_pack_SSA.C1_sd(:,:,s)/C0,'-ok','LineWidth',linewidth,'MarkerSize',markersize);
      plt5 = plot(data_pack_SDPD.x,data_pack_SDPD.C1(:,:,s)/C0,'-rx','LineWidth',linewidth,'MarkerSize',markersize);
      plt6 = plot(data_pack_fd.x,data_pack_fd.C1(:,s)/C0,'--r','LineWidth',linewidth,'MarkerSize',markersize);      
end
% axis square;
set(gca,'FontSize',12)
xlabel('$x$','fontsize',18);
ylabel('$C/C_0$','fontsize',14);
hold off;
h = [plt1(1);plt2(1);plt3(1);plt4(1);plt5(1);plt6(1)];
lgn = legend(h,{'$C_A^{\rm{SSA}}$','$C_A^{\rm{SDPD}}$','$C_A^{\rm{FD}}$','$C_B^{\rm{SSA}}$','$C_B^{\rm{SDPD}}$','$C_B^{\rm{FD}}$'},'Interpreter','latex');

lgn.FontSize = 14;
% axis([0 -0.1 0 1.1]);
grid on;
hold off;
box on;
lgn.FontSize = 14;
tix=get(gca,'xtick')';
set(gca,'xticklabel',num2str(tix,'%.2f'))
tix=get(gca,'ytick')';
set(gca,'yticklabel',num2str(tix,'%.2f'))


% =========================================================================
% compute erros
% =========================================================================
% Concentration
% Compute average error (based on Frobenius norm) (SDPD)
L2_SDPD = norm(data_pack_SDPD.C0(:,:,end)/C0-data_pack_fd.C0(:,end)/C0,'fro')/((nx-1));
fprintf("L2_SDPD = %.16f \n",L2_SDPD);

% Compute average error (based on Frobenius norm) (SSA)
L2_SSA = norm(data_pack_SSA.C0_mean(:,:,end)/C0-data_pack_fd.C0(:,end)/C0,'fro')/((nx-1));
fprintf("L2_SSA  = %.16f \n",L2_SSA);


end





