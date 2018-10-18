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
markersize = 6;


% =========================================================================
% input parameters
% =========================================================================
nx = 33;                     % number of particles in x-direction
ny = 1;                      % number of particles in y-direction
C0 = 32000;                  % reference concentration
dt = 1e-2;                   % time step
freq_results = 100;          % freq. os saved results
Lx = 1;                      % total length
D = 1e-2;                    % diffusivity
t_steps = [1,2,4,8,16,99];   % read data in these specific time steps
n_realizations = 100;        % number of realizations of the SSA simulation
file_prefix = 'diffusion1d_results';




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
data_pack_exact = analytical_heat_1d(nx,D,C0,t_steps,dt,freq_results,data_pack_SDPD.x,Lx);



% =========================================================================
% compute RMS of erros
% =========================================================================




% =========================================================================
% plot results
% =========================================================================


figure(1);
hold on;
for s=1:length(t_steps)
      plt = errorbar(data_pack_SSA.x,data_pack_SSA.field_mean(:,:,s)/C0,data_pack_SSA.field_sd(:,:,s)/C0,'-ok','LineWidth',linewidth,'MarkerSize',markersize);
      plt = plot(data_pack_SDPD.x,data_pack_SDPD.field(:,:,s)/C0,'-bx','LineWidth',linewidth,'MarkerSize',markersize);
      plt = plot(data_pack_exact.x,data_pack_exact.field(:,s)/C0,'--r','LineWidth',linewidth,'MarkerSize',markersize);
end
% axis square;
set(gca,'FontSize',12)
xlabel('$x$','fontsize',18);
ylabel('$C/C_0$','fontsize',14);
% hold off;
lgn = legend('SSA','exact','SDPD');
% lgn = legend('t1','t2','t3','t4');
lgn.FontSize = 14;
axis([0 1 0 1]);
grid on;
hold off;
box on;
lgn.FontSize = 14;
tix=get(gca,'xtick')';
set(gca,'xticklabel',num2str(tix,'%.2f'))
tix=get(gca,'ytick')';
set(gca,'yticklabel',num2str(tix,'%.2f'))


% figure(2);
% contourf(Y,X,C_SDPD(:,:,end));
% axis square;
% xlabel('x');
% xlabel('y');
% 
% 
% figure(3);
% contourf(Y_exact,X_exact,C_exact);
% axis square;
% xlabel('x');
% xlabel('y');

% 
% % Vertical profile
% figure(5);
% plt = errorbar(data_pack_SSA.v_line_mean/C0,y,data_pack_SSA.v_line_sd/C0,'ko','LineWidth',linewidth,'MarkerSize',markersize);
% set(gca,'FontSize',12)
% hold on
% plt = plot(C_SDPD(ceil(nx/2),:,end)/C0,y,'-b','LineWidth',linewidth,'MarkerSize',markersize);
% plt = plot(data_pack_exact.v_line,data_pack_exact.y,'-r','LineWidth',linewidth,'MarkerSize',markersize);
% axis square;
% grid on;
% xlabel('$C/C_0$','fontsize',14);
% ylabel('$y$','fontsize',18);
% axis([0 1 0 1]);
% lgn = legend('SSA','SDPD','Exact','Location','northwest');
% lgn.FontSize = 14;

%figure('Units','normalized','Position',[0.05 0.02 0.94 0.87]);
% figure(5)
% plot(data_pack_statistics.v_line_mean/C0,y,'ko','linewidth', linewidth, 'MarkerSize', markersize)
% hold on
% plot(C_SDPD(ceil(nx/2),:,end)/C0,y,'-b','linewidth', linewidth, 'MarkerSize', markersize);
% xlabel('$C/C_0$','fontsize',20);
% ylabel('$y$','fontsize',20);
% grid on
% axis([0 1 0 1]);
% h = legend('SSA', 'SDPD', 'Exact', 'NorthEast');
% set(h, 'Interpreter', 'Latex')
% set(gca, 'fontsize', 10, 'xtick', 'ytick');
% hold off




% =========================================================================
% compute erros
% =========================================================================

% Compute average error (based on Frobenius norm) (SDPD)
s=6;
L2_SDPD = norm(data_pack_SDPD.field(:,:,s)/C0-data_pack_exact.field(:,s)/C0,'fro')/((nx-1));
fprintf("L2_SDPD = %.12f \n",L2_SDPD);

% Compute average error (based on Frobenius norm) (SSA)
L2_SSA = norm(data_pack_SSA.field_mean(:,:,s)/C0-data_pack_exact.field(:,s)/C0,'fro')/((nx-1));
fprintf("L2_SSA  = %.12f \n",L2_SSA);




end





