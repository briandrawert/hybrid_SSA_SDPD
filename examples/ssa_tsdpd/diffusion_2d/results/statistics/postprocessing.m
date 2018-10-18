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
nx = 33;      % number of particles in x-direction
ny = 33;      % number of particles in y-direction
nt = 101;     % number of files to be read
C0 = 1089000; % reference concentration

% =========================================================================
% read files
% =========================================================================
prefix = 'diffusion2d_results_';
mesh = num2str(nx);
extension = '.csv';
data_pack = read_files(prefix,mesh,extension,nt);


% =========================================================================
% read particle coordinates and create mesh
% =========================================================================
[X,Y,x,y] = read_particle_coordinates(nx,ny,data_pack);


% =========================================================================
% read C_SSA and C_SDPD
% =========================================================================
[C_SSA,C_SDPD] = read_results(nx,ny,nt,data_pack);



% =========================================================================
% compute temporal statistics of SSA realizations over interval [t1,t2]
% =========================================================================
t1 = 50;
t2 = 100;
data_pack_statistics = temporal_statistics(t1,t2,nx,ny,C_SSA);


% =========================================================================
% compute steady-state exact solution
% =========================================================================
data_pack_exact = analytical_heat_2d(nx,ny);




% =========================================================================
% plot results
% =========================================================================

% figure(1);
% [C,h] = contour(Y,X,data_pack_statistics.field_mean/C0,':');
% h.LineColor = 'k';
% h.LineWidth = 2;
% hold on
% [C,h] = contour(data_pack_exact.Y,data_pack_exact.X,data_pack_exact.field);;
% axis square;
% set(gca,'FontSize',12)
% h.LineColor = 'k';
% h.LineWidth = 1.0;
% clabel(C,'manual','FontSize',11,'Color','k','Interpreter','latex')
% xlabel('$x$','FontSize',18);
% ylabel('$y$','FontSize',18);
% lgn = legend('SSA','Exact','Location','southwest');
% lgn.FontSize = 14;
% grid on


% colorbar;

% figure(2);
% contour(Y,X,C_SDPD(:,:,end)/C0);
% axis square;
% xlabel('x');
% xlabel('y');
% colorbar;


% figure(3);
% contour(data_pack_exact.Y,data_pack_exact.X,data_pack_exact.field);
% axis square;
% xlabel('x');
% xlabel('y');
% colorbar;


% Vertical profile
% subplot(2,1,1)
figure(5);
plt = errorbar(data_pack_statistics.v_line_mean/C0,y,data_pack_statistics.v_line_sd/C0,'k-o','LineWidth',linewidth,'MarkerSize',markersize);
set(gca,'FontSize',12)
hold on
plt = plot(C_SDPD(ceil(nx/2),:,end)/C0,y,'-bx','LineWidth',linewidth,'MarkerSize',markersize);
plt = plot(data_pack_exact.v_line,data_pack_exact.y,'--r','LineWidth',linewidth,'MarkerSize',markersize);
% axis square;
grid on;
xlabel('$C/C_0$','fontsize',14);
ylabel('$y$','fontsize',18);
axis([0 1 0 1]);
lgn = legend('SSA','SDPD','Exact','Location','northwest');
lgn.FontSize = 14;
tix=get(gca,'xtick')';
set(gca,'xticklabel',num2str(tix,'%.2f'))
tix=get(gca,'ytick')';
set(gca,'yticklabel',num2str(tix,'%.2f'))



% Horizontal profile
figure(6);
plt = errorbar(x,data_pack_statistics.h_line_mean/C0,data_pack_statistics.h_line_sd/C0,'k-o','LineWidth',linewidth,'MarkerSize',markersize);
set(gca,'FontSize',12)
hold on
plt = plot(x,C_SDPD(:,ceil(ny/2),end)/C0,'-bx','LineWidth',linewidth,'MarkerSize',markersize);
plt = plot(data_pack_exact.x,data_pack_exact.h_line,'--r','LineWidth',linewidth,'MarkerSize',markersize);
% axis square;
grid on;
ylabel('$C/C_0$','fontsize',14);
xlabel('$x$','fontsize',18);
% axis([0 1 0 0.3]);
lgn = legend('SSA','SDPD','Exact','Location','northwest');
lgn.FontSize = 14;
tix=get(gca,'xtick')';
set(gca,'xticklabel',num2str(tix,'%.2f'))
tix=get(gca,'ytick')';
set(gca,'yticklabel',num2str(tix,'%.2f'))






% =========================================================================
% compute erros
% =========================================================================

% Compute average error (based on Frobenius norm) (SDPD)
L2_SDPD = norm(C_SDPD(2:nx-1,2:ny-1,end)/C0-data_pack_exact.field(2:nx-1,2:ny-1),'fro')/((nx-1)*(ny-1));
fprintf("L2_SDPD = %.12f \n",L2_SDPD);

% Compute average error (based on Frobenius norm) (SSA)
L2_SSA = norm(data_pack_statistics.field_mean(2:nx-1,2:ny-1,end)/C0-data_pack_exact.field(2:nx-1,2:ny-1),'fro')/((nx-1)*(ny-1));
fprintf("L2_SSA = %.12f \n",L2_SSA);





end


