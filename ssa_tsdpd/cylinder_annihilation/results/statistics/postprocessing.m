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
C0 = 2601000;                % reference concentration
Lx = 1;                      % total length
D = 1e-2;                    % diffusivity


% =========================================================================
% read files
% =========================================================================

% Read SSA results
cols = 1:11;
rows = 58:94;
file_name = 'results_Ra10E4_SSA.csv';
data_pack_SSA_Ra10E4 = read_files(cols,rows,file_name);

file_name = 'results_Ra10E5_SSA.csv';
data_pack_SSA_Ra10E5 = read_files(cols,rows,file_name);

file_name = 'results_Ra10E6_SSA.csv';
data_pack_SSA_Ra10E6 = read_files(cols,rows,file_name);


% Read SDPD results
cols = 1:9;
rows = 59:94;
file_name = 'results_Ra10E4_SDPD.csv';
data_pack_SDPD_Ra10E4 = read_files(cols,rows,file_name);

file_name = 'results_Ra10E5_SDPD.csv';
data_pack_SDPD_Ra10E5 = read_files(cols,rows,file_name);

file_name = 'results_Ra10E6_SDPD.csv';
data_pack_SDPD_Ra10E6 = read_files(cols,rows,file_name);


% Read reference results
%Ra = 10^4
file_name = 'C_Ra_10E4_Moukalled_Acharya_1996.mat';
aux = load(file_name);
x = aux.C_Ra_10E4_Moukalled_Acharya_1996(:,1);
C = aux.C_Ra_10E4_Moukalled_Acharya_1996(:,2);
data_pack_reference_C_Ra10E4 = [x C];

file_name = 'vy_Ra_10E4_Moukalled_Acharya_1996.mat';
aux = load(file_name);
x = aux.vy_Ra_10E4_Moukalled_Acharya_1996(:,1);
vy = aux.vy_Ra_10E4_Moukalled_Acharya_1996(:,2);
data_pack_reference_v_Ra10E4 = [x vy];

%Ra = 10^5
file_name = 'C_Ra_10E5_Moukalled_Acharya_1996.mat';
aux = load(file_name);
x = aux.C_Ra_10E5_Moukalled_Acharya_1996(:,1);
C = aux.C_Ra_10E5_Moukalled_Acharya_1996(:,2);
data_pack_reference_C_Ra10E5 = [x C];

file_name = 'vy_Ra_10E5_Moukalled_Acharya_1996.mat';
aux = load(file_name);
x = aux.vy_Ra_10E5_Moukalled_Acharya_1996(:,1);
vy = aux.vy_Ra_10E5_Moukalled_Acharya_1996(:,2);
data_pack_reference_v_Ra10E5 = [x vy];

% %Ra = 10^6
file_name = 'C_Ra_10E6_Moukalled_Acharya_1996.mat';
aux = load(file_name);
x = aux.C_Ra_10E6_Moukalled_Acharya_1996(:,1);
C = aux.C_Ra_10E6_Moukalled_Acharya_1996(:,2);
data_pack_reference_C_Ra10E6 = [x C];

file_name = 'vy_Ra_10E6_Moukalled_Acharya_1996.mat';
aux = load(file_name);
x = aux.vy_Ra_10E6_Moukalled_Acharya_1996(:,1);
vy = aux.vy_Ra_10E6_Moukalled_Acharya_1996(:,2);
data_pack_reference_v_Ra10E6 = [x vy];



% =========================================================================
% plot results
% =========================================================================

% Concentration profile
figure(1);
hold on;
plt = errorbar(data_pack_SSA_Ra10E4(:,11)+0.5,data_pack_SSA_Ra10E4(:,3)/C0,data_pack_SSA_Ra10E4(:,4)/C0,'-ko','LineWidth',linewidth,'MarkerSize',markersize);
plt = plot(data_pack_SDPD_Ra10E4(:,9)+0.5,data_pack_SDPD_Ra10E4(:,1),'-bx','LineWidth',linewidth,'MarkerSize',markersize);
plt = plot(data_pack_reference_C_Ra10E4(1:4:end,1),data_pack_reference_C_Ra10E4(1:4:end,2),'rs','LineWidth',linewidth,'MarkerSize',markersize);


plt = errorbar(data_pack_SSA_Ra10E5(:,11)+0.5,data_pack_SSA_Ra10E5(:,3)/C0,data_pack_SSA_Ra10E5(:,4)/C0,'-ko','LineWidth',linewidth,'MarkerSize',markersize);
plt = plot(data_pack_SDPD_Ra10E5(:,9)+0.5,data_pack_SDPD_Ra10E5(:,1),'-bx','LineWidth',linewidth,'MarkerSize',markersize);
plt = plot(data_pack_reference_C_Ra10E5(1:4:end,1),data_pack_reference_C_Ra10E5(1:4:end,2),'rs','LineWidth',linewidth,'MarkerSize',markersize);


plt = errorbar(data_pack_SSA_Ra10E6(:,11)+0.5,data_pack_SSA_Ra10E6(:,3)/C0,data_pack_SSA_Ra10E6(:,4)/C0,'-ko','LineWidth',linewidth,'MarkerSize',markersize);
plt = plot(data_pack_SDPD_Ra10E6(:,9)+0.5,data_pack_SDPD_Ra10E6(:,1),'-bx','LineWidth',linewidth,'MarkerSize',markersize);
plt = plot(data_pack_reference_C_Ra10E6(1:4:end,1),data_pack_reference_C_Ra10E6(1:4:end,2),'rs','LineWidth',linewidth,'MarkerSize',markersize);


set(gca,'FontSize',12)
xlabel('$x^*$','fontsize',18);
ylabel('$C^*$','fontsize',14);
% hold off;
lgn = legend('SSA diffusion + SDPD advection','SDPD diffusion + SDPD advection','Moukalled and Acharya (1996)');
% lgn = legend('t1','t2','t3','t4');
lgn.FontSize = 14;
axis([0.6 1.0 0 1]);
grid on;
hold off;
box on;
lgn.FontSize = 12;
tix=get(gca,'xtick')';
set(gca,'xticklabel',num2str(tix,'%.2f'))
tix=get(gca,'ytick')';
set(gca,'yticklabel',num2str(tix,'%.2f'))



% Velocity profile
figure(2);
hold on;

fv = sqrt(10^4/0.70); %velocity's normalization factor
plt = errorbar(data_pack_SSA_Ra10E4(:,11)+0.5,data_pack_SSA_Ra10E4(:,6)*fv,data_pack_SSA_Ra10E4(:,9)*fv,'-ko','LineWidth',linewidth,'MarkerSize',markersize);
plt = plot(data_pack_SDPD_Ra10E4(:,9)+0.5,data_pack_SDPD_Ra10E4(:,4)*fv,'-bx','LineWidth',linewidth,'MarkerSize',markersize);
plt = plot(data_pack_reference_v_Ra10E4(1:2:end,1),-data_pack_reference_v_Ra10E4(1:2:end,2),'rs','LineWidth',linewidth,'MarkerSize',markersize);

fv = sqrt(10^5/0.70); %velocity's normalization factor
plt = errorbar(data_pack_SSA_Ra10E5(:,11)+0.5,data_pack_SSA_Ra10E5(:,6)*fv,data_pack_SSA_Ra10E5(:,9)*fv,'-ko','LineWidth',linewidth,'MarkerSize',markersize);
plt = plot(data_pack_SDPD_Ra10E5(:,9)+0.5,data_pack_SDPD_Ra10E5(:,4)*fv,'-bx','LineWidth',linewidth,'MarkerSize',markersize);
plt = plot(data_pack_reference_v_Ra10E5(1:2:end,1),-data_pack_reference_v_Ra10E5(1:2:end,2),'rs','LineWidth',linewidth,'MarkerSize',markersize);

fv = sqrt(10^6/0.70); %velocity's normalization factor
plt = errorbar(data_pack_SSA_Ra10E6(:,11)+0.5,data_pack_SSA_Ra10E6(:,6)*fv,data_pack_SSA_Ra10E6(:,9)*fv,'-ko','LineWidth',linewidth,'MarkerSize',markersize);
plt = plot(data_pack_SDPD_Ra10E6(:,9)+0.5,data_pack_SDPD_Ra10E6(:,4)*fv,'-bx','LineWidth',linewidth,'MarkerSize',markersize);
plt = plot(data_pack_reference_v_Ra10E6(1:2:end,1),-data_pack_reference_v_Ra10E6(1:2:end,2),'rs','LineWidth',linewidth,'MarkerSize',markersize);



set(gca,'FontSize',12)
xlabel('$x^*$','fontsize',18);
ylabel('$v_y^*$','fontsize',18);
lgn = legend('SSA diffusion + SDPD advection','SDPD diffusion + SDPD advection','Moukalled and Acharya (1996)');
lgn.FontSize = 14;
axis([0.58 1.0 -300 300]);
grid on;
hold off;
box on;
lgn.FontSize = 12;
tix=get(gca,'xtick')';
set(gca,'xticklabel',num2str(tix,'%.2f'))
tix=get(gca,'ytick')';
set(gca,'yticklabel',num2str(tix,'%.2f'))




% =========================================================================
% compute erros
% =========================================================================

% Compute average error (based on Frobenius norm) (SDPD)
% s=6;
% L2_SDPD = norm(data_pack_SDPD.field(:,:,s)/C0-data_pack_exact.field(:,s)/C0,'fro')/((nx-1));
% fprintf("L2_SDPD = %.12f \n",L2_SDPD);
% 
% % Compute average error (based on Frobenius norm) (SSA)
% L2_SSA = norm(data_pack_SSA.field_mean(:,:,s)/C0-data_pack_exact.field(:,s)/C0,'fro')/((nx-1));
% fprintf("L2_SSA  = %.12f \n",L2_SSA);




end
