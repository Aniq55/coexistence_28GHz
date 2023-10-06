clear all;
clc;

%% Parameters
p_s = 160;      % power satellite EIRP: 57.2 dBW, EIRP = 1.64*ERP, ERP = G*P_Tx
p_c = 16;       % power cellular [42 dBm]
p_b = 39.8;     % power backhaul [46 dBm]
zeta = 4;       % pathloss exponent

% Antenna pattern for backhaul point:
G1 = 1e2;
G0 = 1e-1;
Theta = pi/12;

vec_alpha = linspace(0, pi);
d_alpha = vec_alpha(2) - vec_alpha(1);

% Automation:
N_vec = [100, 500, 1000, 1500, 2000]; % number of satellites
h_vec = [500, 1000, 1500, 2000]*1e3; % orbit altitude

r_e= 6371e3;

lambda_c = 1e-6;
lambda_s = 1e-6;
lambda_b = 1e-7;

% Nakagami
m = 3;
mu = 1;
beta = (factorial(m))^(-1/m);

h = h_vec(2);
N = N_vec(5);

% B_cI = 1e9;     % 1 GHz
% B_cnI = 100e6;  % 100 MHz
% B_bI = 1e9;     % 1 GHz
% B_bnI = 100e6;  % 100 MHz

B = 20e6;
BW_factor = 5;

B_cI = BW_factor*B;     % 40 MHz
B_cnI = B;  % 20 MHz
B_bI = BW_factor*B;     % 40 MHz
B_bnI = B;  % 20 MHz

%% Variables

syms Rs Rb positive

tau = -10; %dB

%% load coefficients for slope and intercept

C_cI_s  = csvread("D:\Satellites\28GHz\data\fits\P_c_I_slope_fit.csv");
C_cI_i  = csvread("D:\Satellites\28GHz\data\fits\P_c_I_intercept_fit.csv");

C_cnI_s = csvread("D:\Satellites\28GHz\data\fits\P_c_notI_slope_fit.csv");
C_cnI_i = csvread("D:\Satellites\28GHz\data\fits\P_c_notI_intercept_fit.csv");

C_bI_s  = csvread("D:\Satellites\28GHz\data\fits\P_b_I_slope_fit.csv");
C_bI_i  = csvread("D:\Satellites\28GHz\data\fits\P_b_I_intercept_fit.csv");

C_bnI_s = csvread("D:\Satellites\28GHz\data\fits\P_b_notI_slope_fit.csv");
C_bnI_i = csvread("D:\Satellites\28GHz\data\fits\P_b_notI_intercept_fit.csv");

%% create slope and intercept functions

f_cI_s = @(Rs, Rb) C_cI_s(1).*Rs + C_cI_s(2).*Rs.*Rs ...
    + C_cI_s(3).*Rb + C_cI_s(4).*Rb.*Rb + C_cI_s(5);

f_cI_i = @(Rs, Rb) C_cI_i(1).*Rs + C_cI_i(2).*Rs.*Rs ...
    + C_cI_i(3).*Rb + C_cI_i(4).*Rb.*Rb + C_cI_i(5);

f_cnI_s = @(Rs, Rb) C_cnI_s(1).*Rs + C_cnI_s(2).*Rs.*Rs ...
    + C_cnI_s(3).*Rb + C_cnI_s(4).*Rb.*Rb + C_cnI_s(5);

f_cnI_i = @(Rs, Rb) C_cnI_i(1).*Rs + C_cnI_i(2).*Rs.*Rs ...
    + C_cnI_i(3).*Rb + C_cnI_i(4).*Rb.*Rb + C_cnI_i(5);

f_bI_s = @(Rs, Rb) C_bI_s(1).*Rs + C_bI_s(2).*Rs.*Rs ...
    + C_bI_s(3).*Rb + C_bI_s(4).*Rb.*Rb + C_bI_s(5);

f_bI_i = @(Rs, Rb) C_bI_i(1).*Rs + C_bI_i(2).*Rs.*Rs ...
    + C_bI_i(3).*Rb + C_bI_i(4).*Rb.*Rb + C_bI_i(5);

f_bnI_s = @(Rs, Rb) C_bnI_s(1).*Rs + C_bnI_s(2).*Rs.*Rs ...
    + C_bnI_s(3).*Rb + C_bnI_s(4).*Rb.*Rb + C_bnI_s(5);

f_bnI_i = @(Rs, Rb) C_bnI_i(1).*Rs + C_bnI_i(2).*Rs.*Rs ...
    + C_bnI_i(3).*Rb + C_bnI_i(4).*Rb.*Rb + C_bnI_i(5);

%% create the argument function for sigmoid

arg_cI      = @(Rs, Rb) tau.*f_cI_s(Rs, Rb)     + f_cI_i(Rs, Rb);
arg_cnI     = @(Rs, Rb) tau.*f_cnI_s(Rs, Rb)    + f_cnI_i(Rs, Rb);
arg_bI      = @(Rs, Rb) tau.*f_bI_s(Rs, Rb)     + f_bI_i(Rs, Rb);
arg_bnI     = @(Rs, Rb) tau.*f_bnI_s(Rs, Rb)    + f_bnI_i(Rs, Rb);

%% create the coverage functions

P_cI    = @(Rs, Rb) 1./(1 + exp(arg_cI(Rs, Rb)));
P_cnI   = @(Rs, Rb) 1./(1 + exp(arg_cnI(Rs, Rb)));
P_bI    = @(Rs, Rb) 1./(1 + exp(arg_bI(Rs, Rb)));
P_bnI   = @(Rs, Rb) 1./(1 + exp(arg_bnI(Rs, Rb)));

ex_s = @(Rs) exp(-pi.*lambda_s.*Rs.*Rs);
ex_b = @(Rb) exp(-pi.*lambda_b.*Rb.*Rb);

%% create the average data rate functions

D_c = @(Rs, Rb) ( ex_s(Rs).*ex_b(Rb).*(B_cI.*P_cI(Rs, Rb) ...
    - B_cnI.*P_cnI(Rs, Rb)) + B_cnI.*P_cnI(Rs, Rb) ).*log2(1 + 10.^(tau./10));

D_b = @(Rs, Rb) ( ex_s(Rs).*(B_bI.*P_bI(Rs, Rb) ...
    - B_bnI.*P_bnI(Rs, Rb)) + B_bnI.*P_bnI(Rs, Rb) ).*log2(1 + 10.^(tau./10));

% D_sum = @(Rs, Rb) D_c(Rs, Rb)+D_b(Rs, Rb);
D_prod = @(Rs, Rb) D_c(Rs, Rb).*D_b(Rs, Rb);

D_sum = @(Rs, Rb) D_c(Rs, Rb);


%%
R_s = [1:10:1e3];
R_b = [1:10:1e3];
n1 = int32(length(R_s));
n2 = int32(length(R_b));

F_Dc = zeros(n1,n2);
F_Db = zeros(n1,n2);
F_Dsum = zeros(n1,n2);
F_Dprod = zeros(n1,n2);

for i = 1:n1
    F_Dc(i,:) = feval(D_c, R_s(i), R_b);
    F_Db(i,:) = feval(D_b, R_s(i), R_b);
    F_Dsum(i,:) = feval(D_sum, R_s(i), R_b);
    F_Dprod(i,:) = feval(D_prod, R_s(i), R_b);
end

%%
% figure;
% hold on; grid on;
% s = surf(R_b, R_s, F_Dprod, 'FaceAlpha',0.5);
% s.EdgeColor = 'none';
% xlabel('R_b')
% ylabel('R_s')
%%
% 
% figure;
% contour(R_b, R_s, F_Dc)
% xlabel('R_b')
% ylabel('R_s')

%%
syms x
x0 = [100, 100];
f_Dc = @(x) -1.*D_c(x(1), x(2));
f_Db = @(x) -1.*D_b(x(1), x(2));
f_Dsum = @(x) -1.*D_sum(x(1), x(2));
f_Dprod = @(x) -1.*D_prod(x(1), x(2));

[x_Dc, fval_Dc] = fminsearchbnd(f_Dc, [100; 100], [0,0]);
[x_Db, fval_Db] = fminsearchbnd(f_Db, [100; 100], [0,0]);
[x_Dsum, fval_Dsum] = fminsearchbnd(f_Dsum, [100; 100], [0,0]);
[x_Dprod, fval_Dprod] = fminsearchbnd(f_Dprod, [100; 100], [0,0]);

%%
figure;
set(gcf,'Position',[100 100 1200 500])
% % tiledlayout(2,2, 'Padding', 'none', 'TileSpacing', 'compact'); 
% 
% subplot(2,2,1)
% % nexttile
% hold on; grid on; box on;
% contourf(R_b, R_s, F_Dc*1e-6, 10)
% colorbar
% scatter(x_Dc(2), x_Dc(1), 'xr')
% xlabel('R_b [m]')
% ylabel('R_s [m]')
% % xlim([0,1e5])
% % ylim([0,1e5])
% title('D_c')
% 
% 
% subplot(2,2,2)
% % nexttile
% hold on; grid on; box on;
% contourf(R_b, R_s, F_Db*1e-6, 10)
% colorbar
% scatter(x_Db(2), x_Db(1), 'xr')
% xlabel('R_b [m]')
% ylabel('R_s [m]')
% % xlim([0,1e5])
% % ylim([0,1e5])
% title('D_b')
% 
subplot(1,2,1)
% nexttile
hold on; grid on; box on;
contourf(R_b, R_s, F_Dsum, 10)
colorbar
plot(x_Dsum(2), x_Dsum(1), '*r', 'MarkerSize', 10)
xlabel('R_b [m]')
ylabel('R_s [m]')
% xlim([0,1e5])
% ylim([0,1e5])
title('f_\Sigma = D_c + D_b')
 
subplot(1,2,2)
% nexttile
hold on; grid on; box on;
contourf(R_b, R_s, F_Dprod, 10)
colorbar
plot(x_Dprod(2), x_Dprod(1), '*r', 'MarkerSize', 10)
xlabel('R_b [m]')
ylabel('R_s [m]')
% xlim([0,1e5])
% ylim([0,1e5])
title('f_\Pi = D_c \cdot D_b')

%%

D_c_zero =  feval(D_c, 0, 0);
D_b_zero =  feval(D_b, 0, 0);
Dzero_val = [D_c_zero, D_b_zero];

D_c_Dcmax = feval(D_c, x_Dc(1), x_Dc(2));
D_b_Dcmax = feval(D_b, x_Dc(1), x_Dc(2));
Dcmax_val = [D_c_Dcmax, D_b_Dcmax];

D_c_Dbmax = feval(D_c, x_Db(1), x_Db(2));
D_b_Dbmax = feval(D_b, x_Db(1), x_Db(2));
Dbmax_val = [D_c_Dbmax, D_b_Dbmax];

D_c_Dsummax = feval(D_c, x_Dsum(1), x_Dsum(2));
D_b_Dsummax = feval(D_b, x_Dsum(1), x_Dsum(2));
Dsummax_val = [D_c_Dsummax, D_b_Dsummax];

D_c_Dprodmax = feval(D_c, x_Dprod(1), x_Dprod(2));
D_b_Dprodmax = feval(D_b, x_Dprod(1), x_Dprod(2));
Dprodmax_val = [D_c_Dprodmax, D_b_Dprodmax];


%%

% figure;
% bar([Dzero_val; Dsummax_val; Dprodmax_val ]);
% box on; grid on;
% ylabel('Data rate [Mbps]')
% xlabel('Maximization Objective')
% xticks([1 2 3 4 5])
% xticklabels({'no exclusion zones', 'D_c+D_b', 'D_c \times D_b'})
% legend('cellular', 'backhaul')


%% 

BW_factor
sum_gains = 100*(Dsummax_val - Dzero_val)./Dzero_val
prod_gains = 100*(Dprodmax_val - Dzero_val)./Dzero_val

%%
figure;
plot(R_s, F_Dsum(1,:), 'LineWidth', 2)
xlabel('R_s [m]')
ylabel('D_c')
box on; grid on;

