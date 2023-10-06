clc;
clear all;

%% PARAMETERS

lambda_c = 1e-6;
lambda_s = 1e-6;
lambda_b = 1e-7;

B = 20e6;
BW_factor = 5;

B_cI = BW_factor*B;     % 40 MHz
B_cnI = B;  % 20 MHz
B_bI = BW_factor*B;     % 40 MHz
B_bnI = B;  % 20 MHz

%% SIGMOID

% Variables

syms Rs Rb positive

tau = -30; %dB
tau_n = 31;

% load coefficients for slope and intercept

C_cI_s  = csvread("D:\Satellites\28GHz\data\fits\P_c_I_slope_fit.csv");
C_cI_i  = csvread("D:\Satellites\28GHz\data\fits\P_c_I_intercept_fit.csv");

C_cnI_s = csvread("D:\Satellites\28GHz\data\fits\P_c_notI_slope_fit.csv");
C_cnI_i = csvread("D:\Satellites\28GHz\data\fits\P_c_notI_intercept_fit.csv");

C_bI_s  = csvread("D:\Satellites\28GHz\data\fits\P_b_I_slope_fit.csv");
C_bI_i  = csvread("D:\Satellites\28GHz\data\fits\P_b_I_intercept_fit.csv");

C_bnI_s = csvread("D:\Satellites\28GHz\data\fits\P_b_notI_slope_fit.csv");
C_bnI_i = csvread("D:\Satellites\28GHz\data\fits\P_b_notI_intercept_fit.csv");

% create slope and intercept functions

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

% create the argument function for sigmoid

arg_cI      = @(Rs, Rb) tau.*f_cI_s(Rs, Rb)     + f_cI_i(Rs, Rb);
arg_cnI     = @(Rs, Rb) tau.*f_cnI_s(Rs, Rb)    + f_cnI_i(Rs, Rb);
arg_bI      = @(Rs, Rb) tau.*f_bI_s(Rs, Rb)     + f_bI_i(Rs, Rb);
arg_bnI     = @(Rs, Rb) tau.*f_bnI_s(Rs, Rb)    + f_bnI_i(Rs, Rb);

% create the coverage functions

P_cI    = @(Rs, Rb) 1./(1 + exp(arg_cI(Rs, Rb)));
P_cnI   = @(Rs, Rb) 1./(1 + exp(arg_cnI(Rs, Rb)));
P_bI    = @(Rs, Rb) 1./(1 + exp(arg_bI(Rs, Rb)));
P_bnI   = @(Rs, Rb) 1./(1 + exp(arg_bnI(Rs, Rb)));

ex_s = @(Rs) exp(-pi.*lambda_s.*Rs.*Rs);
ex_b = @(Rb) exp(-pi.*lambda_b.*Rb.*Rb);

% create the average data rate functions

D_c = @(Rs, Rb) ( ex_s(Rs).*ex_b(Rb).*(B_cI.*P_cI(Rs, Rb) ...
    - B_cnI.*P_cnI(Rs, Rb)) + B_cnI.*P_cnI(Rs, Rb) ).*log2(1 + 10.^(tau./10));

D_b = @(Rs, Rb) ( ex_s(Rs).*(B_bI.*P_bI(Rs, Rb) ...
    - B_bnI.*P_bnI(Rs, Rb)) + B_bnI.*P_bnI(Rs, Rb) ).*log2(1 + 10.^(tau./10));

% D_sum = @(Rs, Rb) D_c(Rs, Rb)+D_b(Rs, Rb);
D_prod = @(Rs, Rb) D_c(Rs, Rb).*D_b(Rs, Rb);

D_sum = @(Rs, Rb) D_c(Rs, Rb);

%% SIMULATION


tau_vec_dB = [-60:1:60];
Radius_list = [10: 10: 500];
n_radius = length(Radius_list);

%%

s_Pc = zeros(n_radius, n_radius);
s_Pb = zeros(n_radius, n_radius);
s_Pcn = zeros(n_radius, n_radius);
s_Pbn = zeros(n_radius, n_radius);
E_c = zeros(n_radius, n_radius);
E_b = zeros(n_radius, n_radius);
D_c_dot = zeros(n_radius, n_radius);
D_b_dot = zeros(n_radius, n_radius);

for i_Rs = 1:n_radius
    for i_Rb = 1:n_radius
    
    Rs = Radius_list(i_Rs);
    Rb = Radius_list(i_Rb);
    
    A = csvread( 'D:\Satellites\28GHz\data\regression\Rs_'+string(int32(Rs))+'_Rb_'+string(int32(Rb))+'.csv');
    a = A(tau_n,:); % tau = 0 dB <-> 61
    
    s_Pc(i_Rs, i_Rb) = a(3);
    s_Pb(i_Rs, i_Rb) = a(4);
    s_Pcn(i_Rs, i_Rb) = a(5);
    s_Pbn(i_Rs, i_Rb) = a(6);
    
    E_c(i_Rs, i_Rb) = exp(-pi*lambda_s*Rs^2 - pi*lambda_b*Rb^2);
    E_b(i_Rs, i_Rb) = exp(-pi*lambda_s*Rs^2);
    
    D_c_dot(i_Rs, i_Rb) = D_c(Rs, Rb);
    D_b_dot(i_Rs, i_Rb) = D_b(Rs, Rb);

    end
    
end

%%
% 
% D_c_sim = zeros(n_radius, n_radius);
% D_b_sim = zeros(n_radius, n_radius);

D_c_sim = (E_c.*(B_cI.*s_Pc - B_cnI.*s_Pcn) + B_cnI.*s_Pcn).*log2(1 + 10.^(tau./10));
D_b_sim = (E_b.*(B_bI.*s_Pb - B_bnI.*s_Pbn) + B_bnI.*s_Pbn).*log2(1 + 10.^(tau./10));



% %%
% figure;
% surf(Radius_list, Radius_list, D_c_sim);
% 
% figure;
% surf(Radius_list, Radius_list, D_b_sim);

%%
figure;
hold on; grid on; box on;
s1 = surf(Radius_list, Radius_list, D_c_sim*1e-6, 'FaceAlpha',0.5);
s2 = surf(Radius_list, Radius_list, D_c_dot*1e-6, 'FaceColor', 'k');
legend('simulation', 'sigmoid fit');
xlabel('R_s [m]')
ylabel('R_b [m]')
zlabel('D_c [Mbps]')

%%




