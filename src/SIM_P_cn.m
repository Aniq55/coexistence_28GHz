clear all
clc

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

%% Generate PPP for ES and BPs
L = 10e3;
Area = 4*L*L;

iter_count= 1e4;
SIR_vec_c = zeros(1, iter_count);

valid_iter = 0;

for iter = 1:iter_count
        
    % BS
    N_BS = poissrnd(lambda_c*Area);
    BS_x = -L + 2*L*rand(1, N_BS);
    BS_y = -L + 2*L*rand(1, N_BS);
    BS_r = abs(BS_x + 1i*BS_y);
    BS_fading = gamrnd(m, 1/m, 1, N_BS);
    [r_c, index] = min(BS_r);
    I_c_vec = p_c.*BS_fading.*max(1,BS_r).^(-zeta);
    S_c = I_c_vec(index);
    I_c = sum(I_c_vec);

    SIR_vec_c(iter) = S_c/(I_c - S_c);
end

%% Coverage
tau_vec_dB = [-60:1:60];
N_records = length(tau_vec_dB);
tau_vec_abs = 10.^(tau_vec_dB./10.0);
P_c_vec = zeros(1, length(tau_vec_abs));

for i=1:length(tau_vec_dB)
    tau = tau_vec_abs(i);
    P_c_vec(i) = sum(SIR_vec_c > tau)./length(SIR_vec_c);
end

csvwrite('D:\Satellites\28GHz\data\regression\P_cn.csv', P_c_vec);
