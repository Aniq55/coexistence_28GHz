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
G_avg = G1*Theta/pi + G0*(1- Theta/pi);

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
L = 50e3;
Area = 4*L*L;

iter_count= 1e4;
SIR_vec_b = zeros(1, iter_count);


for iter = 1:iter_count
  
    % BP
    N_BP = poissrnd(lambda_b*Area);
    BP_x = -L + 2*L*rand(1, N_BP);
    BP_y = -L + 2*L*rand(1, N_BP);
    BP_r = abs(BP_x + 1i*BP_y);
  
    % BP
    BP_gain_vec = [G0, G1];
    BP_gain = BP_gain_vec(1 + int32( rand(1, N_BP) < Theta/pi ));
    BP_fading = gamrnd(m, 1/m, 1, N_BP);
    [r_b, index] = min(BP_r);
    I_b_vec = p_b.*BP_gain.*BP_fading.*max(1,BP_r).^(-zeta);
    S_b = (I_b_vec(index)/BP_gain(index))*G1;
    I_b = sum(I_b_vec) - I_b_vec(index);

    SIR_vec_b(iter) = S_b/I_b;
end


%% Coverage
tau_vec_dB = [-60:1:60];
N_records = length(tau_vec_dB);
tau_vec_abs = 10.^(tau_vec_dB./10.0);
P_b_vec = zeros(1, length(tau_vec_abs));

for i=1:length(tau_vec_dB)
    tau = tau_vec_abs(i);
    P_b_vec(i) = sum(SIR_vec_b > tau)./length(SIR_vec_b);
end

csvwrite('D:\Satellites\28GHz\data\regression\P_bn.csv', P_b_vec);
