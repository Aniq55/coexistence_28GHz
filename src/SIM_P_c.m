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

h = h_vec(2);
N = N_vec(5);

%%
f_alpha_vec = real(csvread('D:\Satellites\28GHz\data\alpha0\f_alpha_'+string(int32(h.*1e-3)) ...
            +'_'+string(N)+'.csv'));
        
x_alpha_vec = real(csvread('D:\Satellites\28GHz\data\alpha0\x_alpha_'+string(int32(h.*1e-3)) ...
            +'_'+string(N)+'.csv'));
        
f_alpha_vec = f_alpha_vec./sum(f_alpha_vec); % PMF

F_alpha = cumsum(f_alpha_vec); %CDF

G_vec_read = csvread("D:\Satellites\28GHz\data\gain.csv"); % dB scale or dBm? [DOUBT]
G_vec = 10.^(G_vec_read./10.0); % absolute scale


%% Generate PPP for ES and BPs
L = 10e3;
Area = 4*L*L;

%% Exclusion Zone Radii
R_s = 500;
R_b = 500;

iter_count= 1e4;
SIR_vec_c = zeros(1, iter_count);

valid_iter = 0;

for iter = 1:iter_count

    % ES
    N_ES = poissrnd(lambda_s*Area);
    ES_x = -L + 2*L*rand(1, N_ES);
    ES_y = -L + 2*L*rand(1, N_ES);
    ES_r = abs(ES_x + 1i*ES_y);
    
    % BP
    N_BP = poissrnd(lambda_b*Area);
    BP_x = -L + 2*L*rand(1, N_BP);
    BP_y = -L + 2*L*rand(1, N_BP);
    BP_r = abs(BP_x + 1i*BP_y);
    
    
    min_ES_dist = min(ES_r);
    min_BP_dist = min(BP_r);
    
    if min_ES_dist>R_s && min_BP_dist>R_b
        
        valid_iter = valid_iter + 1;
        
        % ES
        ES_gain = zeros(1, N_ES);
        for i=1:N_ES
            [min_val, index] = min(abs(F_alpha - rand));
            ES_gain(i) = G_vec(index);
        end
        ES_fading = gamrnd(m, 1/m, 1, N_ES);
        I_s_vec = p_s.*ES_gain.*ES_fading.*max(1,ES_r).^(-zeta);
        I_s = sum(I_s_vec);
        
        % BP
        BP_gain_vec = [G0, G1];
        BP_gain = BP_gain_vec(1 + int32( 2*pi*rand(1, N_BP) < 2*Theta ));
        BP_fading = gamrnd(m, 1/m, 1, N_BP);
        I_b_vec = p_b.*BP_gain.*BP_fading.*max(1,BP_r).^(-zeta);
        I_b = sum(I_b_vec);
        
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
        
        SIR_vec_c(valid_iter) = S_c/(I_c - S_c + I_b + I_s);
        
    end
end


%% Coverage
tau_vec_dB = [-60:1:60];
N_records = length(tau_vec_dB);
tau_vec_abs = 10.^(tau_vec_dB./10.0);
P_c_vec = zeros(1, length(tau_vec_abs));
SIR_vec_c_valid = SIR_vec_c(1:valid_iter);

for i=1:length(tau_vec_dB)
    tau = tau_vec_abs(i);
    P_c_vec(i) = sum(SIR_vec_c_valid > tau)./length(SIR_vec_c_valid);
end


csvwrite( 'D:\Satellites\28GHz\data\regression\P_c_Rs_'...
            +string(int32(R_s))+'_Rb_'+string(int32(R_b))+'.csv', P_c_vec);
