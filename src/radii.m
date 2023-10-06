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

h = h_vec(1);
N = N_vec(5);

display(h+"  "+N)
% PDF [OK]
f_alpha_vec = real(csvread('D:\Satellites\28GHz\data\alpha0\f_alpha_'+string(int32(h.*1e-3)) ...
    +'_'+string(N)+'.csv'));

% Check: The PMF must sum to 1
d_alpha = 1/sum(f_alpha_vec); % 0.0326

G_vec_read = csvread("D:\Satellites\28GHz\data\gain.csv"); % dB scale or dBm? [DOUBT]
G_vec = 10.^(G_vec_read./10.0); % absolute scale

G_bar = sum(f_alpha_vec.*G_vec.*d_alpha);

chi_b = 2*pi*lambda_b*p_b*(G0 + (Theta/pi)*(G1-G0) )/(zeta-2);
chi_s = 2*pi*lambda_s*p_s*G_bar/(zeta-2);

I_dB = [-80: 1: -40];
I_abs = 10.0.^(I_dB./10.0);
n_I = length(I_dB);

R_b_vec = zeros(1, n_I);
R_s_vec = zeros(1, n_I);

syms R_s R_b mu real positive

for k=1:n_I    
    I_0 = I_abs(k);

    L1 = 2*pi*lambda_b*R_b + mu*(2-zeta)*chi_b*R_b^(1-zeta);
    L2 = 2*pi*lambda_s*R_s + mu*(2-zeta)*chi_s*R_s^(1-zeta);
    L3 = chi_b*R_b^(2-zeta) + chi_s*R_s^(2-zeta) - I_0;

    [a,b,c] = solve(L1==0, L2==0, L3==0);

    R_b_vec(k) = a;
    R_s_vec(k) = b;
end

%%

color_vec = ["219ebc", "fb8500" ];
c_1 = [hex2rgb(color_vec(1))];
c_2 = [hex2rgb(color_vec(2))];

figure;
semilogy(I_dB, R_b_vec, 'Color' , c_1, 'LineWidth', 1.5)
hold on;
semilogy(I_dB, R_s_vec, 'Color' , c_2, 'LineWidth', 1.5)
box on; grid on;
legend(['R_b'; 'R_s']);
ylabel('Exclusion Zone Radius [m]');
xlabel('Interference power [dB]');




