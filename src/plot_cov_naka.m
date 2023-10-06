clc;
clear all;

% Plotting
N_vec = [100, 500, 1000, 1500, 2000]; % number of satellites
h_vec = [500, 1000, 1500, 2000]*1e3; % orbit altitude

h = h_vec(2);
N = N_vec(5);

R_s = 100;
R_b = 100;

% SIMULATION:
Pc_sim = csvread('D:\Satellites\28GHz\data\regression\P_c_Rs_'...
    +string(int32(R_s))+'_Rb_'+string(int32(R_b))+'.csv');

Pb_sim = csvread('D:\Satellites\28GHz\data\regression\P_b_Rs_'...
    +string(int32(R_s))+'_Rb_'+string(int32(R_b))+'.csv');

Pcn_sim = csvread('D:\Satellites\28GHz\data\regression\P_cn.csv');

Pbn_sim = csvread('D:\Satellites\28GHz\data\regression\P_bn.csv');


% THEORY:
P_c = csvread('D:\Satellites\28GHz\data\coverage\P_c_'...
    +string(int32(h*1e-3))+'_'+string(N) ...
    +'_'+string(R_s) + '_'+string(R_b) +'.csv');
        
P_b = csvread('D:\Satellites\28GHz\data\coverage\P_b_'...
    +string(int32(h*1e-3))+'_'+string(N) ...
    +'_'+string(R_s) + '_'+string(R_b) +'.csv');

P_cn = csvread('D:\Satellites\28GHz\data\coverage\P_cn_'...
    +string(int32(h*1e-3))+'_'+string(N)+'.csv');

P_bn = csvread('D:\Satellites\28GHz\data\coverage\P_bn_'...
    +string(int32(h*1e-3))+'_'+string(N)+'.csv');



tau_vec_dB_sim = [-60:1:60];
tau_vec_dB = [-50:5:60];
%%

color_vec = ["219ebc", "023047", "ffb703", "fb8500" ];
c_1 = [hex2rgb(color_vec(1))];
c_2 = [hex2rgb(color_vec(2))];
c_3 = [hex2rgb(color_vec(3))];
c_4 = [hex2rgb(color_vec(4))];

figure(1);
hold on;
plot(tau_vec_dB_sim, Pc_sim, '-', 'Color' , c_1, 'LineWidth', 1.5)
plot(tau_vec_dB_sim, Pb_sim, '-', 'Color' , c_2, 'LineWidth', 1.5)
plot(tau_vec_dB_sim, Pcn_sim, '-', 'Color' , c_3, 'LineWidth', 1.5)
plot(tau_vec_dB_sim, Pbn_sim, '-', 'Color' , c_4, 'LineWidth', 1.5)

plot(tau_vec_dB, P_c, 's', 'Color' , c_1, 'LineWidth', 1.5)
plot(tau_vec_dB, P_b, 'x', 'Color' , c_2, 'LineWidth', 1.5) 
plot(tau_vec_dB, P_cn, 'v', 'Color' , c_3, 'LineWidth', 1.5)
plot(tau_vec_dB, P_bn, 'o', 'Color' , c_4, 'LineWidth', 1.5) 
xlabel('SIR [dB]')
ylabel('Coverage Probability')
ylim([0 1])
xlim([-60 60])
box on; grid on;
legend('P_c (simulation)', 'P_b (simulation)', 'P_{cn} (simulation)', 'P_{bn} (simulation)', ...
    'P_c (theory)', 'P_b (theory)', 'P_{cn} (theory)', 'P_{bn} (theory)');

%%
% 
% D_b_1 = P_b_1.*log2(1+10.^(tau_vec_dB./10.0));
% D_c_1 = P_c_1.*log2(1+10.^(tau_vec_dB./10.0));
% D_b_2 = P_b_2.*log2(1+10.^(tau_vec_dB./10.0));
% D_c_2 = P_c_2.*log2(1+10.^(tau_vec_dB./10.0));
% 
% figure(2);
% hold on;
% plot(tau_vec_dB, D_b_1, 's-', 'Color' , c_1, 'LineWidth', 1.5)
% plot(tau_vec_dB, D_b_2, '^-', 'Color' , c_1, 'LineWidth', 1.5)
% plot(tau_vec_dB, D_c_1, 's-', 'Color' , c_2, 'LineWidth', 1.5)
% plot(tau_vec_dB, D_c_2, '^-', 'Color' , c_2, 'LineWidth', 1.5) 
% xlabel('SIR [dB]')
% ylabel('Spectral Efficiency [bit/s/Hz]')
% xlim([-20 20])
% box on; grid on;
% legend('Backhaul point: h = '+ string(int32(h.*1e-3)) +' km, N = '+string(N), ...
%     'Backhaul point: h = '+ string(int32(h.*1e-3)) +' km, N = '+string(N),...
%     'Cellular user: h = '+ string(int32(h.*1e-3)) +' km, N = '+string(N), ...
%     'Cellular user: h = '+ string(int32(h.*1e-3)) +' km, N = '+string(N));
% 
% %% 
% 
% c_1 = [hex2rgb(color_vec(1))];
% c_2 = [hex2rgb(color_vec(2))];
% c_3 = [hex2rgb(color_vec(3))];
% c_4 = [hex2rgb(color_vec(4))];
% 
% figure(3);
% 
% subplot(121)
% hold on;
% plot(tau_vec_dB, P_b_1, 's-', 'Color' , c_1, 'LineWidth', 1.5)
% plot(tau_vec_dB, P_b_2, '^-', 'Color' , c_1, 'LineWidth', 1.5)
% plot(tau_vec_dB, P_c_1, 's-', 'Color' , c_2, 'LineWidth', 1.5)
% plot(tau_vec_dB, P_c_2, '^-', 'Color' , c_2, 'LineWidth', 1.5)
% xlabel('SIR [dB]')
% ylabel('Coverage Probability')
% xlim([-20 20])
% box on; grid on;
% 
% subplot(122)
% hold on;
% plot(tau_vec_dB, D_b_1, 's-', 'Color' , c_1, 'LineWidth', 1.5)
% plot(tau_vec_dB, D_b_2, '^-', 'Color' , c_1, 'LineWidth', 1.5)
% plot(tau_vec_dB, D_c_1, 's-', 'Color' , c_2, 'LineWidth', 1.5)
% plot(tau_vec_dB, D_c_2, '^-', 'Color' , c_2, 'LineWidth', 1.5) 
% xlabel('SIR [dB]')
% ylabel('Spectral Efficiency [bit/s/Hz]')
% xlim([-20 20])
% box on; grid on;
% 
% legend('Backhaul point: h = '+ string(int32(h.*1e-3)) +' km, N = '+string(N), ...
%     'Backhaul point: h = '+ string(int32(h.*1e-3)) +' km, N = '+string(N),...
%     'Cellular user: h = '+ string(int32(h.*1e-3)) +' km, N = '+string(N), ...
%     'Cellular user: h = '+ string(int32(h.*1e-3)) +' km, N = '+string(N));

%%

% Pc_0_sim = csvread('D:\Satellites\28GHz\data\regression\P_c_Rs_0_Rb_0.csv');
% Pc_100_sim = csvread('D:\Satellites\28GHz\data\regression\P_c_Rs_300_Rb_0.csv');
% Pc_500_500_sim = csvread('D:\Satellites\28GHz\data\regression\P_c_Rs_300_Rb_300.csv');
% 
% Pcn_sim = csvread('D:\Satellites\28GHz\data\regression\P_cn.csv');
% 
% figure(1);
% hold on;
% plot(tau_vec_dB_sim, Pc_0_sim, '-', 'Color' , c_1, 'LineWidth', 1.5)
% plot(tau_vec_dB_sim, Pc_100_sim, '-', 'Color' , c_2, 'LineWidth', 1.5)
% % plot(tau_vec_dB_sim, Pc_500_500_sim, '-', 'Color' , c_3, 'LineWidth', 1.5)
% 
% xlabel('SIR [dB]')
% ylabel('Coverage Probability of Cellular Users')
% ylim([0 1])
% xlim([-60 60])
% box on; grid on;
% legend('R_s= 0 m, R_b= 0 m', 'R_s= 100 m, R_b = 0 m');
