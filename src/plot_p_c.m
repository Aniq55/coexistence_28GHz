tau_vec_dB = [-20:5:20];
% P_c_vec = zeros(4, length(tau_vec_dB));

% % Automation:
% N_vec = [100, 500, 1000, 1500, 2000]; % number of satellites
% h_vec = [500, 1000, 1500, 2000]*1e3; % orbit altitude
% 
% h = h_vec(1);
% for j = 1:5
%     N = N_vec(j);
%     P_c_vec(j,:) = real(csvread('D:\Satellites\28GHz\data\coverage\P_c_'+string(int32(h.*1e-3)) ...
%     +'_'+string(N)+'.csv'));
% end
% 
% %%
% figure('Position', [20 20 400 300])
% plot(tau_vec_dB, P_c_vec)
% legend(string(N_vec))

%%

P_c1 = real(csvread('D:\Satellites\28GHz\data\coverage\P_c_1000_100.csv'));
P_c2 = real(csvread('D:\Satellites\28GHz\data\coverage\P_c_1000_2000.csv'));


color_vec = ["219ebc", "fb8500" ];
c_1 = [hex2rgb(color_vec(1))];
c_2 = [hex2rgb(color_vec(2))];

figure;
plot([-20:5:20], P_c1, 'Color' , c_1, 'LineWidth', 1.5);
hold on;
plot([-20:5:5], P_c2, 'Color' , c_2, 'LineWidth', 1.5);
% xlim([-10 30])
lg = legend(["100", "2000"])
title(lg, "N")
grid on
box on
ylabel("Coverage Probability")
xlabel("SIR [dB]")






