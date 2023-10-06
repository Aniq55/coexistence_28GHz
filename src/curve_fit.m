clear all;
clc;

slopes = readtable("D:\Satellites\28GHz\data\fits\P_c_I_slopes.csv")

R_s = slopes.R_s;
R_b = slopes.R_b;
m = length(R_s);

R = [10:50:500];
n = length(R);

C_1 = csvread("D:\Satellites\28GHz\data\fits\P_c_I_slope_fit.csv");
C_2 = csvread("D:\Satellites\28GHz\data\fits\P_c_I_intercept_fit.csv");

X = [R_s, R_s.^2, R_b, R_b.^2, ones(m,1)];
%%
f_1 = X*C_1;
f_2 = X*C_2;

f_1_mat = reshape(f_1, n, n);
f_2_mat = reshape(f_2, n, n);

%% slope
figure;
hold on; grid on;
s_1 = surf(R, R, f_1_mat, 'FaceAlpha',0.5);
s_1.EdgeColor = 'none';
scatter3(slopes.R_s, slopes.R_b, slopes.slope, 's')
xlabel('R_s')
ylabel('R_b')
zlabel("weight")
legend('$f_w(R_s, R_b, \mathbf{a}_1)$', '$w_1(R_s, R_b)$', 'Interpreter', 'latex')

%% intercept
figure;
hold on; grid on;
s_2 = surf(R, R, f_2_mat, 'FaceAlpha',0.5);
s_2.EdgeColor = 'none';
scatter3(slopes.R_s, slopes.R_b, slopes.intercept, 's')
xlabel('R_s')
ylabel('R_b')
zlabel("bias")
legend('$f_w(R_s, R_b, \mathbf{a}_0)$', '$w_0(R_s, R_b)$', 'Interpreter', 'latex')


%%

F = 1./(1 + exp(f_1*(-10) + f_2)); % tau in dB
F_mat = reshape(F, n, n);

figure;
hold on; grid on;
s_2 = surf(R, R, F_mat, 'FaceAlpha',0.5);
s_2.EdgeColor = 'none';
xlabel('R_s [m]')
ylabel('R_b [m]')
zlabel('cellular coverage probability')


%% Generate coverage plots
tau_vec_dB = [-60:1:60];
t_len = length(tau_vec_dB);
R_list = [10:50:500];

Rs = 110;
Rb = 10;

s_ = slopes.slope([slopes.R_s == Rs & slopes.R_b == Rb]);
i_ = slopes.intercept([slopes.R_s == Rs & slopes.R_b == Rb]);

X_vec = [Rs, Rs.^2, Rb, Rb.^2, 1.0];

P_c_actual = csvread( 'D:\Satellites\28GHz\data\regression\Rs_'+string(int32(Rs))+'_Rb_'+string(int32(Rb))+'.csv');
P_c = zeros(1, t_len);
P_c_fit = zeros(1, t_len);
for i = 1:t_len
    tau = tau_vec_dB(i);
    P_c(i) = 1/(1 + exp(X_vec*(C_1*tau + C_2)));
    P_c_fit(i) = 1/(1 + exp(s_*tau + i_));
end


figure;
hold on; grid on; box on;
plot(tau_vec_dB, P_c_actual(:,3), '-', 'LineWidth', 4)
plot(tau_vec_dB, P_c_fit, 'g-', 'LineWidth', 1.5)
plot(tau_vec_dB, P_c, 'k--', 'LineWidth', 1.5)
title("R_s =" + string(Rs)+", R_b = "+string(Rb))
xlabel('\tau [dB]')
ylabel('Coverage Probability')
legend('P_c (simulated)', 'p_1 (sigmoid-I fit)', 'p_2 (sigmoid-II fit)')

%%
sum((P_c' - P_c_actual(:,6)).^2)./length(P_c)
sum((P_c_fit' - P_c_actual(:,6)).^2)./length(P_c)
sum((P_c - P_c_fit).^2)./length(P_c)


%% Automate Error in Fit:
E = zeros(n*n, 7);
k = 1;
for Rs = R_list
    for Rb = R_list
        s_ = slopes.slope([slopes.R_s == Rs & slopes.R_b == Rb]);
        i_ = slopes.intercept([slopes.R_s == Rs & slopes.R_b == Rb]);

        X_vec = [Rs, Rs.^2, Rb, Rb.^2, 1.0];

        P_c_actual_read = csvread( 'D:\Satellites\28GHz\data\regression\Rs_'+ ...
            string(int32(Rs))+'_Rb_'+string(int32(Rb))+'.csv');
        P_c_actual = P_c_actual_read(:,3)';
        P_c = zeros(1, t_len);
        P_c_fit = zeros(1, t_len);
        for i = 1:t_len
            tau = tau_vec_dB(i);
            P_c(i) = 1/(1 + exp(X_vec*(C_1*tau +C_2)));
            P_c_fit(i) = 1/(1 + exp(s_*tau + i_));
        end
        [h,p] = kstest2(P_c, P_c_actual);
        E(k,:) = [Rs, Rb, ...
            sqrt(sum((P_c - P_c_actual).^2)./length(P_c)), ...
            sqrt(sum((P_c_fit - P_c_actual).^2)./length(P_c)), ...
            sqrt(sum((P_c - P_c_fit).^2)./length(P_c)), ...
            h, p];
        
        k= k + 1;
    end
end


%%
figure;
hold on; box on; grid on;
plot(E(:,1), E(:,3), 'x');
plot(E(:,1), E(:,4), '.');
plot(E(:,1), E(:,5), 's');
xlabel("R_s [m]");
ylabel("RMSE");
legend('p_2 vs. P_c', ...
    'p_1 vs. P_c', ...
    'p_1 vs. p_2' );

%%
figure;
hold on; box on; grid on;
plot(E(:,2), E(:,3), 'x');
plot(E(:,2), E(:,4), '.');
plot(E(:,2), E(:,5), 's');
xlabel("R_b [m]");
ylabel("RMSE");
legend('p_2 vs. P_c', ...
    'p_1 vs. P_c', ...
    'p_1 vs. p_2' );

%% Histogram of RMSE

figure;
hold on; box on;
histogram(E(:,3), 10)
histogram(E(:,4), 10)
histogram(E(:,5), 10)
xlabel("RMSE");
ylabel("count");
legend('p_2 vs. P_c', ...
    'p_1 vs. P_c', ...
    'p_1 vs. p_2' );

figure;
hold on; box on; grid on;
ecdf(E(:,3))
ecdf(E(:,4))
ecdf(E(:,5))
xlabel("RMSE");
ylabel("CDF");
legend('p_2 vs. P_c', ...
    'p_1 vs. P_c', ...
    'p_1 vs. p_2' );


%% KS test
figure;
scatter3(E(:,1), E(:,2), E(:,7))
xlabel('R_s')
ylabel('R_b')
zlabel('KS-Test')

%%
data = F_mat;
[X,Y] = meshgrid(1:size(data,2), 1:size(data,1));
%// Define a finer grid of points
[X2,Y2] = meshgrid(1:0.01:size(data,2), 1:0.01:size(data,1));
%// Interpolate the data and show the output
outData = interp2(X, Y, data, X2, Y2, 'linear');

figure;
imagesc(outData);
colorbar
xlabel('R_s')
ylabel('R_b')
zlabel('KS-Test')

%%
figure;
[M,c] = contour(R_list, R_list, data);
c.LineWidth = 2;
xlabel('R_s [m]');
ylabel('R_b [m]');
title("Cellular Coverage Probability");

%%
scatter3(E(:,1), E(:,2), E(:,4), ...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0 .75 .75])
