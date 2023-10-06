clear all
clc

%% Parameters
p_s = 160;     % power satellite EIRP: 57.2 dBW, EIRP = 1.64*ERP, ERP = G*P_Tx
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

R_s = 100;
R_b = 100;

% for h = h_vec
%     for N  = N_vec
        display(h+"  "+N)
        % PDF [OK]
        f_alpha_vec = real(csvread('D:\Satellites\28GHz\data\alpha0\f_alpha_'+string(int32(h.*1e-3)) ...
            +'_'+string(N)+'.csv'));

        % Check: The PMF must sum to 1
        d_alpha = 1/sum(f_alpha_vec); % 0.0326

        G_vec_read = csvread("D:\Satellites\28GHz\data\gain.csv"); % dB scale or dBm? [DOUBT]
        G_vec = 10.^(G_vec_read./10.0); % absolute scale
        
        
        %% ALTERNATE DEFINITIONS:
        lambda_c_tilde = lambda_c*exp( -(pi*lambda_c*(R_s - 0.5/sqrt(lambda_c))^2)*...
            double(R_s > 0.5/sqrt(lambda_c)) ...
            -(pi*lambda_b*(R_b - 0.5/sqrt(lambda_c))^2)*...
            double(R_b > 0.5/sqrt(lambda_c)));

%         lambda_c_tilde = lambda_c*exp( -pi*lambda_c*R_s^2 -pi*lambda_b*R_b^2 );
%         lambda_b_tilde = lambda_b*exp( -pi*lambda_b*R_b^2 );

%         lambda_c_tilde = lambda_c;
%         lambda_b_tilde = lambda_b;

        syms tau v x G z

        

        %% Interference from Cellular Base Stations
        S_2_1 =  (1 - (1 + (1/mu).*beta.*tau.*(v./x).^zeta ).^-m ).*x; % [OK]
        S_2_2 =  (1 - (1 + (2/mu).*beta.*tau.*(v./x).^zeta ).^-m ).*x; % [OK]
        S_2_3 =  (1 - (1 + (3/mu).*beta.*tau.*(v./x).^zeta ).^-m ).*x; % [OK]
        
        I_3_1 = vpaintegral(S_2_1, x, [v Inf]); % [OK]
        I_3_2 = vpaintegral(S_2_2, x, [v Inf]); % [OK]
        I_3_3 = vpaintegral(S_2_3, x, [v Inf]); % [OK]
        
        L_Ic_1 = exp(-2.*pi.*lambda_c.*I_3_1); % [OK]
        L_Ic_2 = exp(-2.*pi.*lambda_c.*I_3_2); % [OK]
        L_Ic_3 = exp(-2.*pi.*lambda_c.*I_3_3); % [OK]

        %% Coverage probability:

        % Evaluate P_c
        % Check: P_c(0) should be 1 [OK]

        v_val_vec = [1.0: 1.0: 1e6];
        F_v = @(z) exp(-pi.*lambda_c.*z.^2);
        F_v_vec = F_v(v_val_vec);

        iter = 1e3;

        tau_vec_dB = [-50:5:60];
        tau_vec_abs = 10.^(tau_vec_dB./10.0);
        P_c_vec = zeros(1, length(tau_vec_abs));

        for j = 1: length(tau_vec_abs)
            tau_val = tau_vec_abs(j);
            
            S_vec_1 = zeros(1, iter);
            S_vec_2 = zeros(1, iter);
            S_vec_3 = zeros(1, iter);
            
            tic
            for i= 1:iter
                u = rand(1);
                [min_val, index] = min(abs(F_v_vec - u));
                try
                    L_Ic_1_mc = vpa(subs(L_Ic_1, {tau, v}, {tau_val, v_val_vec(index)}));
                    L_Ic_2_mc = vpa(subs(L_Ic_2, {tau, v}, {tau_val, v_val_vec(index)}));
                    L_Ic_3_mc = vpa(subs(L_Ic_3, {tau, v}, {tau_val, v_val_vec(index)}));

                    S_vec_1(i) = L_Ic_1_mc;
                    S_vec_2(i) = L_Ic_2_mc;
                    S_vec_3(i) = L_Ic_3_mc;
                    
                catch
                    S_vec_1(i) = NaN;
                    S_vec_2(i) = NaN;
                    S_vec_3(i) = NaN;
                end
            end
            toc

            p1 = mean(S_vec_1(~isnan(S_vec_1)));
            p2 = mean(S_vec_2(~isnan(S_vec_2)));
            p3 = mean(S_vec_3(~isnan(S_vec_3)));
            
            P_c_vec(j) = 3*p1 - 3*p2 + p3;
        end

        %% Save P_c_vec
        csvwrite('D:\Satellites\28GHz\data\coverage\P_cn_'...
            +string(int32(h*1e-3))+'_'+string(N) +'.csv', P_c_vec);

%     end
% end


% %% Figure
% figure;
% plot(tau_vec_dB, P_c_vec, 's-', 'LineWidth', 1.5)
% xlabel('SIR [dB]')
% ylabel('Coverage Probability')
% title('h = '+ string(int32(h.*1e-3)) +' km, N = '+string(N))
% ylim([0 1])
% xlim([-20 20])

