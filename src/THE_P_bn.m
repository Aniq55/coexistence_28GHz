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
        

        syms tau v x G
        %% Interference from Backhaul Points
        S_b_1_g1 = (1 + 1.*G1.*(p_b./(p_b.*G1)).*beta.*tau.*(v./x).^zeta ).^-m;
        S_b_1_g0 = (1 + 1.*G0.*(p_b./(p_b.*G1)).*beta.*tau.*(v./x).^zeta ).^-m;
        
        S_b_2_g1 = (1 + 2.*G1.*(p_b./(p_b.*G1)).*beta.*tau.*(v./x).^zeta ).^-m;
        S_b_2_g0 = (1 + 2.*G0.*(p_b./(p_b.*G1)).*beta.*tau.*(v./x).^zeta ).^-m;
        
        S_b_3_g1 = (1 + 3.*G1.*(p_b./(p_b.*G1)).*beta.*tau.*(v./x).^zeta ).^-m;
        S_b_3_g0 = (1 + 3.*G0.*(p_b./(p_b.*G1)).*beta.*tau.*(v./x).^zeta ).^-m;
        
        
        I_b_1_g1 =  vpaintegral( (1 - S_b_1_g1).*x, x, [v Inf]); % [OK]
        I_b_1_g0 =  vpaintegral( (1 - S_b_1_g0).*x, x, [v Inf]); % [OK]
        
        I_b_2_g1 =  vpaintegral( (1 - S_b_2_g1).*x, x, [v Inf]); % [OK]
        I_b_2_g0 =  vpaintegral( (1 - S_b_2_g0).*x, x, [v Inf]); % [OK]
        
        I_b_3_g1 =  vpaintegral( (1 - S_b_3_g1).*x, x, [v Inf]); % [OK]
        I_b_3_g0 =  vpaintegral( (1 - S_b_3_g0).*x, x, [v Inf]); % [OK]
        
        L_Ib_1 = exp(-2.*pi.*(Theta/pi).*lambda_b.*I_b_1_g1).* exp(-2.*pi.*(1 - Theta/pi).*lambda_b.*I_b_1_g0); % [OK]
        L_Ib_2 = exp(-2.*pi.*(Theta/pi).*lambda_b.*I_b_2_g1).* exp(-2.*pi.*(1 - Theta/pi).*lambda_b.*I_b_2_g0); % [OK]
        L_Ib_3 = exp(-2.*pi.*(Theta/pi).*lambda_b.*I_b_3_g1).* exp(-2.*pi.*(1 - Theta/pi).*lambda_b.*I_b_3_g0); % [OK]

        %% Coverage probability:

        % Evaluate P_b
        % Check: P_b(0) should be 1 [OK]

        v_val_vec = [1.0: 1.0: 1e6];
        F_v = @(z) exp(-pi.*lambda_b.*z.^2);
        F_v_vec = F_v(v_val_vec);

        iter = 1e3;

        tau_vec_dB = [-50:5:60];
        tau_vec_abs = 10.^(tau_vec_dB./10.0);
        P_b_vec = zeros(1, length(tau_vec_abs));

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
                    L_Ib_1_mc = vpa(subs(L_Ib_1, {tau, v}, {tau_val, v_val_vec(index)}));
                    L_Ib_2_mc = vpa(subs(L_Ib_2, {tau, v}, {tau_val, v_val_vec(index)}));
                    L_Ib_3_mc = vpa(subs(L_Ib_3, {tau, v}, {tau_val, v_val_vec(index)}));

                    S_vec_1(i) = L_Ib_1_mc;
                    S_vec_2(i) = L_Ib_2_mc;
                    S_vec_3(i) = L_Ib_3_mc;
                    
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
            
            P_b_vec(j) = 3*p1 - 3*p2 + p3;
        end

        %% Save P_b_vec
        csvwrite('D:\Satellites\28GHz\data\coverage\P_bn_'...
            +string(int32(h*1e-3))+'_'+string(N)+'.csv', P_b_vec);
      
%     end
% end


%% Figure
% figure;
% plot(tau_vec_dB, P_b_vec, 's-', 'LineWidth', 1.5)
% xlabel('SIR [dB]')
% ylabel('Coverage Probability')
% title('h = '+ string(int32(h.*1e-3)) +' km, N = '+string(N))
% ylim([0 1])
% % xlim([-20 20])


