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

        %% Interference from Earth Stations
        % [OK]
        S_1_1 = ( 1 + (1./mu).*G.*beta.*tau.*(p_s./p_c).*(v./x).^zeta )^-m;    
        S_1_2 = ( 1 + (2./mu).*G.*beta.*tau.*(p_s./p_c).*(v./x).^zeta )^-m;        
        S_1_3 = ( 1 + (3./mu).*G.*beta.*tau.*(p_s./p_c).*(v./x).^zeta )^-m;        
        
        
        %% 
        I_2_1 =  vpaintegral( (1 - S_1_1).*x, x, [0 Inf]); % [OK]
        I_2_2 =  vpaintegral( (1 - S_1_2).*x, x, [0 Inf]); % [OK]
        I_2_3 =  vpaintegral( (1 - S_1_3).*x, x, [0 Inf]); % [OK]
        
        %%       
        L_pre_Is_1 = exp(-2.*pi.*lambda_s.*I_2_1); % [OK]
        L_pre_Is_2 = exp(-2.*pi.*lambda_s.*I_2_2); % [OK]
        L_pre_Is_3 = exp(-2.*pi.*lambda_s.*I_2_3); % [OK]
        
        %% Take average over G
        
        L_Is_1 = vpa(subs(L_pre_Is_1, {G}, {G_vec(1)})).*f_alpha_vec(1).*d_alpha;
        L_Is_2 = vpa(subs(L_pre_Is_2, {G}, {G_vec(1)})).*f_alpha_vec(1).*d_alpha;
        L_Is_3 = vpa(subs(L_pre_Is_3, {G}, {G_vec(1)})).*f_alpha_vec(1).*d_alpha;
        
        for j=2:length(G_vec)
            L_Is_1 = L_Is_1 + vpa(subs(L_pre_Is_1, {G}, {G_vec(j)})).*f_alpha_vec(j).*d_alpha;
            L_Is_2 = L_Is_2 + vpa(subs(L_pre_Is_2, {G}, {G_vec(j)})).*f_alpha_vec(j).*d_alpha;
            L_Is_3 = L_Is_3 + vpa(subs(L_pre_Is_3, {G}, {G_vec(j)})).*f_alpha_vec(j).*d_alpha;
        end
        
        %% Interference from Backhaul Points
        S_b_1_g1 = (1 + 1.*G1.*(p_b./p_c).*beta.*tau.*(v./x).^zeta ).^-m;
        S_b_1_g0 = (1 + 1.*G0.*(p_b./p_c).*beta.*tau.*(v./x).^zeta ).^-m;
        
        S_b_2_g1 = (1 + 2.*G1.*(p_b./p_c).*beta.*tau.*(v./x).^zeta ).^-m;
        S_b_2_g0 = (1 + 2.*G0.*(p_b./p_c).*beta.*tau.*(v./x).^zeta ).^-m;
        
        S_b_3_g1 = (1 + 3.*G1.*(p_b./p_c).*beta.*tau.*(v./x).^zeta ).^-m;
        S_b_3_g0 = (1 + 3.*G0.*(p_b./p_c).*beta.*tau.*(v./x).^zeta ).^-m;
        
        
        I_b_1_g1 =  vpaintegral( (1 - S_b_1_g1).*x, x, [v Inf]); % [OK]
        I_b_1_g0 =  vpaintegral( (1 - S_b_1_g0).*x, x, [v Inf]); % [OK]
        
        I_b_2_g1 =  vpaintegral( (1 - S_b_2_g1).*x, x, [v Inf]); % [OK]
        I_b_2_g0 =  vpaintegral( (1 - S_b_2_g0).*x, x, [v Inf]); % [OK]
        
        I_b_3_g1 =  vpaintegral( (1 - S_b_3_g1).*x, x, [v Inf]); % [OK]
        I_b_3_g0 =  vpaintegral( (1 - S_b_3_g0).*x, x, [v Inf]); % [OK]
        
        L_Ib_1 = (Theta/pi).*exp(-2.*pi.*lambda_b.*I_b_1_g1) ...
            + (1 - Theta/pi).*exp(-2.*pi.*lambda_b.*I_b_1_g0); % [OK]
        
        L_Ib_2 = (Theta/pi).*exp(-2.*pi.*lambda_b.*I_b_2_g1) ...
            + (1 - Theta/pi).*exp(-2.*pi.*lambda_b.*I_b_2_g0); % [OK]
        
        L_Ib_3 = (Theta/pi).*exp(-2.*pi.*lambda_b.*I_b_3_g1) ...
            + (1 - Theta/pi).*exp(-2.*pi.*lambda_b.*I_b_3_g0); % [OK]

        %% Interference from Cellular Base Stations
        S_2_1 =  (1 - (1 + (1/mu).*(p_c./(p_b.*G1)).*beta.*tau.*(v./x).^zeta ).^-m ).*x; % [OK]
        S_2_2 =  (1 - (1 + (2/mu).*(p_c./(p_b.*G1)).*beta.*tau.*(v./x).^zeta ).^-m ).*x; % [OK]
        S_2_3 =  (1 - (1 + (3/mu).*(p_c./(p_b.*G1)).*beta.*tau.*(v./x).^zeta ).^-m ).*x; % [OK]
        
        I_3_1 = vpaintegral(S_2_1, x, [0 Inf]); % [OK]
        I_3_2 = vpaintegral(S_2_2, x, [0 Inf]); % [OK]
        I_3_3 = vpaintegral(S_2_3, x, [0 Inf]); % [OK]
        
        L_Ic_1 = exp(-2.*pi.*lambda_c.*I_3_1); % [OK]
        L_Ic_2 = exp(-2.*pi.*lambda_c.*I_3_2); % [OK]
        L_Ic_3 = exp(-2.*pi.*lambda_c.*I_3_3); % [OK]

        %% Coverage probability:

        % Evaluate P_b
        % Check: P_b(0) should be 1 [OK]

        v_val_vec = [1.0: 1.0: 1e6];
        F_v = @(z) exp(-pi.*lambda_b.*z.^2);
        F_v_vec = F_v(v_val_vec);

        iter = 1e4;

        tau_vec_dB = [-50:10:20];
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
                    L_Is_1_mc = vpa(subs(L_Is_1, {tau, v}, {tau_val, v_val_vec(index)}));
                    L_Ic_1_mc = vpa(subs(L_Ic_1, {tau, v}, {tau_val, v_val_vec(index)}));
                    L_Ib_1_mc = vpa(subs(L_Ib_1, {tau, v}, {tau_val, v_val_vec(index)}));
                    
                    L_Is_2_mc = vpa(subs(L_Is_2, {tau, v}, {tau_val, v_val_vec(index)}));
                    L_Ic_2_mc = vpa(subs(L_Ic_2, {tau, v}, {tau_val, v_val_vec(index)}));
                    L_Ib_2_mc = vpa(subs(L_Ib_2, {tau, v}, {tau_val, v_val_vec(index)}));
                    
                    L_Is_3_mc = vpa(subs(L_Is_3, {tau, v}, {tau_val, v_val_vec(index)}));
                    L_Ic_3_mc = vpa(subs(L_Ic_3, {tau, v}, {tau_val, v_val_vec(index)}));
                    L_Ib_3_mc = vpa(subs(L_Ib_3, {tau, v}, {tau_val, v_val_vec(index)}));

                    S_vec_1(i) = L_Is_1_mc*L_Ic_1_mc*L_Ib_1_mc;
                    S_vec_2(i) = L_Is_2_mc*L_Ic_2_mc*L_Ib_2_mc;
                    S_vec_3(i) = L_Is_3_mc*L_Ic_3_mc*L_Ib_3_mc;
                    
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
        csvwrite('D:\Satellites\28GHz\data\coverage\P_b_'+string(int32(h*1e-3))+'_'+string(N)+'.csv', ...
                    P_b_vec);

      
%     end
% end


%% Figure
figure;
plot(tau_vec_dB, P_b_vec, 's-', 'LineWidth', 1.5)
xlabel('SIR [dB]')
ylabel('Coverage Probability')
title('h = '+ string(int32(h.*1e-3)) +' km, N = '+string(N))
ylim([0 1])
% xlim([-20 20])

