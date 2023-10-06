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
Radius_list = [10: 10: 1000];

R_s = 100;
R_b = 100;

% for R_s = Radius_list
%     for R_b = Radius_list
        
        tic

        lambda_c_tilde = lambda_c*exp( -(pi*lambda_c*(R_s - 0.5/sqrt(lambda_c))^2)*...
            double(R_s > 0.5/sqrt(lambda_c)) ...
            -(pi*lambda_b*(R_b - 0.5/sqrt(lambda_c))^2)*...
            double(R_b > 0.5/sqrt(lambda_c)));

        lambda_b_tilde = lambda_b*exp( -(pi*lambda_b*( R_b - 0.5/sqrt(lambda_b))^2)*...
            double(R_b > 0.5/sqrt(lambda_b)));

        %% Simulation

        iter_count= 1e4;
        SIR_vec_c = zeros(1, iter_count);
        SIR_vec_b = zeros(1, iter_count);
        SIR_vec_c_new = zeros(1, iter_count);
        SIR_vec_b_new = zeros(1, iter_count);

        for iter = 1:iter_count

            %ES
            N_ES = poissrnd(lambda_s*Area);
            ES_x = -L + 2*L*rand(1, N_ES);
            ES_y = -L + 2*L*rand(1, N_ES);
            ES_r = abs(ES_x + 1i*ES_y);
            ES_gain = zeros(1, N_ES);
            for i=1:N_ES
                [min_val, index] = min(abs(F_alpha - rand));
                ES_gain(i) = G_vec(index);
            end
            ES_fading = gamrnd(m, 1/m, 1, N_ES);

            I_s_vec = p_s.*ES_gain.*ES_fading.*max(1,ES_r).^(-zeta);
            I_s = sum(I_s_vec);
            ES_R_s_vec = ES_r > R_s;
            I_s_R_s = sum(I_s_vec.*ES_R_s_vec);

            % BP
            N_BP = poissrnd(lambda_b*Area);
            BP_x = -L + 2*L*rand(1, N_BP);
            BP_y = -L + 2*L*rand(1, N_BP);
            BP_r = abs(BP_x + 1i*BP_y);
            BP_gain_vec = [G0, G1];
            BP_gain = BP_gain_vec(1 + int32( 2*pi*rand(1, N_BP) < 2*Theta ));
            BP_fading = gamrnd(m, 1/m, 1, N_BP);
            [r_b, index] = min(BP_r);
            I_b_vec = p_b.*BP_gain.*BP_fading.*max(1,BP_r).^(-zeta);
            S_b = I_b_vec(index);
            I_b = sum(I_b_vec);

            % BP_tilde
            N_BP_tilde = poissrnd(lambda_b_tilde*Area);
            BP_t_x = -L + 2*L*rand(1, N_BP_tilde);
            BP_t_y = -L + 2*L*rand(1, N_BP_tilde);
            BP_t_r = abs(BP_t_x + 1i*BP_t_y);
            BP_t_gain = BP_gain_vec(1 + int32( 2*pi*rand(1, N_BP_tilde) < 2*Theta ));
            BP_t_fading = gamrnd(m, 1/m, 1, N_BP_tilde);
            [r_b_t, index] = min(BP_t_r);
            I_b_t_vec = p_b.*BP_t_gain.*BP_t_fading.*max(1,BP_t_r).^(-zeta);
            S_b_t = (I_b_t_vec(index)/BP_t_gain(index))*G1;
            I_b_t = sum(I_b_t_vec) - I_b_t_vec(index) + S_b_t;
            BP_R_b_vec = BP_t_r > R_b;
            I_b_R_b = sum(I_b_t_vec.*BP_R_b_vec);

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

            % BS_tilde
            N_BS_tilde = poissrnd(lambda_c_tilde*Area);
            BS_t_x = -L + 2*L*rand(1, N_BS_tilde);
            BS_t_y = -L + 2*L*rand(1, N_BS_tilde);
            BS_t_r = abs(BS_t_x + 1i*BS_t_y);
            BS_t_gain = BP_gain_vec(1 + int32( 2*pi*rand(1, N_BS_tilde) < 2*Theta ));
            BS_t_fading = gamrnd(m, 1/m, 1, N_BS_tilde);
            [r_c_t, index] = min(BS_t_r);
            I_c_t_vec = p_c.*BS_t_gain.*BS_t_fading.*max(1,BS_t_r).^(-zeta);
            S_c_t = I_c_t_vec(index);
            I_c_t = sum(I_c_t_vec);


            SIR_vec_b(iter) = S_b_t./(I_s_R_s + I_b_t + I_c_t - S_b_t); % b: I U III
            SIR_vec_c(iter) = S_c_t./(I_s_R_s + I_b_R_b + I_c_t - S_c_t); % c: I
            SIR_vec_b_new(iter) = S_b./(I_b - S_b); % b: II U IV
            SIR_vec_c_new(iter) = S_c./(I_c - S_c); % c: I'

        end

        %% Coverage
        tau_vec_dB = [-80:1:40];
        N_records = length(tau_vec_dB);
        tau_vec_abs = 10.^(tau_vec_dB./10.0);
        P_c_vec = zeros(1, length(tau_vec_abs));
        P_b_vec = zeros(1, length(tau_vec_abs));
        P_c_new_vec = zeros(1, length(tau_vec_abs));
        P_b_new_vec = zeros(1, length(tau_vec_abs));

        for i=1:length(tau_vec_dB)
            tau = tau_vec_abs(i);
            P_c_vec(i) = sum(SIR_vec_c > tau)./length(SIR_vec_c);
            P_b_vec(i) = sum(SIR_vec_b > tau)./length(SIR_vec_b);
            P_c_new_vec(i) = sum(SIR_vec_c_new > tau)./length(SIR_vec_c_new);
            P_b_new_vec(i) = sum(SIR_vec_b_new > tau)./length(SIR_vec_b_new);
        end

        R_s_vec = R_s*ones(N_records,1); 
        R_b_vec = R_b*ones(N_records,1); 

        DATA = [R_s_vec, R_b_vec, P_c_vec', P_b_vec', P_c_new_vec', P_b_new_vec'];

        csvwrite( 'D:\Satellites\28GHz\data\regression\Rs_'...
            +string(int32(R_s))+'_Rb_'+string(int32(R_b))+'.csv', DATA);

        toc
%     end
%     
% end

%% Plotting
% figure;
% hold on;
% box on;
% grid on;
% 
% plot(tau_vec_dB, P_c_new_vec);
% plot(tau_vec_dB, P_b_new_vec);
% plot(tau_vec_dB, P_c_vec);
% plot(tau_vec_dB, P_b_vec);
% 
% ylim([0,1]);
% xlabel('SIR [dB]')
% ylabel('Coverage Probability')
% title('h = '+ string(int32(h.*1e-3)) +' km, N = '+string(N));
% legend('$P_c^{\bar{I}}$', '$P_b^{II \cup IV}$' , '$P_c^I$ (28-GHz)', '$P_b^{I \cup III}$ (28-GHz)', ...
%     'interpreter', 'latex')

%% Poisson Process Plot
% figure('Position', [10 10 500 500]);
% hold on;
% scatter(ES_x, ES_y, 'x');
% scatter(BP_x, BP_y, 'o');
% box on;






        
        