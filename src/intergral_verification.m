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
        
        lambda_c_tilde = lambda_c*exp( -(pi*lambda_c*(R_s - 0.5/sqrt(lambda_c))^2)*...
            double(R_s > 0.5/sqrt(lambda_c)) ...
            -(pi*lambda_b*(R_b - 0.5/sqrt(lambda_c))^2)*...
            double(R_b > 0.5/sqrt(lambda_c)));

        lambda_b_tilde = lambda_b*exp( -(pi*lambda_b*( R_b - 0.5/sqrt(lambda_b))^2)*...
            double(R_b > 0.5/sqrt(lambda_b)));

        syms tau v x G

        %% Interference from Earth Stations
        % [OK]
        S_1_1 = ( 1 + (1./mu).*G.*beta.*tau.*(p_s./(p_b.*G1)).*(v./x).^zeta )^-m;    
        S_1_2 = ( 1 + (2./mu).*G.*beta.*tau.*(p_s./(p_b.*G1)).*(v./x).^zeta )^-m;        
        S_1_3 = ( 1 + (3./mu).*G.*beta.*tau.*(p_s./(p_b.*G1)).*(v./x).^zeta )^-m;        
        
        
        %% 
        I_2_1 =  vpaintegral( (1 - S_1_1).*x, x, [R_s Inf]); % [OK]
        I_2_2 =  vpaintegral( (1 - S_1_2).*x, x, [R_s Inf]); % [OK]
        I_2_3 =  vpaintegral( (1 - S_1_3).*x, x, [R_s Inf]); % [OK]
        
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
        
        
%%
        
        L_pre_Is_1_ = vpa(subs(I_2_1, {G}, {G_vec(1)})).*f_alpha_vec(1).*d_alpha;
        L_pre_Is_2_ = vpa(subs(I_2_2, {G}, {G_vec(1)})).*f_alpha_vec(1).*d_alpha;
        L_pre_Is_3_ = vpa(subs(I_2_3, {G}, {G_vec(1)})).*f_alpha_vec(1).*d_alpha;
        
        
        for j=2:length(G_vec)
            L_Is_1_ = L_Is_1 + vpa(subs(L_pre_Is_1_, {G}, {G_vec(j)})).*f_alpha_vec(j).*d_alpha;
            L_Is_2_ = L_Is_2 + vpa(subs(L_pre_Is_2_, {G}, {G_vec(j)})).*f_alpha_vec(j).*d_alpha;
            L_Is_3_ = L_Is_3 + vpa(subs(L_pre_Is_3_, {G}, {G_vec(j)})).*f_alpha_vec(j).*d_alpha;
        end
        
        L_Is_1_new = exp(-2.*pi.*lambda_s.*L_Is_1_);
        L_Is_2_new = exp(-2.*pi.*lambda_s.*L_Is_2_);
        L_Is_3_new = exp(-2.*pi.*lambda_s.*L_Is_3_);
        
        
%%


%% Coverage probability:

% Evaluate P_c
% Check: P_c(0) should be 1 [OK]

v_val_vec = [1.0: 1.0: 1e6];
F_v = @(z) exp(-pi.*lambda_c_tilde.*z.^2);
F_v_vec = F_v(v_val_vec);

iter = 1e2;

tau_vec_dB = [-10:30:60];
tau_vec_abs = 10.^(tau_vec_dB./10.0);
P_c_vec = zeros(1, length(tau_vec_abs));
P_c_vec_new = zeros(1, length(tau_vec_abs));

for j = 1: length(tau_vec_abs)
    tau_val = tau_vec_abs(j);

    S_vec_1 = zeros(1, iter);
    S_vec_2 = zeros(1, iter);
    S_vec_3 = zeros(1, iter);
    
    S_vec_1_new = zeros(1, iter);
    S_vec_2_new = zeros(1, iter);
    S_vec_3_new = zeros(1, iter);

    tic
    for i= 1:iter
        u = rand(1);
        [min_val, index] = min(abs(F_v_vec - u));
        try
            L_Is_1_mc = vpa(subs(L_Is_1, {tau, v}, {tau_val, v_val_vec(index)}));
            L_Is_2_mc = vpa(subs(L_Is_2, {tau, v}, {tau_val, v_val_vec(index)}));
            L_Is_3_mc = vpa(subs(L_Is_3, {tau, v}, {tau_val, v_val_vec(index)}));
            
            S_vec_1(i) = L_Is_1_mc;
            S_vec_2(i) = L_Is_2_mc;
            S_vec_3(i) = L_Is_3_mc;
            
            L_Is_1_mc_new = vpa(subs(L_Is_1_new, {tau, v}, {tau_val, v_val_vec(index)}));
            L_Is_2_mc_new = vpa(subs(L_Is_2_new, {tau, v}, {tau_val, v_val_vec(index)}));
            L_Is_3_mc_new = vpa(subs(L_Is_3_new, {tau, v}, {tau_val, v_val_vec(index)}));
            
            S_vec_1_new(i) = L_Is_1_mc_new;
            S_vec_2_new(i) = L_Is_2_mc_new;
            S_vec_3_new(i) = L_Is_3_mc_new;

        catch
            S_vec_1(i) = NaN;
            S_vec_2(i) = NaN;
            S_vec_3(i) = NaN;
            
            S_vec_1_new(i) = NaN;
            S_vec_2_new(i) = NaN;
            S_vec_3_new(i) = NaN;
        end
    end
    toc

    Svec_1 = mean(S_vec_1(~isnan(S_vec_1)));
    Svec_2 = mean(S_vec_2(~isnan(S_vec_2)));
    Svec_3 = mean(S_vec_3(~isnan(S_vec_3)));
    
    Svec_1_new = mean(S_vec_1_new(~isnan(S_vec_1_new)));
    Svec_2_new = mean(S_vec_2_new(~isnan(S_vec_2_new)));
    Svec_3_new = mean(S_vec_3_new(~isnan(S_vec_3_new)));
    
    P_c_vec(j) = Svec_1;
    P_c_vec_new(j) = Svec_1_new;
    
end
        
%%
figure;
hold on;
plot(tau_vec_dB, P_c_vec, 's-');
plot(tau_vec_dB, P_c_vec_new, 'o-');
legend('E_G exp', 'exp E_G')

        
        
        
        
        