clc;
clear all;

%%

lambda_b = 1e-7;
lambda_c = 1e-6;
rho = 100; %m

sigma_b = 1e-11;
sigma_c = 1e-11;

p_c = 1; % W
p_b = 10; % W

alpha = pi./3;
G_0 = 0.1;
G_1 = 10;

beta = 4;
half_beta = 2./beta;

tau_vec_dB = [-10:5:50];
tau_vec_abs = 10.^(tau_vec_dB./10.0);
tau_vec_size = length(tau_vec_dB);

P_c_vec = zeros(1, tau_vec_size);
P_b_vec = zeros(1, tau_vec_size);

for j=1: tau_vec_size

    tau = tau_vec_abs(j);

    %% Coverage probability [Cellular]

    I_0_pre_int = @(x) 1./(1 + x.^(1./half_beta));
    I_0 = integral(@(x) I_0_pre_int(x), tau.^(-half_beta), Inf);

    I_1_pre_int = @(r) exp( -(sigma_c.^2).*tau.*(r.^beta)./p_c ...
        - 0.5.*( ((tau./p_c).^half_beta).*lambda_b.*( ...
        alpha.*G_1.^half_beta + (2.*pi - alpha).*G_0.^half_beta )./sinc(half_beta)).*(p_b.^half_beta).*r.^2 ...
        - pi.*lambda_c.*exp(-pi.*lambda_b.*rho.^2).*(1 + 0.5.*(tau.^half_beta).*I_0 ).*r.^2 ).*r;

    I_1 = integral( @(r) I_1_pre_int(r), 0, Inf);

    P_c = 2.*pi.*lambda_c.*exp(-pi.*lambda_b.*rho.^2).*I_1;

    P_c_vec(j) = P_c;

    %% Coverage probability [Backhaul]

    I_3_pre_int = @(x, r) x./(1.0 + ((x./r).^beta)./(tau.*G_0./G_1) );
    I_3 = @(r) integral(@(x) I_3_pre_int(x,r), r, Inf);

    I_4_pre_int = @(x,r) x./(1.0 + (1./tau).*(x./r).^beta );
    I_4 = @(r) integral( @(x) I_4_pre_int(x,r), r, Inf);

    I_5_pre_int = @(x,r) x./( 1.0 + (G_1.*(x./r).^beta)./(tau.*p_c./p_b) );
    I_5 = @(r) integral( @(x) I_5_pre_int(x,r), rho, Inf);

    %% Monte-carlo

    f = @(r) exp( ... 
        - sigma_b.^2.*(tau.*r.^beta)./(p_b.*G_1) ...
        - (2.*pi - alpha).*lambda_b.*I_3(r) - alpha.*lambda_b.*I_4(r) ...
        - 2.*pi.*lambda_c.*exp(- pi.*lambda_b.*rho.^2).*I_5(r) );

    v_val_vec = [1.0: 1.0: 1e6];
    F_v = @(z) exp(-pi.*lambda_c.*z.^2);
    F_v_vec = F_v(v_val_vec);

    iter = 1e3;

    S_vec = zeros(1, iter);

    for i= 1:iter
        u = rand(1);
        [min_val, index] = min(abs(F_v_vec - u));
        try   
            S_vec(i) = f(v_val_vec(index));
        catch
            S_vec(i) = NaN;
        end
    end

    P_b = mean(S_vec(~isnan(S_vec)));
    P_b_vec(j) = P_b;

end

%%
color_vec = ["219ebc", "fb8500" ];
c_1 = [hex2rgb(color_vec(1))];
c_2 = [hex2rgb(color_vec(2))];

figure;
hold on;
plot(tau_vec_dB, P_c_vec, 'Color' , c_1, 'LineWidth', 1.5);
plot(tau_vec_dB, P_b_vec, 'Color' , c_2, 'LineWidth', 1.5);
box on; grid on;
legend(['P_c'; 'P_b']);
ylabel('Coverage Probability');
xlabel('SINR [dB]');






