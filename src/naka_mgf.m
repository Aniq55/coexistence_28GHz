clear all
clc

m = 3;
beta = (factorial(m))^(-1/m);
v = 1;
zeta = 4;

iters = 1e4;
H = gamrnd(m,1/m,[1,iters]);

syms tau x
MC = @(x, tau) exp(-beta.*1e2.*m.*tau.*(x./v).^(-zeta).*H);

MGF = @(x, tau) (1 + beta.*1e2*tau.*(x/v)^(-zeta))^(-m);

%%

tau_dB_vec = [-5:5:30];
tau_abs_vec = 10.^(tau_dB_vec./10);

T = length(tau_abs_vec);
MC_vec = zeros(1,T);
MGF_vec = zeros(1,T);

for j = 1:T
    tau_val = tau_abs_vec(j);
    
    MC_vec(j) = mean(vpa(subs(MC, {tau, x}, {tau_val, 10})));
    MGF_vec(j) = vpa(subs(MGF, {tau, x}, {tau_val, 10}));   
    
end

%%

figure;
hold on;

plot(tau_dB_vec, MC_vec, 'LineWidth', 1.5);
plot(tau_dB_vec, MGF_vec, 'x--' , 'LineWidth', 1.5);

legend('Monte Carlo', 'MGF')
box on; grid on;
ylabel('E[exp(-sI_s)], s = \beta m \tau v^\zeta / p_c')
xlim([-5 25])
ylim([0 1])
title('Laplace Transform Gap')



