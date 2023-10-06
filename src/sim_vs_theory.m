N_vec = [100, 500, 1000, 1500, 2000]; % number of satellites
h_vec = [500, 1000, 1500, 2000]*1e3; % orbit altitude
color_vec = ["8ecae6", "219ebc", "023047", "ffb703", "fb8500" ];

%%
N = N_vec(1);
h = h_vec(1);

psi_sim = csvread('D:\Satellites\28GHz\data\psi0\'+string(int32(h*1e-3))+'_'+string(N)+'.csv');
[F_psi_sim, psi_vec_sim] = ecdf(psi_sim);

F_psi = csvread('D:\Satellites\28GHz\data\psi0\F_psi_'+string(int32(h*1e-3))+'_'+string(N)+'.csv');
psi_vec = csvread('D:\Satellites\28GHz\data\psi0\psi_vec_'+string(int32(h*1e-3))+'_'+string(N)+'.csv');


figure;
plot(psi_vec_sim, F_psi_sim, ...
    'rx-', 'LineWidth', 1, 'MarkerIndices', 1:3000:length(psi_vec_sim));
hold on;
plot(rad2deg(psi_vec), F_psi, 'k-', 'LineWidth', 1.5)

legend('Simulation', 'Theoretical', 'Location',"southeast")
xlabel('\psi [degree]')
ylabel('F_{\psi_0}(\psi)')

xlim([0 90])
ylim([0 1])
grid on;
box on;

%% PSI: Fix h, vary N

h = h_vec(3);

figure;
hold on;
for j=1: length(N_vec)
    N= N_vec(j);
    c_this = [hex2rgb(color_vec(j))];
    F_psi = csvread('D:\Satellites\28GHz\data\psi0\F_psi_'+string(int32(h*1e-3))+'_'+string(N)+'.csv');
    psi_vec = csvread('D:\Satellites\28GHz\data\psi0\psi_vec_'+string(int32(h*1e-3))+'_'+string(N)+'.csv');

    plot(rad2deg(psi_vec), F_psi, 'Color', c_this, 'LineWidth', 1.5)

end

plot([180,181], [0,1], 'k-');
plot([180,181], [0,1], 'kx');


for j=1: length(N_vec)
    N= N_vec(j);
    c_this = [hex2rgb(color_vec(j))];
    psi_sim = csvread('D:\Satellites\28GHz\data\psi0\'+string(int32(h*1e-3))+'_'+string(N)+'.csv');
    [F_psi_sim, psi_vec_sim] = ecdf(psi_sim);

    plot(psi_vec_sim, F_psi_sim, ...
        'x',  'Color', c_this, 'LineWidth', 1, 'MarkerIndices', 1:3000:length(psi_vec_sim));
end

xlabel('\psi [degree]')
ylabel('F_{\psi_0}(\psi)')

xlim([0 30])
ylim([0 1])
grid on;
box on;

lg = legend([string(N_vec), "theory", "simulation"], 'Location',"southeast");
title(lg, "N");

%% THETA: Fix h, vary N

h = h_vec(3);

figure;
hold on;
for j=1: length(N_vec)
    N= N_vec(j);
    c_this = [hex2rgb(color_vec(j))];

    F_theta = csvread('D:\Satellites\28GHz\data\theta0\F_theta_'+string(int32(h*1e-3))+'_'+string(N)+'.csv');
    theta_vec = csvread('D:\Satellites\28GHz\data\theta0\x_theta_'+string(int32(h*1e-3))+'_'+string(N)+'.csv');
    
    plot(rad2deg(theta_vec), F_theta, 'Color', c_this, 'LineWidth', 1.5)
end

plot([180,181], [0,1], 'k-');
plot([180,181], [0,1], 'kx');

for j=1: length(N_vec)
    N= N_vec(j);
    c_this = [hex2rgb(color_vec(j))];
    
    theta_sim = csvread('D:\Satellites\28GHz\data\theta0\'+string(int32(h*1e-3))+'_'+string(N)+'.csv');
    [F_theta_sim, theta_vec_sim] = ecdf(theta_sim);

    plot(theta_vec_sim, F_theta_sim, ...
        'x',  'Color', c_this, 'LineWidth', 1, 'MarkerIndices', 1:3000:length(theta_vec_sim));
   
end

xlabel('\theta [degree]')
ylabel('F_{\theta_0}(\theta)')

xlim([0 90])
ylim([0 1])
grid on;
box on;

lg = legend([string(N_vec), "theory", "simulation"], 'Location',"northwest");
title(lg, "N");

%% ALPHA: Fix h, vary N

figure;
hold on;
for j=1: length(N_vec)
    N= N_vec(j);
    c_this = [hex2rgb(color_vec(j))];
    
    f_alpha = csvread('D:\Satellites\28GHz\data\alpha0\f_alpha_'+string(int32(h*1e-3))+'_'+string(N)+'.csv');
    alpha_vec = csvread('D:\Satellites\28GHz\data\alpha0\x_alpha_'+string(int32(h*1e-3))+'_'+string(N)+'.csv');
    
    plot(rad2deg(alpha_vec), f_alpha, 'Color', c_this, 'LineWidth', 1.5); 
end

xlabel('\alpha [degree]');
ylabel('f_{\alpha_0}(\alpha)');

xlim([0 180]);
grid on;
box on;
lg=legend(string(N_vec), 'Location',"northwest");
title(lg, "N");
