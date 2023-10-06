% exclusion zone

N_vec = [100, 1000]; % number of satellites
h_vec = [500, 1000, 1500, 2000]*1e3; % orbit altitude

G_vec_read = csvread("D:\Satellites\28GHz\data\gain.csv"); % dB scale or dBm? [DOUBT]
G_vec = 10.^(G_vec_read./10.0); % absolute scale
vec_alpha = linspace(0, pi);

h = h_vec(2);

R_vec_final = [];

for N= N_vec
    f_alpha_vec = real(csvread('D:\Satellites\28GHz\data\alpha0\f_alpha_'+string(int32(h.*1e-3)) ...
                +'_'+string(N)+'.csv'));

    f_alpha_vec = f_alpha_vec./sum(f_alpha_vec);

    G_mean = dot(G_vec, f_alpha_vec);

    lambda_s = 1e-6;
    zeta = 4;
    p_s = 19.9;    % power satellite [43 dBm]


    %%

    I_s_level_vec = [-80:1:-50];
    R_solved_vec = [];

    for i= 1:length(I_s_level_vec)
        I_level = 10^(I_s_level_vec(i)/10.0);
        R_solved_vec = [R_solved_vec; ((zeta-2)*I_level/(2*pi*lambda_s*p_s*G_mean))^(1/(2-zeta))] ;
    end

    R_vec_final = [R_vec_final; R_solved_vec'];
    
 end
%%

color_vec = ["219ebc", "fb8500"];

figure;
hold on;
for i=1:length(N_vec)
    c_this = [hex2rgb(color_vec(i))];
    plot(I_s_level_vec, R_vec_final(i,:), 'Color' , c_this, 'LineWidth', 1.5);
end
xlabel('Interference I_s(R)  [dB]');
ylabel('Exclusion Zone radius R [m]');
grid on;
box on;
lg = legend(string(N_vec))
title(lg, "N")