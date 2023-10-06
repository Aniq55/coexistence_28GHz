% Plot average Interference vs alpha_0 for different R =||x_s - u|| values

lambda_s = 1e-7;
zeta = 4;
p_s = 19.9;    % power satellite [43 dBm]

G_vec_read = csvread("D:\Satellites\28GHz\data\gain.csv"); % dB scale or dBm? [DOUBT]
G_vec = 10.^(G_vec_read./10.0); % absolute scale
vec_alpha = linspace(0, pi);

h = 1000*1e3; N = 2000;
f_alpha_vec = real(csvread('D:\Satellites\28GHz\data\alpha0\f_alpha_'+string(int32(h.*1e-3)) ...
                +'_'+string(N)+'.csv'));
f_alpha_vec = f_alpha_vec./sum(f_alpha_vec);
            
G_mean = dot(G_vec, f_alpha_vec);

R_vec = [50, 100, 250, 500, 1000];
I_s_vec = [];

for i = 1:length(R_vec)
    R= R_vec(i);
    
    fn_int = @(x) exp(-pi.*lambda_s.*x.^2).*x.^(1-zeta);
    int_val = integral(fn_int, R, Inf);

%     I_s = 2*pi*lambda_s*p_s*G_mean*int_val ...
%         + p_s.*G_vec.*R^(-zeta);
    
    I_s = p_s.*G_vec.*R^(-zeta);
    
    I_s_vec = [I_s_vec; I_s];
    
end

%%
color_vec = ["8ecae6", "219ebc", "023047", "ffb703", "fb8500" ];

figure;
hold on;
for i=1:length(R_vec)
    c_this = [hex2rgb(color_vec(i))];
    plot(rad2deg(vec_alpha), 10*log10(I_s_vec(i,:)), 'Color' , c_this, 'LineWidth', 1.5);
end

grid on;
box on;
xlim([40, 90])
lg = legend(string(R_vec))
title(lg, "R [m]")
xlabel("\alpha [degree]")
ylabel("I_s(\alpha_0, R)   [dB]")