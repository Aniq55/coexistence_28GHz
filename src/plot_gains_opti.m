d = csvread("D:/Satellites/28GHz/data/bandwidth_gains.csv");

factor = d(:,1);

gain_Dc_sum = d(:,2);
gain_Db_sum = d(:,3);
gain_Dc_prod = d(:,4);
gain_Db_prod = d(:,5);

%%

color_vec = ["219ebc", "023047", "ffb703", "fb8500" ];
c_1 = [hex2rgb(color_vec(1))];
c_2 = [hex2rgb(color_vec(2))];
c_3 = [hex2rgb(color_vec(3))];
c_4 = [hex2rgb(color_vec(4))];

figure;
hold on;
plot(factor, gain_Dc_sum,  's-', 'Color' , c_4, 'LineWidth', 1);
plot(factor, gain_Dc_prod, '*-', 'Color' , c_1, 'LineWidth', 1);
plot(factor, gain_Db_sum, 'v-', 'Color' , c_4, 'LineWidth', 1);
    plot(factor, gain_Db_prod, 'o-', 'Color' , c_1, 'LineWidth', 1);
box on;
grid on;
legend(["$ D_c, f_\Sigma $", "$ D_c, f_\Pi $", ...
    "$ D_b, f_\Sigma $", "$ D_b, f_\Pi $"], 'Interpreter', 'latex')
xlabel("B/B_{cn}")
ylabel("Data rate improvement [%]")