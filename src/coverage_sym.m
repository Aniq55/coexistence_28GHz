clear all
clc

%% Parameters
p_s = 19.9;    % power satellite [43 dBm]
p_c = 0.398;   % power cellular [26 dBm]
zeta = 4;      % pathloss exponent

vec_alpha = linspace(0, pi);
d_alpha = vec_alpha(2) - vec_alpha(1);

% Automation:
N_vec = [100, 500, 1000, 1500, 2000]; % number of satellites
h_vec = [500, 1000, 1500, 2000]*1e3; % orbit altitude

r_e= 6371e3;

lambda_c = 1e-6;
lambda_s = 1e-7;

h = h_vec(1);
% for h = h_vec
    for N  = N_vec
        display(h+"  "+N)
        % PDF [OK]
        f_alpha_vec = real(csvread('D:\Satellites\28GHz\data\alpha0\f_alpha_'+string(int32(h.*1e-3)) ...
            +'_'+string(N)+'.csv'));

        % Check: The PMF must sum to 1
        d_alpha = 1/sum(f_alpha_vec); % 0.0326

        G_vec_read = csvread("D:\Satellites\28GHz\data\gain.csv"); % dB scale or dBm? [DOUBT]
        G_vec = 10.^(G_vec_read./10.0); % absolute scale

        syms tau v x

        %% Interference from Earth Stations
        % [OK]
        S_1 = (f_alpha_vec(1) .* d_alpha )./( 1 + G_vec(1).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(2) .* d_alpha )./( 1 + G_vec(2).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(3) .* d_alpha )./( 1 + G_vec(3).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(4) .* d_alpha )./( 1 + G_vec(4).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(5) .* d_alpha )./( 1 + G_vec(5).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(6) .* d_alpha )./( 1 + G_vec(6).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(7) .* d_alpha )./( 1 + G_vec(7).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(8) .* d_alpha )./( 1 + G_vec(8).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(9) .* d_alpha )./( 1 + G_vec(9).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(10) .* d_alpha )./( 1 + G_vec(10).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(11) .* d_alpha )./( 1 + G_vec(11).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(12) .* d_alpha )./( 1 + G_vec(12).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(13) .* d_alpha )./( 1 + G_vec(13).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(14) .* d_alpha )./( 1 + G_vec(14).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(15) .* d_alpha )./( 1 + G_vec(15).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(16) .* d_alpha )./( 1 + G_vec(16).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(17) .* d_alpha )./( 1 + G_vec(17).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(18) .* d_alpha )./( 1 + G_vec(18).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(19) .* d_alpha )./( 1 + G_vec(19).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(20) .* d_alpha )./( 1 + G_vec(20).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(21) .* d_alpha )./( 1 + G_vec(21).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(22) .* d_alpha )./( 1 + G_vec(22).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(23) .* d_alpha )./( 1 + G_vec(23).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(24) .* d_alpha )./( 1 + G_vec(24).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(25) .* d_alpha )./( 1 + G_vec(25).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(26) .* d_alpha )./( 1 + G_vec(26).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(27) .* d_alpha )./( 1 + G_vec(27).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(28) .* d_alpha )./( 1 + G_vec(28).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(29) .* d_alpha )./( 1 + G_vec(29).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(30) .* d_alpha )./( 1 + G_vec(30).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(31) .* d_alpha )./( 1 + G_vec(31).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(32) .* d_alpha )./( 1 + G_vec(32).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(33) .* d_alpha )./( 1 + G_vec(33).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(34) .* d_alpha )./( 1 + G_vec(34).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(35) .* d_alpha )./( 1 + G_vec(35).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(36) .* d_alpha )./( 1 + G_vec(36).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(37) .* d_alpha )./( 1 + G_vec(37).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(38) .* d_alpha )./( 1 + G_vec(38).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(39) .* d_alpha )./( 1 + G_vec(39).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(40) .* d_alpha )./( 1 + G_vec(40).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(41) .* d_alpha )./( 1 + G_vec(41).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(42) .* d_alpha )./( 1 + G_vec(42).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(43) .* d_alpha )./( 1 + G_vec(43).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(44) .* d_alpha )./( 1 + G_vec(44).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(45) .* d_alpha )./( 1 + G_vec(45).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(46) .* d_alpha )./( 1 + G_vec(46).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(47) .* d_alpha )./( 1 + G_vec(47).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(48) .* d_alpha )./( 1 + G_vec(48).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(49) .* d_alpha )./( 1 + G_vec(49).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(50) .* d_alpha )./( 1 + G_vec(50).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(51) .* d_alpha )./( 1 + G_vec(51).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(52) .* d_alpha )./( 1 + G_vec(52).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(53) .* d_alpha )./( 1 + G_vec(53).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(54) .* d_alpha )./( 1 + G_vec(54).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(55) .* d_alpha )./( 1 + G_vec(55).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(56) .* d_alpha )./( 1 + G_vec(56).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(57) .* d_alpha )./( 1 + G_vec(57).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(58) .* d_alpha )./( 1 + G_vec(58).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(59) .* d_alpha )./( 1 + G_vec(59).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(60) .* d_alpha )./( 1 + G_vec(60).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(61) .* d_alpha )./( 1 + G_vec(61).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(62) .* d_alpha )./( 1 + G_vec(62).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(63) .* d_alpha )./( 1 + G_vec(63).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(64) .* d_alpha )./( 1 + G_vec(64).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(65) .* d_alpha )./( 1 + G_vec(65).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(66) .* d_alpha )./( 1 + G_vec(66).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(67) .* d_alpha )./( 1 + G_vec(67).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(68) .* d_alpha )./( 1 + G_vec(68).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(69) .* d_alpha )./( 1 + G_vec(69).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(70) .* d_alpha )./( 1 + G_vec(70).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(71) .* d_alpha )./( 1 + G_vec(71).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(72) .* d_alpha )./( 1 + G_vec(72).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(73) .* d_alpha )./( 1 + G_vec(73).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(74) .* d_alpha )./( 1 + G_vec(74).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(75) .* d_alpha )./( 1 + G_vec(75).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(76) .* d_alpha )./( 1 + G_vec(76).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(77) .* d_alpha )./( 1 + G_vec(77).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(78) .* d_alpha )./( 1 + G_vec(78).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(79) .* d_alpha )./( 1 + G_vec(79).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(80) .* d_alpha )./( 1 + G_vec(80).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(81) .* d_alpha )./( 1 + G_vec(81).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(82) .* d_alpha )./( 1 + G_vec(82).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(83) .* d_alpha )./( 1 + G_vec(83).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(84) .* d_alpha )./( 1 + G_vec(84).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(85) .* d_alpha )./( 1 + G_vec(85).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(86) .* d_alpha )./( 1 + G_vec(86).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(87) .* d_alpha )./( 1 + G_vec(87).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(88) .* d_alpha )./( 1 + G_vec(88).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(89) .* d_alpha )./( 1 + G_vec(89).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(90) .* d_alpha )./( 1 + G_vec(90).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(91) .* d_alpha )./( 1 + G_vec(91).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(92) .* d_alpha )./( 1 + G_vec(92).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(93) .* d_alpha )./( 1 + G_vec(93).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(94) .* d_alpha )./( 1 + G_vec(94).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(95) .* d_alpha )./( 1 + G_vec(95).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(96) .* d_alpha )./( 1 + G_vec(96).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(97) .* d_alpha )./( 1 + G_vec(97).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(98) .* d_alpha )./( 1 + G_vec(98).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(99) .* d_alpha )./( 1 + G_vec(99).*tau.*(p_s./p_c).*(x./v).^(-zeta) ) ...
            + (f_alpha_vec(100) .* d_alpha )./( 1 + G_vec(100).*tau.*(p_s./p_c).*(x./v).^(-zeta) );


        %%
        I_2 =  vpaintegral( (1 - S_1).*x, x, [0 Inf]); % [OK]

        %%
        L_Is_pre_exp = -2.*pi.*lambda_s.*I_2;
        L_Is = exp(-2.*pi.*lambda_s.*I_2); % [OK]

        %% Interference from Cellular Base Stations
        S_2 =  (1 - 1./(1 + tau.*(x./v).^(-zeta)) ).*x; % [OK]
        I_3 = vpaintegral(S_2, x, [v Inf]); % [OK]
        L_Ic = exp(-2.*pi.*lambda_c.*I_3); % [OK]

        %% Coverage probability:
        % P_c_pre_int = 2.*pi.*lambda_c.*v.*exp(-pi.*lambda_c.*v.^2).*L_Is.*L_Ic; % [OK]
        % P_c = vpaintegral(P_c_pre_int, 0, Inf) % [OK]

        prob_v = 2.*pi.*lambda_c.*v.*exp(-pi.*lambda_c.*v.^2);

        %% Evaluate P_c
        % Check: P_c(0) should be 1 [OK]

        v_val_vec = [1.0: 1.0: 1e6];
        F_v = @(z) exp(-pi.*lambda_c.*z.^2);
        F_v_vec = F_v(v_val_vec);

        iter = 1e3;

        tau_vec_dB = [-10:5:50];
        tau_vec_abs = 10.^(tau_vec_dB./10.0);
        P_c_vec = zeros(1, length(tau_vec_abs));

        tic

        for j = 1: length(tau_vec_abs)
            tau_val = tau_vec_abs(j);
            S_vec = zeros(1, iter);

            for i= 1:iter
                u = rand(1);
                [min_val, index] = min(abs(F_v_vec - u));
                try
                    L_Is_mc = vpa(subs(L_Is, {tau, v}, {tau_val, v_val_vec(index)}));
                    L_Ic_mc = vpa(subs(L_Ic, {tau, v}, {tau_val, v_val_vec(index)}));   
                    S_vec(i) = L_Is_mc*L_Ic_mc;
                catch
                    S_vec(i) = NaN;
                end
            end

            P_c_vec(j) = mean(S_vec(~isnan(S_vec)));
        end

        toc

        %% Coverage
        % tau_vec_dB = [-10:5:50];
        % tau_vec_abs = 10.^(tau_vec_dB./10.0);
        % P_c_vec = P_c(tau_vec_abs);

        %% Save P_c_vec
        csvwrite('D:\Satellites\28GHz\data\coverage\P_c_'+string(int32(h*1e-3))+'_'+string(N)+'.csv', ...
                    P_c_vec);

      
    end
% end


%% Figure
% figure;
% plot(tau_vec_dB, P_c_vec, 'LineWidth', 2)
% xlabel('SIR [dB]')
% ylabel('Coverage Probability')
% title('h = '+ string(int32(h.*1e-3)) +' km, N = '+string(N))

