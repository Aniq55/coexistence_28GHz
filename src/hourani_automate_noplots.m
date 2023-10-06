clc
clear all
sympref('FloatingPointOutput', true);

N_vec = [100, 500, 1000, 1500, 2000]; % number of satellites
h_vec = [500, 1000, 1500, 2000]*1e3; % orbit altitude

syms theta phi1 pi2 
r_e= 6371e3;

for h = h_vec
    for N = N_vec
        gamma = 1 + h/r_e;
        phi1 = acot( (gamma^2*tan(theta) + sqrt(gamma^2*(sec(theta))^2 -1) )/(gamma^2-1)  );
        phi2 = acot( (gamma^2*tan(theta) - sqrt(gamma^2*(sec(theta))^2 -1) )/(gamma^2-1)  );

        d1 = abs(diff(phi1, theta));
        d2 = abs(diff(phi2, theta));

        f_theta(theta) =  (N/2)*sin(phi1)*exp(-(N/2)*(1- cos(phi1)))*d1;

        f_1_theta = matlabFunction(f_theta);
        x_theta = linspace(0, pi/2);
        m = length(x_theta);

        f_theta_val = zeros(1,m);
        for j = [1:m]
            f_theta_val(j) = f_1_theta(x_theta(j));
        end

        G_hourani = zeros(1, length(x_theta));

        for j=[1:m]
            G_hourani(j) = integral(f_1_theta, 0, x_theta(j));
        end

        epsilon = 1e-3;
        f_int_theta_full = integral(@(theta)f_1_theta(theta), 0,  pi/2, 'ArrayValued',true);
        
        psi_vec = [0:0.01:pi/2];
        F_psi = 1 - exp(-(N/2.0)*(1.0-cos(psi_vec)));

        % psi_vec, x_theta, F_psi, G_hourani
        csvwrite('D:\Satellites\28GHz\data\psi0\psi_vec_'+string(int32(h*1e-3))+'_'+string(N)+'.csv', ...
            psi_vec);
        csvwrite('D:\Satellites\28GHz\data\psi0\F_psi_'+string(int32(h*1e-3))+'_'+string(N)+'.csv', ...
            F_psi);
        csvwrite('D:\Satellites\28GHz\data\psi0\x_theta_'+string(int32(h*1e-3))+'_'+string(N)+'.csv', ...
            x_theta);
        csvwrite('D:\Satellites\28GHz\data\psi0\F_theta_'+string(int32(h*1e-3))+'_'+string(N)+'.csv', ...
            G_hourani);

        syms v alpha

        f_alpha_pre_int(v, alpha) = sin(alpha)*f_theta(acos(v))/(pi*sqrt( (1-v^2)*(v^2 - (cos(alpha))^2 ) ));
        f_alpha_pre_int(cos(0.1), 0.1);
        f_alpha_pre_int_mat = @(v1, a1) matlabFunction(f_alpha_pre_int(v1, a1));

        x_alpha = linspace(0, pi);
        m = length(x_alpha);
        f_alpha_val = zeros(1,m);

        for j = [1:m]
            f_alpha_val(j) = double(vpaintegral(f_alpha_pre_int(v, x_alpha(j)),...
                v, abs(cos(x_alpha(j))), 1));
        end

        % x_alpha, f_alpha_val
        csvwrite('D:\Satellites\28GHz\data\alpha0\f_alpha_'+string(int32(h*1e-3))+'_'+string(N)+'.csv', ...
            f_alpha_val);
        csvwrite('D:\Satellites\28GHz\data\alpha0\x_alpha_'+string(int32(h*1e-3))+'_'+string(N)+'.csv', ...
            x_alpha);
    end
end

        