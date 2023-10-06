N_vec = [100, 500, 1000, 1500, 2000]; % number of satellites
h_vec = [500, 1000, 1500, 2000]*1e3; % orbit altitude

r_e = 6371e3;   % radius of the earth
iter = 1e6;

for h = h_vec
    for N_fix = N_vec
        lambda = N_fix/(4*pi*(r_e + h)^2);
        theta_0_vec = zeros(1, iter);
        psi_0_vec = zeros(1, iter);

        j= 1;
        for i =[1:iter]
            N = poissrnd(lambda*4*pi*(r_e + h)^2);
            
            V = unifrnd(-1, 1, 1, N);
            THETA = unifrnd(0, 2*pi, 1, N);

            V_tx = unifrnd(-1,1);
            THETA_tx = unifrnd(0, 2*pi);

            x_sat = (r_e + h)*sqrt(1-V.^2).*cos(THETA);
            y_sat = (r_e + h)*sqrt(1-V.^2).*sin(THETA);
            z_sat = (r_e + h)*V;

            tx_x = r_e*sqrt(1-V_tx.^2).*cos(THETA_tx);
            tx_y = r_e*sqrt(1-V_tx.^2).*sin(THETA_tx);
            tx_z = r_e*V_tx;


            % user vector
            u = [tx_x; tx_y; tx_z];

            % satellite vector
            s = [x_sat; y_sat; z_sat];

            % find the nearest satellite:
            s_minus_u = sum((s - u).^2);
            [minval, idx] = min(s_minus_u);

            s_0 = s(:, idx);
            r_0 = sqrt(s_minus_u(idx));
            psi_0 = acos( dot(u,s_0)/(norm(u)*norm(s_0)) );

            theta_0 = atan( cot(psi_0) - (r_e/(r_e + h))*csc(psi_0) );

            if r_0 <= sqrt( (r_e + h)^2 - r_e^2 )
                psi_0_vec(j) = psi_0;
                theta_0_vec(j) = theta_0;
                j= j+1;
            end
        end

        psi_0_vec_final = psi_0_vec(1:j-1);
        theta_0_vec_final = theta_0_vec(1:j-1);



        % histogram(rad2deg(theta_0_vec), 100)

        csvwrite('D:\Satellites\28GHz\data\theta0\'+string(int32(h*1e-3))+'_'+string(N_fix)+'.csv', ...
            rad2deg(theta_0_vec_final))
        csvwrite('D:\Satellites\28GHz\data\psi0\'+string(int32(h*1e-3))+'_'+string(N_fix)+'.csv', ...
            rad2deg(psi_0_vec_final))

        
    end
end


%% PLOT GLOBE AND SATELLITES
% figure('Position', [10 10 1000 1000]);
% [X,Y,Z] = sphere(25);
% s_plot = surf(X*r_e,Y*r_e,Z*r_e, 'FaceAlpha',0.5);
% s_plot.EdgeColor = '#0f7c91';
% s_plot.FaceColor = '#4fe8f0';
% hold on
% 
% scatter3(x_sat, y_sat, z_sat, 'r.')
% daspect([1 1 1])
% set(gca,'BoxStyle','full','Box','off')
% set(gcf,'Color',[0.90 0.90 0.90])

%%
% 
% plotearth('NEOMap', 'D:\Satellites\28GHz\src\PlotEarthImages\Population')