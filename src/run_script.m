% tic
%     run SIM_P_b.m
%     run SIM_P_bn.m
%     run SIM_P_c.m
%     run SIM_P_cn.m
% toc

tic
    run THE_P_b.m
    run THE_P_bn.m
    run THE_P_c.m
    run THE_P_cn.m
toc

run plot_cov_naka.m
