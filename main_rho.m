clc;
clear all;
close all

%% system parameters
Nt = 32; % number of transmit antennas
M = 4; % number of users
K = 64; % number of subcarriers

%% simulation parameters
n_chans = 12;
J = 32;
rho_vec = 0:0.2:1;
n_rho = length(rho_vec);
HBF = 1;
plot_MSE = 1; plot_rate = 1; plot_tradeoff = 1;
use_ZF = 0;
SNR_dB = 12; Pt = db2pow(SNR_dB);

T = 181;
idx = 1;
markers = {'','','',''};
Nrf = 6; % number of RF chains
F0 = exp(1i*randn(Nt,Nrf));
tau = zeros(n_chans,n_rho);

for ss = 1:n_rho
    ss
    rho = rho_vec(ss);
    % load data
    [Q0, Rate, theta, H, C0, Pd_theta, a] = load_data(Nt,M,K,SNR_dB);
    
    %% start simulations for channels
    parfor nn = 1:n_chans
        %nn
        h = H(:,:,:,nn); % current channel
        
        % get Q and rate
        if use_ZF == 1  % by ZF solution
            use_waterfilling = 1;
            [Q0_nn, Rate_nn] = ZF_beamforming(Nt,M,K,h,Pt,use_waterfilling);
        else % by SCA
            Q0_nn = Q0(:,:,:,nn);
            Rate_nn = Rate(:,nn);
            %rate_comm = mean(Rate_nn)
        end
        
        % choose best subcarriers
        [Omg, Omg_rand] = subcarrier_select(K,J,Q0_nn,C0);
        
        % obtain beampatterns, MSE, and rate
        tau(nn,ss) = JCAS_design_tau(Nt,M,K,C0,Q0_nn,Pt,Omg,rho,F0,HBF);
    end
end
for nn = 1:n_chans
    plot(rho_vec, tau(nn,:), '*','LineWidth',1.5,'MarkerSize',6); hold on;
end
tau_mean = mean(tau,1);
plot(rho_vec, tau_mean, '-r','LineWidth',2,'MarkerSize',8); hold on;
xlabel('$\rho$','fontsize',12,'interpreter','latex');
ylabel('$\tau$','fontsize',12,'interpreter','latex');