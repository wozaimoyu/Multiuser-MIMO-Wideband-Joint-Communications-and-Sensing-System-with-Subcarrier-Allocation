clc;
clear all;
close all

%% system parameters
Nt = 32; % number of transmit antennas
M = 4; % number of users
K = 64; % number of subcarriers

%% simulation parameters
n_chans = 100;
J = 32;
Nrf_vec = [1,2:2:32];
n_Nrf = length(Nrf_vec);
HBF = 1;
plot_MSE = 1; plot_rate = 1; plot_tradeoff = 1;
use_ZF = 0;
SNR_dB = 10; Pt = db2pow(SNR_dB);
CSI_error = 0;

T = 181;
markers = {'','','',''};
rho_vec = 0.4;%[0.2, 0.4, 0.6];
for rr = 1:length(rho_vec)
    rho = rho_vec(rr);
    
    % proposed scheme
    rate_dig_comm = zeros(n_chans,n_Nrf);
    rate_dig_propose = zeros(n_chans,n_Nrf);
    rate_dig_overlap = zeros(n_chans,n_Nrf);
    rate_dig_non_overlap = zeros(n_chans,n_Nrf);
    
    rate_hyb_comm = zeros(n_chans,n_Nrf);
    rate_hyb_propose = zeros(n_chans,n_Nrf);
    rate_hyb_overlap = zeros(n_chans,n_Nrf);
    rate_hyb_non_overlap = zeros(n_chans,n_Nrf);
    
    %% start simulations for SNRs
    for ss = 1:n_Nrf
        ss
        Nrf = Nrf_vec(ss);
        F0 = exp(1i*randn(Nt,Nrf));
       % load data
        [Q0, Rate, theta, H, C0, Pd_theta, a] = load_data(Nt,M,K,SNR_dB);
        
        %% start simulations for channels
        parfor nn = 1:n_chans
            nn
            h = H(:,:,:,nn); % current channel
            
            % get Q and rate
            Q0_nn = Q0(:,:,:,nn);
            Rate_nn = Rate(:,nn);
            
            % choose best subcarriers
            [Omg, Omg_rand] = subcarrier_select(K,J,Q0_nn,C0);
            
            % obtain beampatterns, MSE, and rate
            [~,~, rate_dig_overlap(nn,ss),...
                ~, ~, rate_dig_non_overlap(nn,ss),...
                ~, ~, rate_dig_propose(nn,ss), ...
                ~, ~, rate_hyb_overlap(nn,ss),...
                ~, ~, rate_hyb_non_overlap(nn,ss),...
                ~, ~, rate_hyb_propose(nn,ss),...
                ~, ~,...
                rate_dig_comm(nn,ss),rate_hyb_comm(nn,ss), ~, ~]...
                = JCAS_design(Nt,M,K,C0,Q0_nn,Pt,Pd_theta,Omg,Omg_rand,rho,a,h,T,F0,HBF,CSI_error);
        end
    end
    
    % legends
    schemes = {'Communication only', 'Prop. JCAS', 'Conv. JCAS, overlap', 'Conv. JCAS, nonoverlap'};
    
    %% plot rate ============================================================================================
    %% compute everage rate and MSE
    rate_dig_comm_mean = mean(rate_dig_comm,1);
    rate_dig_propose_mean = mean(rate_dig_propose,1);
    rate_dig_overlap_mean = mean(rate_dig_overlap,1);
    rate_dig_non_overlap_mean = mean(rate_dig_non_overlap,1);
    
    rate_hyb_comm_mean = mean(rate_hyb_comm,1);
    rate_hyb_propose_mean = mean(rate_hyb_propose,1);
    rate_hyb_overlap_mean = mean(rate_hyb_overlap,1);
    rate_hyb_non_overlap_mean = mean(rate_hyb_non_overlap,1);
    
    
    figure
    plot(Nrf_vec, rate_dig_comm_mean, '-k+','LineWidth',2,'MarkerSize',8); hold on;
    plot(Nrf_vec, rate_hyb_comm_mean, '-ks','LineWidth',2,'MarkerSize',8); hold on;
    
    plot(Nrf_vec, rate_dig_propose_mean, '-r*','LineWidth',2,'MarkerSize',8); hold on;
    plot(Nrf_vec, rate_hyb_propose_mean, '-rp','LineWidth',2,'MarkerSize',8); hold on;
    
    plot(Nrf_vec, rate_dig_overlap_mean, ':bo','LineWidth',2,'MarkerSize',8); hold on;
    plot(Nrf_vec, rate_hyb_overlap_mean, ':bh','LineWidth',2,'MarkerSize',8); hold on;
    
    plot(Nrf_vec, rate_dig_non_overlap_mean, '-.^','Color',[0.4660 0.6740 0.1880]','LineWidth',2,'MarkerSize',8); hold on;
    plot(Nrf_vec, rate_hyb_non_overlap_mean, '-.>','Color',[0.4660 0.6740 0.1880]','LineWidth',2,'MarkerSize',8); hold on;
    
    
    
    xlabel('Number of RF chains ($N_{\mathrm{RF}}$)','fontsize',12,'interpreter','latex');
    ylabel('Total achievable rate [bits/s/Hz]','fontsize',12,'interpreter','latex');
    xticks(Nrf_vec)
    legend('DBF, Communication only',...
        'HBF, Communication only',...
        'DBF, Prop. JCAS',...
        'HBF, Prop. JCAS',...
        'DBF, Conv. JCAS, overlap',...
        'HBF, Conv. JCAS, overlap',...
        'DBF, Conv. JCAS, nonoverlap',...
        'HBF, Conv. JCAS, nonoverlap',...
        'Location','Best','fontsize',10,'interpreter','latex')
    xlim([1,32])
    xticks(Nrf_vec)
    grid on
end