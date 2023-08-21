clc;
clear all;
close all

%% system parameters
Nt = 32; % number of transmit antennas
M = 4; % number of users
K = 64; % number of subcarriers

%% simulation parameters
n_chans = 12;
J_vec = [1,8:8:64];
n_J = length(J_vec);
plot_MSE = 1; plot_rate = 1; plot_tradeoff = 1;
CSI_error = 0;
HBF = 1;

SNR_dB = 12; Pt = db2pow(SNR_dB);

T = 181;
Nrf = 6; % number of RF chains
F0 = exp(1i*randn(Nt,Nrf));
rho_vec = [0.2, 0.4, 0.6];
for rr = 1:length(rho_vec)
    rho = rho_vec(rr);
    
    % rate
    rate_hyb_propose = zeros(n_chans,n_J);
    rate_hyb_overlap = zeros(n_chans,n_J);
    rate_hyb_non_overlap = zeros(n_chans,n_J);
    rate_hyb_comm = zeros(n_chans,n_J);
    rate_hyb_propose_rand = zeros(n_chans,n_J);
    % MSE
    MSE_hyb_propose = zeros(n_chans,n_J);
    MSE_hyb_overlap = zeros(n_chans,n_J);
    MSE_hyb_non_overlap = zeros(n_chans,n_J);
    MSE_hyb_comm = zeros(n_chans,n_J);
    
    
    %% start simulations for J
    for ss = 1:n_J
        ss
        J = J_vec(ss);
        % load data
        [Q0, Rate, theta, H, C0, Pd_theta, a] = load_data(Nt,M,K,SNR_dB);
        
        %% start simulations for channels
        parfor nn = 1:n_chans
            h = H(:,:,:,nn); % current channel
            
            % get Q and rate
            Q0_nn = Q0(:,:,:,nn);
            Rate_nn = Rate(:,nn);
            
            % choose best subcarriers
            [Omg, Omg_rand] = subcarrier_select(K,J,Q0_nn,C0);
            
            % obtain beampatterns, MSE, and rate
            [~, ~, ~,...
                ~,~,~,...
                ~,~,~,...
                ~, MSE_hyb_overlap(nn,ss), rate_hyb_overlap(nn,ss),...
                ~, MSE_hyb_non_overlap(nn,ss), rate_hyb_non_overlap(nn,ss),...
                ~, MSE_hyb_propose(nn,ss), rate_hyb_propose(nn,ss),...
                ~, ~,...
                ~,rate_hyb_comm(nn,ss), ~, MSE_hyb_comm(nn,ss), ~]...
                = JCAS_design(Nt,M,K,C0,Q0_nn,Pt,Pd_theta,Omg,Omg_rand,rho,a,h,T,F0,HBF,CSI_error);
        end
    end
    
    % legends
    schemes = {'Communication only', 'Prop. JCAS', 'Conv. JCAS, overlap', 'Conv. JCAS, nonoverlap'};
    
    %% plot rate ============================================================================================
    if plot_rate == 1
        % compute rate
        rate_hyb_comm_mean = mean(rate_hyb_comm,1);% comm only
        rate_hyb_propose_mean = mean(rate_hyb_propose,1);
        rate_hyb_overlap_mean = mean(rate_hyb_overlap,1);
        rate_hyb_non_overlap_mean = mean(rate_hyb_non_overlap,1);
        % plot figure
        figure(1)
        plot(J_vec, rate_hyb_comm_mean, '-k+','LineWidth',2,'MarkerSize',8); hold on;
        plot(J_vec, rate_hyb_propose_mean, '-r*','LineWidth',2,'MarkerSize',8); hold on;
        plot(J_vec, rate_hyb_overlap_mean, ':bo','LineWidth',2,'MarkerSize',8); hold on;
        plot(J_vec, rate_hyb_non_overlap_mean,  '-.^','Color',[0.4660 0.6740 0.1880]','LineWidth',2,'MarkerSize',8); hold on;
        xlabel('Number of JCAS subcarriers $(J)$','fontsize',12,'interpreter','latex');
        ylabel('Achievable rate [bits/s/Hz]','fontsize',12,'interpreter','latex');
        xlim([J_vec(1) J_vec(n_J)])
        xticks(J_vec)
        legend(schemes{1},schemes{2},schemes{3},schemes{4},'Location','Best','fontsize',11,'interpreter','latex')
        grid on
    end
    
    %% plot MSE ============================================================================================
    if plot_MSE == 1
        % compute MSE
        MSE_mean_comm = mean(MSE_hyb_comm,1);
        MSE_mean_hyb_propose = mean(MSE_hyb_propose,1);
        MSE_mean_hyb_overlap = mean(MSE_hyb_overlap,1);
        MSE_mean_hyb_non_overlap = mean(MSE_hyb_non_overlap,1);
        % plot figure
        figure(2)
        plot(J_vec, MSE_mean_comm, '-k+','LineWidth',2,'MarkerSize',8); hold on;
        plot(J_vec, MSE_mean_hyb_propose, '-r*','LineWidth',2,'MarkerSize',8); hold on;
        plot(J_vec, MSE_mean_hyb_overlap, ':bo','LineWidth',2,'MarkerSize',8); hold on;
        plot(J_vec, MSE_mean_hyb_non_overlap,  '-.^','Color',[0.4660 0.6740 0.1880],'LineWidth',2,'MarkerSize',8); hold on;
        xlabel('Number of JCAS subcarriers $(J)$','fontsize',12,'interpreter','latex');
        ylabel('Average MSE','fontsize',12,'interpreter','latex');
        legend(schemes{1},schemes{2},schemes{3},schemes{4},'Location','Best','fontsize',11,'interpreter','latex')
        xlim([J_vec(1) J_vec(n_J)])
        xticks(J_vec)
        grid on
    end
    
    if plot_tradeoff == 1
        figure(3)
        plot(MSE_mean_comm, rate_hyb_comm_mean, '-k+','LineWidth',2,'MarkerSize',8); hold on;
        plot(MSE_mean_hyb_propose, rate_hyb_propose_mean, '-r*','LineWidth',2,'MarkerSize',8); hold on;
        plot(MSE_mean_hyb_overlap, rate_hyb_overlap_mean, ':bo','LineWidth',2,'MarkerSize',8); hold on;
        plot(MSE_mean_hyb_non_overlap, rate_hyb_non_overlap_mean,  '-.^','Color',[0.4660 0.6740 0.1880],'LineWidth',2,'MarkerSize',8); hold on;
        xlabel('Average MSE','fontsize',12,'interpreter','latex');
        ylabel('Achievable rate [bits/s/Hz]','fontsize',12,'interpreter','latex');
        %xticks(J_vec)
        legend(schemes{1},schemes{2},schemes{3},schemes{4},'Location','Best','fontsize',11,'interpreter','latex')
        grid on
    end
end