%% USER GROUPING FOR MULTI-USER MODULAR XL-MIMO COMMUNICATIONS
% 
clc; clear; close all
%% A. Near-Field Channel Model for Multi-User Modular XLMIMO Communications
%
% Section II presents a generic near-field array response of modular XL-ULA for a user/scatterer located at (r, θ).
% Based on such array responses, in this section, we consider an uplink multi-user modular XL-MIMO communication
% system in multi-path setting.

%% General Parameters
c = physconst('LightSpeed');        % Speed of light
fc = 2.4e9;                         % Operating frequency
lambda = c/fc;                      % Signal wavelength
d = lambda/2;                       % Inter-element spacing for antennas within each module
d_bar = d/lambda;

beta_0_dB = 30;
beta_0 = 10.^(beta_0_dB/10);        % beta_0 denoting the channel power at the reference distance d0 = 1 m.

N = 32;                             % Number of modules
M = 4;                              % Number of antenna elements within each module
NM = N*M;                           % Total number of array elements
S = (M - 1)*d;                      % Physical size of each module

Gamma = 13;                         % Modular separation parameter (dependent on the discontinuous surface of practical installation structure) Gamma >= M
Gamma_d = Gamma*d;                  % Inter-module distance between the reference elements (say the center elements of different modules)

D = ((N - 1)*Gamma + (M - 1))*d;    % Total physical size of the modular XL-ULA

NN = -(N-1)/2:(N-1)/2;
MM = -(M-1)/2:(M-1)/2;

% L_k + 1 denote the number of channel paths for user k,
% including one possible line-of-sight (LoS) component and Lk non-line-of-sight (NLoS) components.
L_k = 5;                            % Number of scatterers per user

%% Position of the m-th element within module n
y = zeros(1,NM);
indice = 0;
for n = 1:N
    for m = 1:M
        indice = indice + 1;
        position = (NN(n)*Gamma + MM(m))*d;
        y(indice) = position;
    end
end
w = [zeros(1,length(y)); y];

% Position of the m-th element within module n (sub-array based USW model for distinct AoAs/AoDs)
y_n = zeros(1,length(NN));
for n = 1:length(NN)
    y_n(n) = NN(n)*Gamma*d;
end

%% Users and scatterers locations
% (rk,0, θk,0) denote the location of user k, with k belonging to the user set K = {1, ..,K},
% and (rk,l, θk,l), 1 ≤ l ≤ Lk, denote the location of its l-th scatterer.

% Users location
r_k0(1) = 200;                         % 1st User's location
r_k0(2) = 300;                         % 2nd User's location
r_k0(3) = 250;                         % 3rd User's location
theta_k0(1) = 0;                       % 1st User's location
theta_k0(2) = 0;                       % 2nd User's location
theta_k0(3) = deg2rad(45);             % 3rd User's location

K = length(r_k0);                      % K denote the number of single-antenna users

%% Distance between q and the m-th element in module n
% r_n,m = |q - w|

q_k0 = zeros(K,2);
r_nm = zeros(K,length(y));
for k = 1:K
    q_k0(k,:) = [r_k0(k)*cos(theta_k0(k)), r_k0(k)*sin(theta_k0(k))].';
    for i = 1:length(y)
        r_nm(k,i) = norm(q_k0(k,:) - w(:,i));
    end
end

% Scatterers location
% For multipath scenarios, the location (rk,l, θk,l) of the l-th scatterer of user k follows the distributions
% rk,l ∼ U(0, 200 m) and θk,l ∼ U(−π/2 , π/2)
r_kl = zeros(K,L_k);
theta_kl = zeros(K,L_k);
q_kl = zeros(K,L_k,2);
t_kl = zeros(K,L_k);                                    % tk,l is the distance from scatterer l to the user k
for k = 1:K
    for l = 1:L_k
        r_kl(k,l) = rand*200;
        theta_kl(k,l) = -pi/2 + (pi/2 - (-pi/2)).* rand;

        q_kl(k,l,:) = [r_kl(k,l)*cos(theta_kl(k,l)), r_kl(k,l)*sin(theta_kl(k,l))].';

        t_kl(k,l) = norm(q_k0(k,:) - squeeze(q_kl(k,l,:)).');
    end
end

%% USW-based near-field array response vector
% Observation of the user based in the Sub-array based USW model for distinct AoAs/AoDs (Eqs. 5, 6, 7)
a_ARV_USW = zeros(NM,K);
for k = 1:K
    a_ARV_USW(:,k) =  USW_ArrayResponseVector(r_k0(k), theta_k0(k), y_n, N, M, fc);
end

% observation of the scatterer based in the Sub-array based USW model
a_ARV_USW_scatt = zeros(NM,K,L_k);
for k = 1:K
    for l = 1:L_k
        a_ARV_USW_scatt(:,k,l) = USW_ArrayResponseVector(r_kl(k,l), theta_kl(k,l), y_n, N, M, fc);
    end
end

% USW-based channel vector
alpha_k_0 = zeros(1,K);
for k = 1:K
    % alpha_k_0(k) = lambda/(4*pi*r(k));
    alpha_k_0(k) = 1/(r_k0(k));
end

% σk,l represents the radar cross section (RCS) of scatterer l belonging to uniform distribution over [1, 40] m2.
% wk,l is the phase shift following from uniform distribution over [−π, π)
alpha_k_l = zeros(K,L_k);
for k = 1:K
    for l = 1:L_k
        w_kl = -pi + (pi - (-pi)).* rand;
        sigma_kl = (rand*40);
        % alpha_k_l(k,l) = (rand*40)/((4*pi)^(3/2)*t_kl(k,l)*r_kl(k,l))*exp(-1i*2*pi/lambda*t_kl(k,l) + 1i*w_kl);
        alpha_k_l(k,l) = sigma_kl/(t_kl(k,l)*r_kl(k,l))*exp(-1i*2*pi/lambda*t_kl(k,l) + 1i*w_kl);
    end
end

% Based on the array response vector the near-field channel vector for user k can be modelled as
h_NF = zeros(NM,K);
for k = 1:K
    % h_NF(:,k) = sqrt(beta_0)/r(k)*a_ARV_USW(:,k);
    h_NF_temp = sqrt(beta_0)*alpha_k_0(k)*a_ARV_USW(:,k);
    for l = 1:L_k
        h_NF_temp = h_NF_temp + sqrt(beta_0)*alpha_k_l(k,l)*a_ARV_USW_scatt(:,k,l);
    end
    h_NF(:,k) = h_NF_temp;
end

%% MRC Beamforming based in the near-field CSI
h_BF_NF = h_NF;



%% UPW-based far-field array response vector
% Observation of the user based in the in the UPW model
a_ARV_UPW = zeros(NM, K);
for k = 1:K
    a_ARV_UPW(:,k) = UPW_ArrayResponseVector(r_k0(k), theta_k0(k), N, M, Gamma, fc);
end

% Observation of the scatterer based in the in the UPW model
a_ARV_UPW_scatt = zeros(NM, L_k, K);
for k = 1:K
    for l = 1:L_k
        a_ARV_UPW_scatt(:,k,l) = UPW_ArrayResponseVector(r_kl(k,l), theta_kl(k,l), N, M, Gamma, fc);
    end
end

% UPW-based channel vector
% Based on the array response vector the near-field channel vector for user k can be modelled as
h_FF = zeros(NM,K);
for k = 1:K
    % h_FF(:,k) = sqrt(beta_0)/r(k)*a_ARV_UPW(:,k);
    h_FF_temp = sqrt(beta_0)*alpha_k_0(k)*a_ARV_UPW(:,k);
    for l = 1:L_k
        h_FF_temp = h_FF_temp + sqrt(beta_0)*alpha_k_l(k,l)*a_ARV_UPW_scatt(:,k,l);
    end
    h_FF(:,k) = h_FF_temp;
end

%% MRC Beamforming based in the far-field CSI
h_BF_FF = h_FF;


%% SINR and SumRate using MRC, ZF and MMSE beamformings
P_dB = -20:0.1:30;
P = 10.^(P_dB/10);


%% Zero forzing calculation
A_K_NF = Ak_ZeroForcing_Calculation(h_BF_NF);
A_K_FF = Ak_ZeroForcing_Calculation(h_BF_FF);

%% MMSE calculation
C_USW = C_MMSE_Calculation(h_BF_NF, P);
C_UPW = C_MMSE_Calculation(h_BF_FF, P);


SINR_USW_MRC = zeros(length(P),K);
SINR_UPW_MRC = zeros(length(P),K);
SINR_USW_ZF = zeros(length(P),K);
SINR_UPW_ZF = zeros(length(P),K);
SINR_USW_MMSE = zeros(length(P),K);
SINR_UPW_MMSE = zeros(length(P),K);

SumRate_USW_MRC = zeros(1,length(P));
SumRate_UPW_MRC = zeros(1,length(P));
SumRate_USW_ZF = zeros(1,length(P));
SumRate_UPW_ZF = zeros(1,length(P));
SumRate_USW_MMSE = zeros(1,length(P));
SumRate_UPW_MMSE = zeros(1,length(P));
for pp = 1:length(P)
    SINR_USW_MRC(pp,:) = SINR_MRC(P(pp), h_NF, h_BF_NF);
    SINR_UPW_MRC(pp,:) = SINR_MRC(P(pp), h_NF, h_BF_FF);

    SINR_USW_ZF(pp,:) = SINR_ZF(P(pp), A_K_NF, h_NF, h_BF_NF);
    SINR_UPW_ZF(pp,:) = SINR_ZF(P(pp), A_K_FF, h_NF, h_BF_FF);

    SINR_USW_MMSE(pp,:) = SINR_MMSE(P(pp), C_USW(pp,:,:,:), h_NF, h_BF_NF);
    SINR_UPW_MMSE(pp,:) = SINR_MMSE(P(pp), C_UPW(pp,:,:,:), h_NF, h_BF_FF);
    
    SumRate_temp_USW_MRC = 0;
    SumRate_temp_UPW_MRC = 0;
    SumRate_temp_USW_ZF = 0;
    SumRate_temp_UPW_ZF = 0;
    SumRate_temp_USW_MMSE = 0;
    SumRate_temp_UPW_MMSE = 0;
    for k = 1:K
        SumRate_temp_USW_MRC = SumRate_temp_USW_MRC + log2(1 + SINR_USW_MRC(pp,k));
        SumRate_temp_UPW_MRC = SumRate_temp_UPW_MRC + log2(1 + SINR_UPW_MRC(pp,k));

        SumRate_temp_USW_ZF = SumRate_temp_USW_ZF + log2(1 + SINR_USW_ZF(pp,k));
        SumRate_temp_UPW_ZF = SumRate_temp_UPW_ZF + log2(1 + SINR_UPW_ZF(pp,k));

        SumRate_temp_USW_MMSE = SumRate_temp_USW_MMSE + log2(1 + SINR_USW_MMSE(pp,k));
        SumRate_temp_UPW_MMSE = SumRate_temp_UPW_MMSE + log2(1 + SINR_UPW_MMSE(pp,k));
    end
    SumRate_USW_MRC(pp) = SumRate_temp_USW_MRC;
    SumRate_UPW_MRC(pp) = SumRate_temp_UPW_MRC;

    SumRate_USW_ZF(pp) = SumRate_temp_USW_ZF;
    SumRate_UPW_ZF(pp) = SumRate_temp_UPW_ZF;

    SumRate_USW_MMSE(pp) = SumRate_temp_USW_MMSE;
    SumRate_UPW_MMSE(pp) = SumRate_temp_UPW_MMSE;     
end










%% Figuras

F_UsersLocations = figure;
set(F_UsersLocations, 'Position',  [30 580, 560, 420])
set(F_UsersLocations, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(F_UsersLocations, 'defaultLegendInterpreter','latex');
set(F_UsersLocations, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(F_UsersLocations, 'defaultLineLineWidth',1.5);
set(F_UsersLocations, 'color','w');

F_SINR = figure;
set(F_SINR, 'Position',  [580 480, 660, 520])
set(F_SINR, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(F_SINR, 'defaultLegendInterpreter','latex');
set(F_SINR, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(F_SINR, 'defaultLineLineWidth',1.5);
set(F_SINR, 'color','w');

figure(F_UsersLocations);
xline(0,'k')
hold on
yline(0,'k')


for k = 1:K
    figure(F_UsersLocations);
    h1 = plot(q_k0(k,1),q_k0(k,2),'s','MarkerSize',10);
    set(h1, 'markerfacecolor', get(h1, 'color'));
    plot([0, q_k0(k,1)],[ 0,q_k0(k,2)],'Color',get(h1, 'color'))
    nombre = ['U',num2str(k)];
    text(q_k0(k,1),q_k0(k,2),nombre,'VerticalAlignment','bottom','HorizontalAlignment','right')

    plot(q_kl(k,:,1),q_kl(k,:,2),'*','Color',get(h1, 'color'),'MarkerSize',6);

    figure(F_SINR);
    hold on
    plot(P_dB,10*log10(SINR_USW_MRC(:,k)),'-^','Color',get(h1, 'color'),'MarkerIndices',1:40:length(P_dB))
    plot(P_dB,10*log10(SINR_USW_ZF(:,k)),'-*','Color',get(h1, 'color'),'MarkerIndices',5:40:length(P_dB)-5)
    plot(P_dB,10*log10(SINR_USW_MMSE(:,k)),'-+','Color',get(h1, 'color'),'MarkerIndices',15:40:length(P_dB)-15)

    plot(P_dB,10*log10(SINR_UPW_MRC(:,k)),'--^','Color',get(h1, 'color'),'MarkerIndices',1:40:length(P_dB))    
    plot(P_dB,10*log10(SINR_UPW_ZF(:,k)),'--*','Color',get(h1, 'color'),'MarkerIndices',5:40:length(P_dB)-5)    
    plot(P_dB,10*log10(SINR_UPW_MMSE(:,k)),'--+','Color',get(h1, 'color'),'MarkerIndices',15:40:length(P_dB)-15)    
end

figure(F_UsersLocations);
text(50,1,{'Antenna Array'},'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',11)
% Antenna elements
for m = 1:length(y)
    plot(w(1,m), w(2,m),'k*')
end

axis equal
grid on

% Total physical size of the modular XL-ULA
dim = [.29 .11 .06 .07];
str = {'Total physical size of the',['modular XL-ULA: ', num2str(D), ' m']};
annotation('textbox',dim,'String',str,'verticalalignment','Bottom','HorizontalAlignment','center','FitBoxToText','on','Color','k','EdgeColor','n','interpreter','latex',...
    'FontWeight','bold','FontSize',11,'BackgroundColor','n');




figure(F_SINR);
axis([-inf inf -inf inf])
grid on
xlabel('$\bar{P_r}$ (dB)')
ylabel('SINR (dB)')
legend('User 1 (MRC-USW)','User 1 (ZF-USW)','User 1 (MMSE-USW)','User 1 (MRC-UPW)','User 1 (ZF-UPW)','User 1 (MMSE-UPW)',...
       'User 2 (MRC-USW)','User 2 (ZF-USW)','User 2 (MMSE-USW)','User 2 (MRC-UPW)','User 2 (ZF-UPW)','User 2 (MMSE-UPW)',...
       'User 3 (MRC-USW)','User 3 (ZF-USW)','User 3 (MMSE-USW)','User 3 (MRC-UPW)','User 3 (ZF-UPW)','User 3 (MMSE-UPW)',...
       'Location','southoutside','numcolumns',3)







F_SumRate = figure;
set(F_SumRate, 'Position',  [1250 580, 560, 420])
set(F_SumRate, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(F_SumRate, 'defaultLegendInterpreter','latex');
set(F_SumRate, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(F_SumRate, 'defaultLineLineWidth',1.5);
set(F_SumRate, 'color','w');
plot(P_dB, SumRate_USW_MRC,'-x','MarkerIndices',1:20:length(P_dB))
hold on
plot(P_dB, SumRate_USW_ZF,'-^','MarkerIndices',1:20:length(P_dB))
plot(P_dB, SumRate_USW_MMSE,'-+','MarkerIndices',10:20:length(P_dB)-10)

plot(P_dB, SumRate_UPW_MRC,'-x','MarkerIndices',1:20:length(P_dB))
plot(P_dB, SumRate_UPW_ZF,'-^','MarkerIndices',1:20:length(P_dB))
plot(P_dB, SumRate_UPW_MMSE,'-+','MarkerIndices',10:20:length(P_dB)-10)
grid on
xlabel('$\bar{P_r}$ (dB)')
ylabel('Sum Rate, R (bps/Hz)')
legend('MRC-USW','ZF-USW','MMSE-USW',...
    'MRC-UPW','ZF-UPW','MMSE-UPW','Location','northwest')


