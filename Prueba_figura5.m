clc; clear; close all

%% Parameters
c = physconst('LightSpeed');        % Speed of light
fc = 2.4e9;                         % Operating frequency
lambda = c/fc;                      % Signal wavelength
d = lambda/2;                       % Inter-element spacing for antennas within each module
d_bar = d/lambda;

beta_0_dB = 30;
beta_0 = 10.^(beta_0_dB/10);        % beta_0 denoting the channel power at the reference distance d0 = 1 m.

N = 32;                              % Number of modules
M = 4;                              % Number of antenna elements within each module
NM = N*M;                           % Total number of array elements
S = (M - 1)*d;                      % Physical size of each module

Gamma = 13;                         % Modular separation parameter (dependent on the discontinuous surface of practical installation structure) Gamma >= M
Gamma_d = Gamma*d;                  % Inter-module distance between the reference elements (say the center elements of different modules)

D = ((N - 1)*Gamma + (M - 1))*d;    % Total physical size of the modular XL-ULA is

NN = -(N-1)/2:(N-1)/2;
MM = -(M-1)/2:(M-1)/2;

theta_prima = -pi/2:1e-3:pi/2;                                                             % Variable angle

% Position of the m-th element within module n
y = zeros(1,N*M);
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

%% Beamforming pointing to an hypotetical User location
r = 200;                                                % Distance from the array center

theta = 0;                                              % Angle with respect to the positive x-axis, θ ∈ [-π/2, π/2]

q = [r*cos(theta), r*sin(theta)].';

%% Distance between q and the m-th element in module n
% r_n,m = |q - w|

r_nm = zeros(1,length(y));
for i = 1:length(y)
    r_nm(i) = norm(q - w(:,i));
end

%% Generate the beamforming with the observation of the hypotetical user
r_n0_BF = zeros(length(y_n),length(theta_prima));
sin_theta_n_BF = zeros(length(y_n),length(theta_prima));
a_BF = zeros(length(y),length(theta_prima));
for tv = 1:length(theta_prima)
    b_approx_2 = zeros(length(MM),length(y_n));
    a_BF_temp = [];
    for in = 1:length(y_n)
        r_n0_BF(in,tv) = sqrt(r^2 - 2*r*y_n(in)*sin(theta_prima(tv)) + y_n(in)^2);

        sin_theta_n_BF(in,tv) = (r*sin(theta_prima(tv)) - y_n(in))/r_n0_BF(in,tv);

        for m = 1:length(MM)
            b_approx_2(m,in) = exp(1i*2*pi*MM(m)*d*sin_theta_n_BF(in,tv)/lambda);
        end

        a_BF_temp = [a_BF_temp; (exp(-1i*2*pi*r_n0_BF(in,tv)/lambda)*b_approx_2(:,in))];
    end
    a_BF(:,tv) = a_BF_temp;
end




%% User location
% Generate a grid with every possible locations of the user 
% (rectangular coordinates)
y_pos = 0:2:10;
x_pos = 0:2:40;
[X,Y] = meshgrid(x_pos,y_pos);

[Theta_Var,R_Var] = cart2pol(X,Y);  

f_locations = figure;
set(f_locations, 'Position',  [20, 60, 560, 420*2])
set(f_locations, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_locations, 'defaultLegendInterpreter','latex');
set(f_locations, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_locations,'defaultLineLineWidth',2);
set(f_locations,'color','w');
subplot(3,1,1)
plot(X,Y)
xlabel('X (m)')
ylabel('Y (m)')

subplot(3,1,2)
plot(R_Var,Theta_Var)
xlabel('X (m)')
ylabel('Y (m)')

subplot(3,1,3) 
polarplot(Theta_Var,R_Var)

%% Generate the Array response vector for every possible location of the user
r_n0_ARV = zeros(length(y_n),length(y_pos),length(x_pos));
sin_theta_n_ARV = zeros(length(y_n),length(y_pos),length(x_pos));
a_USW_ARV_temp = [];
for in = 1:length(y_n)
    r_n0_ARV(in,:,:) = sqrt(R_Var.^2 - 2*R_Var.*y_n(in).*sin(Theta_Var) + y_n(in)^2);                              % Eq. 5

    sin_theta_n_ARV(in,:,:) = (R_Var.*sin(Theta_Var) - y_n(in))./squeeze(r_n0_ARV(in,:,:));                                    % Eq. 6
end

r_n0_ARV = permute(r_n0_ARV,[1 3 2]) ;
sin_theta_n_ARV = permute(sin_theta_n_ARV,[1 3 2]) ;

b_ARV = zeros(length(MM),length(y_n),length(x_pos),length(y_pos));
for in = 1:length(y_n)
    for m = 1:length(MM)
        b_ARV(m,in,:,:) = exp(1i*2*pi*MM(m)*d.*sin_theta_n_ARV(in,:,:)/lambda);
    end

    a_USW_ARV_temp = [a_USW_ARV_temp; (exp(-1i*2*pi*r_n0_ARV(in,:,:)/lambda).*squeeze(b_ARV(:,in,:,:)))];         % eq. 7
end

a_ARV = a_USW_ARV_temp;


%% Evaluation of the power reached by the user in every possible location
G_NF_NF = zeros(length(theta_prima),length(x_pos),length(y_pos));
for tv = 1:length(theta_prima)
    G_NF_NF(tv,:,:) = (1/(M*N))*abs(tensorprod(a_BF(:,tv)',a_ARV,2,1));
end

Power = squeeze(10*log10(mean(G_NF_NF.*conj(G_NF_NF)))).';

f_power = figure;
set(f_power, 'Position',  [680, 480, 560, 420])
set(f_power, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_power, 'defaultLegendInterpreter','latex');
set(f_power, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_power,'defaultLineLineWidth',2);
set(f_power,'color','w');
contour(X,Y,Power)


f_beamforming = figure;
set(f_beamforming, 'Position',  [1260, 480, 560, 420])
set(f_beamforming, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_beamforming, 'defaultLegendInterpreter','latex');
set(f_beamforming, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_beamforming,'defaultLineLineWidth',2);
set(f_beamforming,'color','w');
hold on
for tv = 1:length(theta_prima) 
    contour(X,Y,10*log10(squeeze(G_NF_NF(tv,:,:))).')
end