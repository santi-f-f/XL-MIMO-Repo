% NEAR-FIELD BEAM FOCUSING PATTERNS AND GRATING LOBES FOR MODULAR XL-ULA
% Near-Field Observation with Near-Field Beamforming

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

%% User location
% r = 1 + (3-1)* rand(1,1);                             % Distance from the array center, r ∈ [1, 3] m
r = 200;                                                % Making distance of the user to the center of the array > Rayleigh's distance

% theta = -pi/2 + (pi/2 - (-pi/2)).* rand(1,1);         % Angle with respect to the positive x-axis, θ ∈ [-π/2, π/2]
theta = 0;                                              % Making theta = 0 for convenience

q = [r*cos(theta), r*sin(theta)].';

%% Distance between q and the m-th element in module n
% r_n,m = |q - w|

r_nm = zeros(1,length(y));
for i = 1:length(y)
    r_nm(i) = norm(q - w(:,i));
end

%% Rayleigh distance of the whole array and each module
r_Rayleigh_eachModule = 2*(S^2)/lambda;
r_Rayleigh_wholeArray = 2*(D^2)/lambda;

disp(['Given that the user is located ',num2str(r),' meters from the center of the array,'])
disp(['Rayleigh distance of the whole array is ', num2str(r_Rayleigh_wholeArray), ' meters, and '])
disp(['Rayleigh distance of each module is ', num2str(r_Rayleigh_eachModule), ' meters, '])
if r >= r_Rayleigh_wholeArray
    disp('The user/scatterer is located within the far-field region')
elseif (r < r_Rayleigh_wholeArray) && (r >= r_Rayleigh_eachModule)
    disp('The user/scatterer is located within the far-field region of each array module, but within the near-field region of the whole array')
end


%% Near-Field Observation with Near-Field Beamforming
% It is considered the beam focusing pattern at the near-field observation location by using near-field beamforming
%
% When 2S^2/λ ≤ r < 2D^2/λ, the near-field beam focusing pattern for modular XL-ULA in (18)
% depends on the locations (rn, θn) and (r′n , θ′n) viewed from different array modules.



% GNF,NF(r, θ; r′, θ′) --> Beam focusing pattern at the near-field observation location under near-field beamforming

% USW-based near-field array response vector (is fixed)
% Observation of the user based in the Sub-array based USW model for distinct AoAs/AoDs (Eqs. 5, 6, 7)
b_ARV = zeros(length(MM),length(y_n));
r_n0_ARV = zeros(1, length(y_n));
sin_theta_n_ARV = zeros(1, length(y_n));
a_USW_ARV_temp = [];
for in = 1:length(y_n)
    r_n0_ARV(in) = sqrt(r^2 - 2*r*y_n(in)*sin(theta) + y_n(in)^2);                              % Eq. 5

    sin_theta_n_ARV(in) = (r*sin(theta) - y_n(in))/r_n0_ARV(in);                                    % Eq. 6

    for m = 1:length(MM)
        b_ARV(m,in) = exp(1i*2*pi*MM(m)*d*sin_theta_n_ARV(in)/lambda);
    end

    a_USW_ARV_temp = [a_USW_ARV_temp; (exp(-1i*2*pi*r_n0_ARV(in)/lambda)*b_ARV(:,in))];         % eq. 7
end

a_ARV = a_USW_ARV_temp;

% USW-based near-field beamforming (is variable with theta')
% Beamforming based in the Sub-array based USW model for distinct AoAs/AoDs (Eqs. 5, 6, 7)

r_prima = -100:1:1000;                                                             % Variable angle

Delta_r = zeros(1,length(r_prima));
r_n0_BF = zeros(length(y_n),length(r_prima));
sin_theta_n_BF = zeros(length(y_n),length(r_prima));
a_BF = zeros(length(y),length(r_prima));
G_NF_NF = zeros(1,length(r_prima));
for rv = 1:length(r_prima)
    Delta_r(rv) = r - r_prima(rv);
    b_approx_2 = zeros(length(MM),length(y_n));
    a_BF_temp = [];
    for in = 1:length(y_n)
        r_n0_BF(in,rv) = sqrt(r_prima(rv)^2 - 2*r_prima(rv)*y_n(in)*sin(theta) + y_n(in)^2);

        sin_theta_n_BF(in,rv) = (r_prima(rv)*sin(theta) - y_n(in))/r_n0_BF(in,rv);

        for m = 1:length(MM)
            b_approx_2(m,in) = exp(1i*2*pi*MM(m)*d*sin_theta_n_BF(in,rv)/lambda);
        end

        a_BF_temp = [a_BF_temp; (exp(-1i*2*pi*r_n0_BF(in,rv)/lambda)*b_approx_2(:,in))];
    end
    a_BF(:,rv) = a_BF_temp;

    G_NF_NF(rv) = (1/(M*N))*abs(a_BF(:,rv)'*a_ARV);                                             % Eq. 11
end


%% Near-Field Observation with Near-Field Beamforming
%
% Eq. 18
%
% G_NF_NF_18 is Eq. 18

r_n = zeros(1,length(y_n));
sin_theta_n = zeros(1,length(y_n));
for in = 1:length(y_n)
    r_n(in) = sqrt(r^2 - 2*r*y_n(in)*sin(theta) + y_n(in)^2);
    sin_theta_n(in) = (r*sin(theta) - y_n(in))/r_n(in);                                    % Eq. 6
end

r_n_prima = zeros(length(y_n),length(r_prima));
Delta_rn = zeros(length(y_n),length(r_prima));
sin_theta_n_prima = zeros(length(y_n),length(r_prima));
Delta_theta_n = zeros(length(y_n),length(r_prima));
H_M_dbar_DeltaTheta_n = zeros(length(y_n),length(r_prima));
G_NF_NF_18 = zeros(1,length(r_prima));
G_NF_NF_19 = zeros(1,length(r_prima));

Delta_theta = 1;
for rv = 1:length(r_prima)
    %       H_M_dbar(rv) = sin(pi*M*d_bar*Delta_theta)/(M*sin(pi*d_bar*Delta_theta));

    %       θ = θ′ = 0 --> Delta_theta = 0 --> but, H_M_dbar can't be 0 --> for some reason H_M_dbar must be one for every r'
    H_M_dbar(rv) = 1;

    %       H_MN_dbar(rv) = sin(pi*M*N*d_bar*Delta_theta)/(M*N*sin(pi*d_bar*Delta_theta));
    %       G_NF_FF_COL(rv) = abs(H_MN_dbar(rv));
    %       G_NF_FF_COL can't be calculated as before because Delta_theta = 0

    G_NF_NF_18_temp = 0;
    termino_exp_temp = 0;
    for in = 1:length(y_n)
        r_n_prima(in,rv) = sqrt(r_prima(rv)^2 - 2*r_prima(rv)*y_n(in)*sin(theta) + y_n(in)^2);
        Delta_rn(in,rv) =  r_n(in) - r_n_prima(in,rv);

        sin_theta_n_prima(in,rv) = (r_prima(rv)*sin(theta) - y_n(in))/r_n_prima(in,rv);
        Delta_theta_n(in,rv) = sin_theta_n(in) - sin_theta_n_prima(in,rv);

        H_M_dbar_DeltaTheta_n(in,rv) = sin(pi*M*d_bar*Delta_theta_n(in,rv))/(M*sin(pi*d_bar*Delta_theta_n(in,rv)));

        G_NF_NF_18_temp = G_NF_NF_18_temp + exp(-1i*2*pi/lambda*Delta_rn(in,rv))*H_M_dbar_DeltaTheta_n(in,rv);

        termino_exp_temp = termino_exp_temp + exp(-1i*2*pi/lambda*Delta_rn(in,rv));
    end
    G_NF_NF_18(rv) = (1/N)*abs(G_NF_NF_18_temp);
    G_NF_NF_19(rv) = (1/N)*abs(termino_exp_temp)*abs(H_M_dbar(rv));
end


%% Collocated calculation
Gamma_M = M;                         % Modular separation parameter (dependent on the discontinuous surface of practical installation structure) Gamma >= M

% Position of the m-th element within module n
y_Gamma_M = zeros(1,N*M);
indice = 0;
for n = 1:N
    for m = 1:M
        indice = indice + 1;
        position = (NN(n)*Gamma_M + MM(m))*d;
        y_Gamma_M(indice) = position;
    end
end
w_Gamma_M = [zeros(1,length(y_Gamma_M)); y_Gamma_M];

% Position of the m-th element within module n (sub-array based USW model for distinct AoAs/AoDs)
y_n_Gamma_M = zeros(1,length(NN));
for n = 1:length(NN)
    y_n_Gamma_M(n) = NN(n)*Gamma_M*d;
end

% Distance between q and the m-th element in module n
% r_n,m = |q - w|

r_nm_Gamma_M = zeros(1,length(y));
for i = 1:length(y)
    r_nm_Gamma_M(i) = norm(q - w_Gamma_M(:,i));
end

r_n_Gamma_M = zeros(1,length(y_n_Gamma_M));
sin_theta_n_Gamma_M = zeros(1,length(y_n_Gamma_M));
for in = 1:length(y_n_Gamma_M)
    r_n_Gamma_M(in) = sqrt(r^2 - 2*r*y_n_Gamma_M(in)*sin(theta) + y_n_Gamma_M(in)^2);
    sin_theta_n_Gamma_M(in) = (r*sin(theta) - y_n_Gamma_M(in))/r_n_Gamma_M(in);                                    % Eq. 6
end



r_n_prima_Gamma_M = zeros(length(y_n), length(r_prima));
Delta_rn_Gamma_M = zeros(length(y_n), length(r_prima));
sin_theta_n_prima_Gamma_M = zeros(length(y_n), length(r_prima));
Delta_theta_n_Gamma_M = zeros(length(y_n), length(r_prima));
H_M_dbar_DeltaTheta_n_Gamma_M = zeros(length(y_n), length(r_prima));
G_NF_NF_COL = zeros(1, length(r_prima));
for rv = 1:length(r_prima)

    G_NF_NF_COL_temp = 0;
    termino_exp_temp = 0;
    for in = 1:length(y_n)
        r_n_prima_Gamma_M(in,rv) = sqrt(r_prima(rv)^2 - 2*r_prima(rv)*y_n_Gamma_M(in)*sin(theta) + y_n_Gamma_M(in)^2);
        Delta_rn_Gamma_M(in,rv) =  r_n_Gamma_M(in) - r_n_prima_Gamma_M(in,rv);

        sin_theta_n_prima_Gamma_M(in,rv) = (r_prima(rv)*sin(theta) - y_n_Gamma_M(in))/r_n_prima_Gamma_M(in,rv);
        Delta_theta_n_Gamma_M(in,rv) = sin_theta_n(in) - sin_theta_n_prima_Gamma_M(in,rv);

        H_M_dbar_DeltaTheta_n_Gamma_M(in,rv) = sin(pi*M*d_bar*Delta_theta_n_Gamma_M(in,rv))/(M*sin(pi*d_bar*Delta_theta_n_Gamma_M(in,rv)));

        G_NF_NF_COL_temp = G_NF_NF_COL_temp + exp(-1i*2*pi/lambda*Delta_rn_Gamma_M(in,rv))*H_M_dbar_DeltaTheta_n_Gamma_M(in,rv);


    end
    G_NF_NF_COL(rv) = (1/N)*abs(G_NF_NF_COL_temp);

end





%% Figures

% 1st calculation
f_Beam_Patterns_NearField_FarField_1 = figure;
set(f_Beam_Patterns_NearField_FarField_1, 'Position',  [50 580, 560, 420])
set(f_Beam_Patterns_NearField_FarField_1, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_Beam_Patterns_NearField_FarField_1, 'defaultLegendInterpreter','latex');
set(f_Beam_Patterns_NearField_FarField_1, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_Beam_Patterns_NearField_FarField_1, 'defaultLineLineWidth',1.5);
set(f_Beam_Patterns_NearField_FarField_1, 'color','w');
plot(-Delta_r,10*log10(G_NF_NF),'r')
hold on
plot(-Delta_r,10*log10(G_NF_NF_18),'k')
plot(-Delta_r,10*log10(G_NF_NF_19),'b')
plot(-Delta_r,10*log10(G_NF_NF_COL),'g')
plot([0 0],[-10 0],'k--','LineWidth',1)
yline(-3,'k--','LineWidth',1)
axis([-150 800 -10 5])
grid on
xlabel('$\Delta_{r}$ (m)')
ylabel('$G_{NF,NF}(r, \theta; r`, \theta`)$ (dB)')
legend('Modular, (11)','Modular, (18)','Modular, (19)', 'Collocated, $\Gamma$ = M')


% Main Lobe
dim = [.35 .55 .06 .07];
str = {'Main Lobe'};
annotation('textbox',dim,'String',str,'verticalalignment','Bottom','HorizontalAlignment','center','FitBoxToText','on','Color','k','EdgeColor','n','interpreter','latex',...
    'FontWeight','bold','FontSize',14,'BackgroundColor','n');

dim = [.58 .48 .06 .07];
str = {'$BW^{mod}_r = \Delta_{r^+, res}^{mod}(r`, \theta`) + \Delta_{r^-, res}^{mod}(r`, \theta`)$'};
annotation('textbox',dim,'String',str,'verticalalignment','Bottom','HorizontalAlignment','center','FitBoxToText','on','Color','k','EdgeColor','n','interpreter','latex',...
    'FontWeight','bold','FontSize',14,'BackgroundColor','n');

dim = [.08 .45 .06 .07];
str = {'-3'};
annotation('textbox',dim,'String',str,'verticalalignment','Bottom','HorizontalAlignment','center','FitBoxToText','on','Color','k','EdgeColor','n','interpreter','latex',...
    'FontWeight','bold','FontSize',12,'BackgroundColor','n');

xa = [.225 .29];
ya = [.495 .495];
line_Radi = annotation('doublearrow',xa,ya);
line_Radi.Head1Width = 10;
line_Radi.Head2Width = 10;
line_Radi.Head1Length = 10;
line_Radi.Head2Length = 10;

xa = [.225 .225];
ya = [.495 .495];
line_Radi = annotation('doublearrow',xa,ya);
line_Radi.Head1Width = 30;
line_Radi.Head2Width = 30;
line_Radi.Head1Length = 0;
line_Radi.Head2Length = 0;

xa = [.29 .29];
ya = [.495 .495];
line_Radi = annotation('doublearrow',xa,ya);
line_Radi.Head1Width = 30;
line_Radi.Head2Width = 30;
line_Radi.Head1Length = 0;
line_Radi.Head2Length = 0;




%% 1) Spatial Resolution: Equation 21
theta_prima = theta;

C_x = zeros(1,length(r_prima));
S_x = zeros(1,length(r_prima));
F_x = zeros(1,length(r_prima));
G_NF_NF_21 = zeros(1,length(r_prima));

z = zeros(1,length(r_prima));
C_E_x = zeros(1,length(r_prima));
S_E_x = zeros(1,length(r_prima));
F_E_x = zeros(1,length(r_prima));
E_NF_NF = zeros(1,length(r_prima));

for rv = 1:length(r_prima)
    if (1/r - 1/r_prima(rv)) ~= -inf
        x = ((N*Gamma*abs(cos(theta_prima)))/2)*sqrt(pi*d_bar*d*abs(1/r - 1/r_prima(rv)));
    end

    C_x(rv) = integral_C(x);
    S_x(rv) = integral_S(x);

    F_x(rv) = C_x(rv) + 1i*S_x(rv);

    G_NF_NF_21(rv) = abs(F_x(rv)/x);

    z(rv) = 1/r - 1/r_prima(rv);
    if (z(rv)) ~= -inf
        x_E = ((N*Gamma*abs(cos(theta_prima)))/2)*sqrt(pi*d_bar*d*abs(z(rv)));
    end

    C_E_x(rv) = integral_C(x_E);
    S_E_x(rv) = integral_S(x_E);

    F_E_x(rv) = C_E_x(rv) + 1i*S_E_x(rv);

    E_NF_NF(rv) = abs(F_E_x(rv)/x_E);
end

f_Beam_Patterns_NearField_FarField_1 = figure;
set(f_Beam_Patterns_NearField_FarField_1, 'Position',  [680 580, 560, 420])
set(f_Beam_Patterns_NearField_FarField_1, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_Beam_Patterns_NearField_FarField_1, 'defaultLegendInterpreter','latex');
set(f_Beam_Patterns_NearField_FarField_1, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_Beam_Patterns_NearField_FarField_1, 'defaultLineLineWidth',1.5);
set(f_Beam_Patterns_NearField_FarField_1, 'color','w');

plot(-Delta_r,10*log10(G_NF_NF_21),'k')
hold on
plot(-Delta_r,10*log10(E_NF_NF),'b')
axis([-150 800 -10 5])
grid on
xlabel('$\Delta_{r}$ (m)')
ylabel('$G_{NF,NF}(r, \theta; r`, \theta`)$ (dB)')
legend('Eq. 21','Eq. 22')


%% Half power effective distance
r_hp = 0.1*(cos(theta_prima)^2)*2*(D^2)/lambda;

D_approx = Gamma*N*d;
r_hp_approx = 0.1*(cos(theta_prima)^2)*2*(D_approx^2)/lambda;

%% Effective Rayleigh distance
r_eff = 0.367*(cos(theta_prima)^2)*2*(D^2)/lambda;
r_eff_approx = 0.367*(cos(theta_prima)^2)*2*(D_approx^2)/lambda;

disp(['Rayleigh distance of the whole array: ', num2str(r_Rayleigh_wholeArray), ' meters, and '])
disp(['Half power effective distance: ', num2str(r_hp), ' meters, and '])
disp(['Effective Rayleigh distance: ', num2str(r_eff), ' meters, and '])

