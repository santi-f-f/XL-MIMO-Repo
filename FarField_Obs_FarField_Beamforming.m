% NEAR-FIELD BEAM FOCUSING PATTERNS AND GRATING LOBES FOR MODULAR XL-ULA
% Far-Field Observation with Far-Field Beamforming

clc; clear; close all

%% Parameters
c = physconst('LightSpeed');        % Speed of light
fc = 2.4e9;                         % Operating frequency
lambda = c/fc;                      % Signal wavelength
d = lambda/2;                       % Inter-element spacing for antennas within each module
d_bar = d/lambda;

beta_0_dB = 30;
beta_0 = 10.^(beta_0_dB/10);        % beta_0 denoting the channel power at the reference distance d0 = 1 m.

N = 4;                              % Number of modules
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

%% Far-Field Observation with Far-Field Beamforming
%
% Under far-field beamforming, the beamforming vector v only needs to match the desired direction θ′
% based on the UPW model in (4), i.e., vFF = a(θ′).
% When the observation location (r, θ) is within the far-field region of the array,
% the near-field array response vector a(r, θ) in (10) degenerates to (4).
%
% The array response vector and the beamforming vector are designed both based in the UPW model.

theta_prima = -pi/2:0.001:pi/2;                                                 % Variable angle

% UPW-based far-field array response vector (is fixed)
r_nm_ARV = zeros(1,length(y));
p_ARV = zeros(1,length(NN));
b_ARV = zeros(1,length(MM));
for i = 1:length(y)
    r_nm_ARV(i) = r - y(i)*sin(theta);                                          % Unnecessary
end

for n = 1:length(NN)
    p_ARV(n) = exp(1i*2*pi*NN(n)*Gamma*d*sin(theta)/lambda);
end

for m = 1:length(MM)
    b_ARV(m) = exp(1i*2*pi*MM(m)*d*sin(theta)/lambda);
end

a_ARV = exp(-1i*2*pi.*r/lambda)*kron(p_ARV,b_ARV).';                            % Eq. 4

% UPW-based far-field beamforming (is variable with theta')
Delta_theta = zeros(1,length(theta_prima));
G_UPW_FF_FF = zeros(1,length(theta_prima));
r_nm_BF = zeros(length(y),length(theta_prima));
p_BF = zeros(length(NN),length(theta_prima));
b_BF = zeros(length(MM),length(theta_prima));
a_BF = zeros(M*N,length(theta_prima));
for tv = 1:length(theta_prima)
    Delta_theta(tv) = sin(theta)-sin(theta_prima(tv));

    for i = 1:length(y)
        r_nm_BF(i,tv) = r - y(i)*sin(theta_prima(tv));                          % Unnecessary
    end

    for n = 1:length(NN)
        p_BF(n,tv) = exp(1i*2*pi*NN(n)*Gamma*d*sin(theta_prima(tv))/lambda);
    end

    for m = 1:length(MM)
        b_BF(m,tv) = exp(1i*2*pi*MM(m)*d*sin(theta_prima(tv))/lambda);
    end

    a_BF(:,tv) = exp(-1i*2*pi.*r/lambda)*kron(p_BF(:,tv),b_BF(:,tv));           % Eq. 4

    v_UPW_BF = a_BF;

    G_UPW_FF_FF(tv) = (1/(M*N))*abs(v_UPW_BF(:,tv)'*a_ARV);                     % Eq. 10
end




%% A. Far-Field Observation with Far-Field Beamforming
% When the observation location (r, θ) is within the far-field region of the array,
% the near-field array response vector a(r, θ) degenerates to a(θ).
%
% Eqs. 12 & 13

G_FF_FF_deltaTheta_MOD = zeros(1,length(Delta_theta));
H = zeros(1,length(Delta_theta));
G_FF_FF_deltaTheta_COL = zeros(1,length(Delta_theta));
for Dt = 1:length(Delta_theta)
    G_FF_FF_deltaTheta_MOD(Dt) = abs(sin(pi*N*Gamma*d_bar*Delta_theta(Dt))/(N*sin(pi*Gamma*d_bar*Delta_theta(Dt))))*...
        abs(sin(pi*M*d_bar*Delta_theta(Dt))/(M*sin(pi*d_bar*Delta_theta(Dt))));

    H(Dt) = abs(sin(pi*M*d_bar*Delta_theta(Dt))/(M*sin(pi*d_bar*Delta_theta(Dt))));

    G_FF_FF_deltaTheta_COL(Dt) = abs(sin(pi*M*N*d_bar*Delta_theta(Dt))/(M*N*sin(pi*d_bar*Delta_theta(Dt))));
end


%% Spatial angular resolution of the modular ULA
Delta_theta_res_mod = 1/(N*Gamma*d_bar);

%% Figures

% 1st calculation
f_Beam_Patterns_FarField = figure;
set(f_Beam_Patterns_FarField, 'Position',  [20 580, 560, 420])
set(f_Beam_Patterns_FarField, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_Beam_Patterns_FarField, 'defaultLegendInterpreter','latex');
set(f_Beam_Patterns_FarField, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_Beam_Patterns_FarField, 'defaultLineLineWidth',1.5);
set(f_Beam_Patterns_FarField, 'color','w');

plot(Delta_theta,10*log10(G_UPW_FF_FF),'m')
hold on
plot(Delta_theta,10*log10(G_FF_FF_deltaTheta_MOD),'b--')
axis([-1 1 -40 10])
grid on
legend('$G_{\rm NF,FF}(r, \theta; \theta`) = \frac{1}{MN} |{\rm \boldmath{a}}(\theta`)^H {\rm \boldmath{a}}(r,\theta)|$ (Eq. 10)','Modular (Eq. 12)','FontSize',14)
xlabel('$\Delta_{\theta}$')
ylabel('$G_{FF,FF}(\theta; \theta `)$ (dB)')









% 2nd calculation
f_Beam_Patterns_FarField = figure;
set(f_Beam_Patterns_FarField, 'Position',  [600 580, 560, 420])
set(f_Beam_Patterns_FarField, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_Beam_Patterns_FarField, 'defaultLegendInterpreter','latex');
set(f_Beam_Patterns_FarField, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_Beam_Patterns_FarField, 'defaultLineLineWidth',1.5);
set(f_Beam_Patterns_FarField, 'color','w');
plot(Delta_theta, 10*log10(H),'r--')
hold on
plot(Delta_theta,10*log10(G_FF_FF_deltaTheta_MOD),'b')
plot(Delta_theta,10*log10(G_FF_FF_deltaTheta_COL),'g')
axis([-1 1 -40 10])
grid on
legend('$|H_{M,\bar{d}}(\Delta_{\theta})|$ in (12)','Modular: (12)','Collocated: (13)','FontSize',12)
xlabel('$\Delta_{\theta}$')
ylabel('$G_{FF,FF}(\theta; \theta `)$ (dB)')

% Annotations

% Main lobe
dim = [.47 .2 .2 .2];
str = {'$\frac{2}{NM\bar{d}}$'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','k','EdgeColor','n','interpreter','latex',...
    'FontWeight','bold','FontSize',18,'BackgroundColor','n');

dim = [.445 .09 .2 .2];
str = {'Main lobe'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','k','EdgeColor','w','interpreter','latex',...
    'FontWeight','bold','FontSize',14,'BackgroundColor','n');

xa = [.468 .566];
ya = [.43 .43];
line_Radi = annotation('doublearrow',xa,ya);
line_Radi.LineWidth = 0.5;
line_Radi.Head1Width = 10;
line_Radi.Head2Width = 10;
line_Radi.Head1Length = 7.5;
line_Radi.Head2Length = 7.5;


% Grating lobes
dim = [.18 .58 .2 .2];
str = {'Grating Lobes'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','k','EdgeColor','w','interpreter','latex',...
    'FontWeight','bold','FontSize',13,'BackgroundColor','n');

xa = [.28 .23];
ya = [.72 .65];
line_Radi = annotation('arrow',xa,ya);
line_Radi.LineWidth = 0.5;

xa = [.33 .388];
ya = [.72 .68];
line_Radi = annotation('arrow',xa,ya);
line_Radi.LineWidth = 0.5;


% 1/Gamma d
dim = [.4 .65 .2 .2];
str = {'$\frac{1}{\Gamma\bar{d}}$'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','k','EdgeColor','w','interpreter','latex',...
    'FontWeight','bold','FontSize',15,'BackgroundColor','n');

xa = [.4 .458];
ya = [.74 .74];
line_Radi = annotation('doublearrow',xa,ya);
line_Radi.LineWidth = 0.5;
line_Radi.Head1Width = 30;
line_Radi.Head2Width = 30;
line_Radi.Head1Length = 0;
line_Radi.Head2Length = 0;

xa = [.4 .458];
ya = [.74 .74];
line_Radi = annotation('doublearrow',xa,ya);
line_Radi.LineWidth = 0.5;
line_Radi.Head1Width = 10;
line_Radi.Head2Width = 10;
line_Radi.Head1Length = 7.5;
line_Radi.Head2Length = 7.5;


% 2/N Gamma d
dim = [.48 .65 .2 .2];
str = {'$\frac{2}{N\Gamma\bar{d}}$'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','k','EdgeColor','w','interpreter','latex',...
    'FontWeight','bold','FontSize',15,'BackgroundColor','n');

xa = [.505 .505];
ya = [.74 .74];
line_Radi = annotation('doublearrow',xa,ya);
line_Radi.Head1Width = 30;
line_Radi.Head2Width = 30;
line_Radi.Head1Length = 0;
line_Radi.Head2Length = 0;

xa = [.53 .53];
ya = [.74 .74];
line_Radi = annotation('doublearrow',xa,ya);
line_Radi.Head1Width = 30;
line_Radi.Head2Width = 30;
line_Radi.Head1Length = 0;
line_Radi.Head2Length = 0;

xa = [.47 .505];
ya = [.74 .74];
line_Radi = annotation('arrow',xa,ya);
line_Radi.LineWidth = 0.5;

xa = [.565 .53];
ya = [.74 .74];
line_Radi = annotation('arrow',xa,ya);
line_Radi.LineWidth = 0.5;













f_Beam_Patterns_FarField = figure;
set(f_Beam_Patterns_FarField, 'Position',  [600 50, 560, 420])
set(f_Beam_Patterns_FarField, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_Beam_Patterns_FarField, 'defaultLegendInterpreter','latex');
set(f_Beam_Patterns_FarField, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_Beam_Patterns_FarField, 'defaultLineLineWidth',1.5);
set(f_Beam_Patterns_FarField, 'color','w');
plot(Delta_theta, (H),'r--')
hold on
plot(Delta_theta,(G_FF_FF_deltaTheta_MOD),'b')
plot(Delta_theta,(G_FF_FF_deltaTheta_COL),'g')
% axis([-1 1 -40 10])
grid on
legend('$|H_{M,\bar{d}}(\Delta_{\theta})|$ in (12)','Modular: (12)','Collocated: (13)','FontSize',12)
xlabel('$\Delta_{\theta}$')
ylabel('$G_{FF,FF}(\theta; \theta `)$')
    

f_Beam_Patterns_NearField_FarField_Spherical = figure;
set(f_Beam_Patterns_NearField_FarField_Spherical, 'Position',  [1220 60, 660, 520])
set(f_Beam_Patterns_NearField_FarField_Spherical, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_Beam_Patterns_NearField_FarField_Spherical, 'defaultLegendInterpreter','latex');
set(f_Beam_Patterns_NearField_FarField_Spherical, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_Beam_Patterns_NearField_FarField_Spherical, 'defaultLineLineWidth',1.5);
set(f_Beam_Patterns_NearField_FarField_Spherical, 'color','w');

P = polarpattern(Delta_theta*90,(H),'AngleDirection','cw','AngleAtTop',-90,'AngleResolution',30);
add(P,Delta_theta*90,(G_FF_FF_deltaTheta_MOD))
add(P,Delta_theta*90,(G_FF_FF_deltaTheta_COL))
