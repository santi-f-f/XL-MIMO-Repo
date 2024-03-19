% NEAR-FIELD BEAM FOCUSING PATTERNS AND GRATING LOBES FOR MODULAR XL-ULA
% Near-Field Observation with Far-Field Beamforming

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

%% Near-Field Observation with Far-Field Beamforming

theta_prima = -pi/2:1e-3:pi/2;                                                             % Variable angle

% USW-based near-field array response vector (is fixed)
% Observation of the user based in the Sub-array based USW model for distinct AoAs/AoDs (Eqs. 5, 6, 7)
b_ARV = zeros(length(MM),length(y_n));
r_n0_ARV = zeros(1, length(y_n));
sin_theta_n = zeros(1, length(y_n));
a_USW_ARV_temp = [];
for in = 1:length(y_n)
    r_n0_ARV(in) = sqrt(r^2 - 2*r*y_n(in)*sin(theta) + y_n(in)^2);                              % Eq. 5

    sin_theta_n(in) = (r*sin(theta) - y_n(in))/r_n0_ARV(in);                                    % Eq. 6

    for m = 1:length(MM)
        b_ARV(m,in) = exp(1i*2*pi*MM(m)*d*sin_theta_n(in)/lambda);
    end

    a_USW_ARV_temp = [a_USW_ARV_temp; (exp(-1i*2*pi*r_n0_ARV(in)/lambda)*b_ARV(:,in))];         % eq. 7
end

a_ARV = a_USW_ARV_temp;

% UPW-based far-field beamforming (is variable with theta')
Delta_theta = zeros(1,length(theta_prima));
r_nm_BF = zeros(length(y),length(theta_prima));
p_BF = zeros(length(NN),length(theta_prima));
b_BF = zeros(length(MM),length(theta_prima));
a_BF = zeros(length(y),length(theta_prima));
G_UPW_NF_FF = zeros(1,length(theta_prima));
for tv = 1:length(theta_prima)  
    Delta_theta(tv) = sin(theta)-sin(theta_prima(tv));

    for i = 1:length(y)
        r_nm_BF(i,tv) = r - y(i)*sin(theta_prima(tv));                                           % Unnecessary
    end

    for n = 1:length(NN)
        p_BF(n,tv) = exp(1i*2*pi*NN(n)*Gamma*d*sin(theta_prima(tv))/lambda);
    end

    for m = 1:length(MM)
        b_BF(m,tv) = exp(1i*2*pi*MM(m)*d*sin(theta_prima(tv))/lambda);
    end

    a_BF(:,tv) = exp(-1i*2*pi.*r/lambda)*kron(p_BF(:,tv),b_BF(:,tv));                            % Eq. 4

    v_BF = a_BF;

    G_UPW_NF_FF(tv) = (1/(M*N))*abs(v_BF(:,tv)'*a_ARV);                                          % Eq. 10 - Far-field beamforming with Near-field observation 
end



%% Near-Field Observation with Far-Field Beamforming
% Eq. 16 & 17
%
% When the observation location (r, θ) is within the Rayleigh distance of the modular XL-ULA,
% we use the two near-field array models to simplify the beam focusing pattern.
%
% G_NF_FF_16 is Eq. 16
% G_NF_FF_17 is Eq. 17

H_MN_dbar = zeros(1,length(theta_prima));
H_M_dbar = zeros(1,length(theta_prima));
G_NF_FF_COL = zeros(1,length(theta_prima));
G_NF_FF_16 = zeros(1,length(theta_prima));
G_NF_FF_17 = zeros(1,length(theta_prima));
for tv = 1:length(theta_prima)
    H_MN_dbar(tv) = sin(pi*M*N*d_bar*Delta_theta(tv))/(M*N*sin(pi*d_bar*Delta_theta(tv)));
    H_M_dbar(tv) = sin(pi*M*d_bar*Delta_theta(tv))/(M*sin(pi*d_bar*Delta_theta(tv)));

    G_NF_FF_COL(tv) = abs(H_MN_dbar(tv));

    G_NF_FF_temp_16 = 0;
    G_NF_FF_temp_17 = 0;
    for in = 1:length(y_n)
        H_M_dbar_DeltaTheta_n = sin(pi*M*d_bar*(sin_theta_n(in) - sin(theta_prima(tv))))/(M*sin(pi*d_bar*(sin_theta_n(in) - sin(theta_prima(tv)))));

        G_NF_FF_temp_16 = G_NF_FF_temp_16 + exp(-1i*2*pi/lambda*(r_n0_ARV(in) + NN(in)*Gamma*d*(sin(theta_prima(tv)))))*H_M_dbar_DeltaTheta_n;

        G_NF_FF_temp_17 = G_NF_FF_temp_17 + exp(-1i*2*pi/lambda*(r_n0_ARV(in) + NN(in)*Gamma*d*(sin(theta_prima(tv)))));
    end
    G_NF_FF_16(tv) = (1/N)*abs(G_NF_FF_temp_16);
    G_NF_FF_17(tv) = (1/N)*abs(G_NF_FF_temp_17)*abs(H_M_dbar(tv));    
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
plot(Delta_theta,10*log10(G_UPW_NF_FF),'r')
hold on
plot(Delta_theta,10*log10(G_NF_FF_16),'k')
axis([-1 1 -40 10])
grid on
xlabel('$\Delta_{\theta}$')
ylabel('$G_{NF,FF}(r, \theta; \theta`)$ (dB)')
legend('Modular, (10)','Modular, (16)')




f_Beam_Patterns_NearField_FarField_subplot = figure;
set(f_Beam_Patterns_NearField_FarField_subplot, 'Position',  [50 60, 560, 420])
set(f_Beam_Patterns_NearField_FarField_subplot, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_Beam_Patterns_NearField_FarField_subplot, 'defaultLegendInterpreter','latex');
set(f_Beam_Patterns_NearField_FarField_subplot, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_Beam_Patterns_NearField_FarField_subplot, 'defaultLineLineWidth',1.5);
set(f_Beam_Patterns_NearField_FarField_subplot, 'color','w');
subplot(2,1,1)
plot(Delta_theta,10*log10(G_UPW_NF_FF),'r')
axis([-1 1 -40 10])
grid on
xlabel('$\Delta_{\theta}$')
ylabel('$G_{NF,FF}(r, \theta; \theta`)$ (dB)')
legend('Modular, (10)')

subplot(2,1,2)
plot(Delta_theta,10*log10(G_NF_FF_16),'k')
axis([-1 1 -40 10])
grid on
xlabel('$\Delta_{\theta}$')
ylabel('$G_{NF,FF}(r, \theta; \theta`)$ (dB)')
legend('Modular, (16)')


% 2nd calculation
f_Beam_Patterns_NearField_FarField = figure;
set(f_Beam_Patterns_NearField_FarField, 'Position',  [640 580, 560, 420])
set(f_Beam_Patterns_NearField_FarField, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_Beam_Patterns_NearField_FarField, 'defaultLegendInterpreter','latex');
set(f_Beam_Patterns_NearField_FarField, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_Beam_Patterns_NearField_FarField, 'defaultLineLineWidth',1.5);
set(f_Beam_Patterns_NearField_FarField, 'color','w');
plot(Delta_theta,10*log10(G_NF_FF_16),'k')
hold on
plot(Delta_theta,10*log10((G_NF_FF_17)),'b')
plot(Delta_theta,10*log10(G_NF_FF_COL),'g')
axis([-1 1 -40 10])
grid on
xlabel('$\Delta_{\theta}$')
ylabel('$G_{NF,FF}(r, \theta; \theta`)$ (dB)')
legend('Modular (16)','Modular (17)', 'Collocated, $\Gamma$ = M')




% 1st & 2nd calculations
f_Beam_Patterns_NearField_FarField_1 = figure;
set(f_Beam_Patterns_NearField_FarField_1, 'Position',  [640 60, 560, 420])
set(f_Beam_Patterns_NearField_FarField_1, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_Beam_Patterns_NearField_FarField_1, 'defaultLegendInterpreter','latex');
set(f_Beam_Patterns_NearField_FarField_1, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_Beam_Patterns_NearField_FarField_1, 'defaultLineLineWidth',1.5);
set(f_Beam_Patterns_NearField_FarField_1, 'color','w');
plot(Delta_theta,10*log10(G_UPW_NF_FF),'r')
hold on
plot(Delta_theta,10*log10(G_NF_FF_16),'k')
plot(Delta_theta,10*log10(G_NF_FF_17),'b')
plot(Delta_theta,10*log10(G_NF_FF_COL),'g')
axis([-1 1 -40 10])
grid on
xlabel('$\Delta_{\theta}$')
ylabel('$G_{NF,FF}(r, \theta; \theta`)$ (dB)')
legend('Modular, (10)','Modular, (16)','Modular (17)','Collocated, $\Gamma$ = M')

% Energy Spread Effect
dim = [.35 .7 .06 .07];
str = {'Energy Spread','Effect'};
annotation('textbox',dim,'String',str,'verticalalignment','Bottom','HorizontalAlignment','center','FitBoxToText','on','Color','k','EdgeColor','n','interpreter','latex',...
    'FontWeight','bold','FontSize',13,'BackgroundColor','n');

dim = [.37 .58 .06 .07];
str = {''};
annotation('textbox',dim,'String',str,'verticalalignment','Bottom','FitBoxToText','on','Color','k','EdgeColor','k','interpreter','latex',...
    'FontWeight','bold','FontSize',13,'BackgroundColor','n');

dim = [.49 .62 .06 .07];
str = {''};
annotation('textbox',dim,'String',str,'verticalalignment','Bottom','FitBoxToText','on','Color','k','EdgeColor','k','interpreter','latex',...
    'FontWeight','bold','FontSize',13,'BackgroundColor','n');

xa = [.43 .485];
ya = [.72 .69];
line_a = annotation('arrow',xa,ya);
line_a.LineWidth = 0.5;

xb = [.4 .41];
yb = [.71 .66];
line_b = annotation('arrow',xb,yb);
line_b.LineWidth = 0.5;

% Grating Lobes
dim = [.7 .65 .06 .07];
str = {'Grating Lobes'};
annotation('textbox',dim,'String',str,'verticalalignment','Bottom','HorizontalAlignment','center','FitBoxToText','on','Color','k','EdgeColor','n','interpreter','latex',...
    'FontWeight','bold','FontSize',13,'BackgroundColor','n');

xc = [.69 .66];
yc = [.67 .62];
line_c = annotation('arrow',xc,yc);
line_c.LineWidth = 0.5;

xd = [.73 .79];
yd = [.67 .58];
line_d = annotation('arrow',xd,yd);
line_d.LineWidth = 0.5;




% Polar coordinates
f_Beam_Patterns_NearField_FarField_Spherical = figure;
set(f_Beam_Patterns_NearField_FarField_Spherical, 'Position',  [1220 60, 660, 520])
set(f_Beam_Patterns_NearField_FarField_Spherical, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_Beam_Patterns_NearField_FarField_Spherical, 'defaultLegendInterpreter','latex');
set(f_Beam_Patterns_NearField_FarField_Spherical, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_Beam_Patterns_NearField_FarField_Spherical, 'defaultLineLineWidth',1.5);
set(f_Beam_Patterns_NearField_FarField_Spherical, 'color','w');

P = polarpattern(Delta_theta*90,(G_UPW_NF_FF),'AngleDirection','cw','AngleAtTop',-90,'AngleResolution',30);
add(P,Delta_theta*90,(G_NF_FF_16))
add(P,Delta_theta*90,(G_NF_FF_17))
% add(P,Delta_theta*90,(G_NF_FF_COL))































%% Far Field Channel (NUSW-based channel vector)
h_farField = sqrt(beta_0)/r.*a_ARV;

F = a_BF;

h_A = F\h_farField;

f_pruebas = figure;
set(f_pruebas, 'Position',  [1220 580, 560, 420])
set(f_pruebas, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_pruebas, 'defaultLegendInterpreter','latex');
set(f_pruebas, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_pruebas, 'defaultLineLineWidth',1.5);
set(f_pruebas, 'color','w');
plot(h_A)

theta_n = asin(sin_theta_n);
for in = 1:length(y_n)
a_theta_l(in) = (1/length(y_n))*exp(1i*pi*theta_n(in));

theta_n_new(in) = (2*NN(in)-N+1)/N;
a_theta_n(in) = (1/length(y_n))*exp(1i*pi*theta_n_new(in));
FF(in) = a_theta_n(in);
end

a_theta_l = a_theta_l.';

h_farfiel_m = 0;
for in =1:length(y_n)
       h_farfiel_m = h_farfiel_m + exp(-1i*2*pi/lambda*r_n0_ARV)*a_theta_l(in); 
end

h_A_new = FF\h_farfiel_m;

f_pruebas = figure;
set(f_pruebas, 'Position',  [1220 580, 560, 420])
set(f_pruebas, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_pruebas, 'defaultLegendInterpreter','latex');
set(f_pruebas, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_pruebas, 'defaultLineLineWidth',1.5);
set(f_pruebas, 'color','w');
plot(h_A_new)