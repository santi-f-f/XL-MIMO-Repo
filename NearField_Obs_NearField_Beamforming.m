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

theta_prima = -pi/2:1e-3:pi/2;                                                             % Variable angle

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
Delta_theta = zeros(1,length(theta_prima));
r_n0_BF = zeros(length(y_n),length(theta_prima));
sin_theta_n_BF = zeros(length(y_n),length(theta_prima));
a_BF = zeros(length(y),length(theta_prima));
G_NF_NF = zeros(1,length(theta_prima));
for tv = 1:length(theta_prima)
    Delta_theta(tv) = sin(theta)-sin(theta_prima(tv));
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

    G_NF_NF(tv) = (1/(M*N))*abs(a_BF(:,tv)'*a_ARV);                                             % Eq. 11
end





%% Near-Field Observation with Near-Field Beamforming
%
% Eq. 18
%
% G_NF_NF_18 is Eq. 18
r_prima = r;

r_n = zeros(1,length(y_n));
sin_theta_n = zeros(1,length(y_n));
for in = 1:length(y_n)
    r_n(in) = sqrt(r^2 - 2*r*y_n(in)*sin(theta) + y_n(in)^2);
    sin_theta_n(in) = (r*sin(theta) - y_n(in))/r_n(in);                                    % Eq. 6
end

H_M_dbar = zeros(1,length(theta_prima));
H_MN_dbar = zeros(1,length(theta_prima));
G_NF_FF_COL = zeros(1,length(theta_prima));
r_n_prima = zeros(length(y_n),length(theta_prima));
Delta_rn = zeros(length(y_n),length(theta_prima));
sin_theta_n_prima = zeros(length(y_n),length(theta_prima));
Delta_theta_n = zeros(length(y_n),length(theta_prima));
H_M_dbar_DeltaTheta_n = zeros(length(y_n),length(theta_prima));
G_NF_NF_18 = zeros(1,length(theta_prima));
G_NF_NF_19 = zeros(1,length(theta_prima));
for tv = 1:length(theta_prima)
    H_M_dbar(tv) = sin(pi*M*d_bar*Delta_theta(tv))/(M*sin(pi*d_bar*Delta_theta(tv)));

    H_MN_dbar(tv) = sin(pi*M*N*d_bar*Delta_theta(tv))/(M*N*sin(pi*d_bar*Delta_theta(tv)));
    G_NF_FF_COL(tv) = abs(H_MN_dbar(tv));


    G_NF_NF_18_temp = 0;
    termino_exp_temp = 0;
    for in = 1:length(y_n)
        r_n_prima(in,tv) = sqrt(r_prima^2 - 2*r_prima*y_n(in)*sin(theta_prima(tv)) + y_n(in)^2);
        Delta_rn(in,tv) =  r_n(in) - r_n_prima(in,tv);

        sin_theta_n_prima(in,tv) = (r_prima*sin(theta_prima(tv)) - y_n(in))/r_n_prima(in,tv);
        Delta_theta_n(in,tv) = sin_theta_n(in) - sin_theta_n_prima(in,tv);

        H_M_dbar_DeltaTheta_n(in,tv) = sin(pi*M*d_bar*Delta_theta_n(in,tv))/(M*sin(pi*d_bar*Delta_theta_n(in,tv)));

        G_NF_NF_18_temp = G_NF_NF_18_temp + exp(-1i*2*pi/lambda*Delta_rn(in,tv))*H_M_dbar_DeltaTheta_n(in,tv);

        termino_exp_temp = termino_exp_temp + exp(-1i*2*pi/lambda*Delta_rn(in,tv));
    end
    G_NF_NF_18(tv) = (1/N)*abs(G_NF_NF_18_temp);
    G_NF_NF_19(tv) = (1/N)*abs(termino_exp_temp)*abs(H_M_dbar(tv));

end



%% Figures

f_Beam_Patterns_NearField_FarField_1 = figure;
set(f_Beam_Patterns_NearField_FarField_1, 'Position',  [30 580, 560, 420])
set(f_Beam_Patterns_NearField_FarField_1, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_Beam_Patterns_NearField_FarField_1, 'defaultLegendInterpreter','latex');
set(f_Beam_Patterns_NearField_FarField_1, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_Beam_Patterns_NearField_FarField_1, 'defaultLineLineWidth',1.5);
set(f_Beam_Patterns_NearField_FarField_1, 'color','w');
plot(Delta_theta,10*log10(G_NF_NF),'r')
hold on
plot(Delta_theta,10*log10(G_NF_NF_18),'k')
plot(Delta_theta,10*log10(G_NF_NF_19),'b')
plot(Delta_theta,10*log10(G_NF_FF_COL),'g')
axis([-1 1 -40 10])
grid on
xlabel('$\Delta_{\theta}$')
ylabel('$G_{NF,NF}(r, \theta; r`, \theta`)$ (dB)')
legend('Modular, (11)','Modular, (18)','Modular, (19)', 'Collocated, $\Gamma$ = M')

% Annotations
% Energy Spread Effect
dim = [.745 .65 .06 .07];
str = {'Energy Spread Effect'};
annotation('textbox',dim,'String',str,'verticalalignment','Bottom','HorizontalAlignment','center','FitBoxToText','on','Color','k','EdgeColor','n','interpreter','latex',...
    'FontWeight','bold','FontSize',11.5,'BackgroundColor','n');

dim = [.787 .545 .06 .07];
str = {''};
annotation('textbox',dim,'String',str,'verticalalignment','Bottom','FitBoxToText','on','Color','k','EdgeColor','k','interpreter','latex',...
    'FontWeight','bold','FontSize',13,'BackgroundColor','n');

xa = [.77 .8];
ya = [.67 .63];
line_a = annotation('arrow',xa,ya);
line_a.LineWidth = 0.5;

% Grating Lobes
dim = [.35 .57 .06 .07];
str = {'Grating Lobes'};
annotation('textbox',dim,'String',str,'verticalalignment','Bottom','HorizontalAlignment','center','FitBoxToText','on','Color','k','EdgeColor','n','interpreter','latex',...
    'FontWeight','bold','FontSize',11.5,'BackgroundColor','n');

xa = [.33 .39];
ya = [.63 .695];
line_a = annotation('arrow',xa,ya);
line_a.LineWidth = 0.5;

xa = [.42 .45];
ya = [.63 .73];
line_a = annotation('arrow',xa,ya);
line_a.LineWidth = 0.5;

% 1/Gamma d
dim = [.22 .6 .06 .07];
str = {'$\frac{1}{\Gamma\bar{d}}$'};
annotation('textbox',dim,'String',str,'verticalalignment','Bottom','HorizontalAlignment','center','FitBoxToText','on','Color','k','EdgeColor','n','interpreter','latex',...
    'FontWeight','bold','FontSize',14,'BackgroundColor','n');

xa = [.218 .218];
ya = [.595 .595];
line_Radi = annotation('doublearrow',xa,ya);
line_Radi.Head1Width = 30;
line_Radi.Head2Width = 30;
line_Radi.Head1Length = 0;
line_Radi.Head2Length = 0;

xa = [.28 .28];
ya = [.595 .595];
line_Radi = annotation('doublearrow',xa,ya);
line_Radi.Head1Width = 30;
line_Radi.Head2Width = 30;
line_Radi.Head1Length = 0;
line_Radi.Head2Length = 0;

xa = [.218 .28];
ya = [.595 .595];
line_Radi = annotation('doublearrow',xa,ya);
line_Radi.Head1Width = 7.5;
line_Radi.Head2Width = 7.5;

% Main Lobe
dim = [.43 .78 .06 .07];
str = {'Main Lobe'};
annotation('textbox',dim,'String',str,'verticalalignment','Bottom','HorizontalAlignment','center','FitBoxToText','on','Color','k','EdgeColor','n','interpreter','latex',...
    'FontWeight','bold','FontSize',11.5,'BackgroundColor','n');

xa = [.38 .49];
ya = [.81 .75];
line_a = annotation('arrow',xa,ya);
line_a.LineWidth = 0.5;

pos = [-0.05 -4 0.1 5];
rectangle('Position',pos,'Curvature',[1 1],'LineStyle','-','LineWidth',0.8)

% ZOOM
axes('Position',[.16 .7 .21 .21])
box on
plot(Delta_theta,10*log10(G_NF_NF),'r')
hold on
plot(Delta_theta,10*log10(G_NF_NF_18),'k')
plot(Delta_theta,10*log10(G_NF_NF_19),'b')
plot(Delta_theta,10*log10(G_NF_FF_COL),'g')
axis([-0.11 0.11 -4.5 1])
grid on

% 2/N Gamma d
dim = [.18 .8 .06 .07];
str = {'$\frac{2}{N\Gamma\bar{d}}$'};
annotation('textbox',dim,'String',str,'verticalalignment','Bottom','HorizontalAlignment','center','FitBoxToText','on','Color','k','EdgeColor','n','interpreter','latex',...
    'FontWeight','bold','FontSize',14,'BackgroundColor','n');

xa = [.262 .262];
ya = [.855 .855];
line_Radi = annotation('doublearrow',xa,ya);
line_Radi.Head1Width = 20;
line_Radi.Head2Width = 20;
line_Radi.Head1Length = 0;
line_Radi.Head2Length = 0;

xa = [.268 .268];
ya = [.855 .855];
line_Radi = annotation('doublearrow',xa,ya);
line_Radi.Head1Width = 20;
line_Radi.Head2Width = 20;
line_Radi.Head1Length = 0;
line_Radi.Head2Length = 0;

xa = [.24 .262];
ya = [.855 .855];
line_Radi = annotation('arrow',xa,ya);
line_Radi.LineWidth = 0.5;
line_Radi.HeadWidth = 7.5;
line_Radi.HeadLength = 7.5;

xa = [.29 .268];
ya = [.855 .855];
line_Radi = annotation('arrow',xa,ya);
line_Radi.LineWidth = 0.5;
line_Radi.HeadWidth = 7.5;
line_Radi.HeadLength = 7.5;

% 2/N M d
dim = [.295 .74 .06 .07];
str = {'$\frac{2}{N M \bar{d}}$'};
annotation('textbox',dim,'String',str,'verticalalignment','Bottom','HorizontalAlignment','center','FitBoxToText','on','Color','k','EdgeColor','n','interpreter','latex',...
    'FontWeight','bold','FontSize',14,'BackgroundColor','n');

xa = [.252 .252];
ya = [.735 .735];
line_Radi = annotation('doublearrow',xa,ya);
line_Radi.Head1Width = 20;
line_Radi.Head2Width = 20;
line_Radi.Head1Length = 0;
line_Radi.Head2Length = 0;

xa = [.278 .278];
ya = [.735 .735];
line_Radi = annotation('doublearrow',xa,ya);
line_Radi.Head1Width = 20;
line_Radi.Head2Width = 20;
line_Radi.Head1Length = 0;
line_Radi.Head2Length = 0;

xa = [.23 .252];
ya = [.735 .735];
line_Radi = annotation('arrow',xa,ya);
line_Radi.LineWidth = 0.5;
line_Radi.HeadWidth = 7.5;
line_Radi.HeadLength = 7.5;

xa = [.3 .278];
ya = [.735 .735];
line_Radi = annotation('arrow',xa,ya);
line_Radi.LineWidth = 0.5;
line_Radi.HeadWidth = 7.5;
line_Radi.HeadLength = 7.5;









f_Beam_Patterns_NearField_FarField_1 = figure;
set(f_Beam_Patterns_NearField_FarField_1, 'Position',  [600 580, 560, 420])
set(f_Beam_Patterns_NearField_FarField_1, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_Beam_Patterns_NearField_FarField_1, 'defaultLegendInterpreter','latex');
set(f_Beam_Patterns_NearField_FarField_1, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_Beam_Patterns_NearField_FarField_1, 'defaultLineLineWidth',1.5);
set(f_Beam_Patterns_NearField_FarField_1, 'color','w');
plot(Delta_theta,(G_NF_NF),'r')
hold on
plot(Delta_theta,(G_NF_NF_18),'k')
plot(Delta_theta,(G_NF_NF_19),'b')
plot(Delta_theta,(G_NF_FF_COL),'g')
grid on
xlabel('$\Delta_{\theta}$')
ylabel('$G_{NF,NF}(r, \theta; r`, \theta`)$')
legend('Modular, (11)','Modular, (18)','Modular, (19)', 'Collocated, $\Gamma$ = M')









% Polar coordinates
f_Beam_Patterns_NearField_FarField_Spherical = figure;
set(f_Beam_Patterns_NearField_FarField_Spherical, 'Position',  [1220 60, 660, 520])
set(f_Beam_Patterns_NearField_FarField_Spherical, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_Beam_Patterns_NearField_FarField_Spherical, 'defaultLegendInterpreter','latex');
set(f_Beam_Patterns_NearField_FarField_Spherical, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_Beam_Patterns_NearField_FarField_Spherical, 'defaultLineLineWidth',1.5);
set(f_Beam_Patterns_NearField_FarField_Spherical, 'color','w');

P = polarpattern(Delta_theta*90,(G_NF_NF),'AngleDirection','cw','AngleAtTop',-90,'AngleResolution',30);
add(P,Delta_theta*90,(G_NF_NF_18))
add(P,Delta_theta*90,(G_NF_NF_19))
add(P,Delta_theta*90,(G_NF_FF_COL))





















%% Pruebas
f_Beam_Patterns_NearField_FarField_1 = figure;
set(f_Beam_Patterns_NearField_FarField_1, 'Position',  [600 50, 560, 420])
set(f_Beam_Patterns_NearField_FarField_1, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_Beam_Patterns_NearField_FarField_1, 'defaultLegendInterpreter','latex');
set(f_Beam_Patterns_NearField_FarField_1, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_Beam_Patterns_NearField_FarField_1, 'defaultLineLineWidth',1.5);
set(f_Beam_Patterns_NearField_FarField_1, 'color','w');
plot((G_NF_NF),Delta_theta,'r')
hold on
plot((G_NF_NF_18),Delta_theta,'k')
plot((G_NF_NF_19),Delta_theta,'b')
plot((G_NF_FF_COL),Delta_theta,'g')
grid on
xlabel('$\Delta_{\theta}$')
ylabel('$G_{NF,NF}(r, \theta; r`, \theta`)$')
legend('Modular, (11)','Modular, (18)','Modular, (19)', 'Collocated, $\Gamma$ = M')
