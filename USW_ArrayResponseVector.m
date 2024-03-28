function a_ARV = USW_ArrayResponseVector(r, theta, y_n, N, M, fc)
% USW-based near-field array response vector: Observation of the user based in the Sub-array based USW model for distinct AoAs/AoDs (Eqs. 5, 6, 7)
%
% r is the distance of the user to the center of the array
% theta is the angle with respect to the positive x-axis
% y_n is the position of the m-th element within module n (sub-array based USW model for distinct AoAs/AoDs)
% N is the number of modules
% M is the number of antenna elements within each module
% fc is the operating frequency

c = physconst('LightSpeed');        % Speed of light
lambda = c/fc;                      % Signal wavelength
d = lambda/2;                       % Inter-element spacing for antennas within each module

NM = N*M;
MM = -(M-1)/2:(M-1)/2;

r_n0 = zeros(1, length(y_n));
sin_theta_n = zeros(1, length(y_n));
a_ARV = zeros(NM,1);

b = zeros(length(MM),length(y_n));
a_ARV_temp = [];
for in = 1:length(y_n)
    r_n0(1,in) = sqrt(r^2 - 2*r*y_n(in)*sin(theta) + y_n(in)^2);                                    % Eq. 5

    sin_theta_n(1,in) = (r*sin(theta) - y_n(in))/r_n0(1,in);                                        % Eq. 6

    for m = 1:length(MM)
        b(m,in) = exp(1i*2*pi*MM(m)*d*sin_theta_n(1,in)/lambda);
    end

    a_ARV_temp = [a_ARV_temp; (exp(-1i*2*pi*r_n0(1,in)/lambda)*b(:,in))];               % eq. 7
end

a_ARV(:,1) = a_ARV_temp;

