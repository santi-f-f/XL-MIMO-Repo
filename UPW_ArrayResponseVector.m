function a_ARV = UPW_ArrayResponseVector(r, theta, N, M, Gamma, fc)
% UPW-based far-field array response vector: Observation of the user based in the in the UPW model
%
% r is the distance of the user to the center of the array
% theta is the angle with respect to the positive x-axis
% y is the position of the m-th element within module n
% N is the number of modules
% M is the number of antenna elements within each module
% Gamma is the modular separation parameter (dependent on the discontinuous surface of practical installation structure) Gamma >= M
% fc is the operating frequency

c = physconst('LightSpeed');        % Speed of light
lambda = c/fc;                      % Signal wavelength
d = lambda/2;                       % Inter-element spacing for antennas within each module

NN = -(N-1)/2:(N-1)/2;
MM = -(M-1)/2:(M-1)/2;

p = zeros(1,length(NN));
for n = 1:length(NN)
    p(n) = exp(1i*2*pi*NN(n)*Gamma*d*sin(theta)/lambda);
end

b = zeros(1,length(MM));
for m = 1:length(MM)
    b(m) = exp(1i*2*pi*MM(m)*d*sin(theta)/lambda);
end

a_ARV(:,1) = exp(-1i*2*pi.*r/lambda)*kron(p,b);                            % Eq. 4