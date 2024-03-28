function SINR = SINR_MMSE(P, C, h, h_BF)

K = length(h_BF(1,:));          % Number of users

% Preallocation
SINR = zeros(1,K);
num_SINR = zeros(1,K);
den_SINR = zeros(1,K);

for k = 1:K
    temp_C = squeeze(C(1,k,:,:));
    num_SINR(k) = P*(abs(h_BF(:,k)'*(temp_C^-1)'*h(:,k))^2)/(norm(temp_C^-1*h_BF(:,k))^2);

    den_SINR_temp = 0;
    for i = 1:K
        if i~=k
            den_SINR_temp = den_SINR_temp + P*(abs(h_BF(:,k)'*(temp_C^-1)'*h(:,i))^2)/(norm(temp_C^-1*h_BF(:,k))^2);
        end
    end
    den_SINR(k) = den_SINR_temp;

    SINR(1,k) = num_SINR(k)/(den_SINR(k) + 1);
end