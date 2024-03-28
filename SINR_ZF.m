function SINR = SINR_ZF(P, A_K, h, h_BF)
% Eq. 43 from "https://arxiv.org/abs/2308.11289"

K = length(h_BF(1,:));              % Number of users
NM = length(h_BF(:,1));             % Total number of array elements

% Preallocation
SINR = zeros(1,K);
num_SINR = zeros(1,K);
den_SINR = zeros(1,K);
for k = 1:K
    temp_A_K = squeeze(A_K(k,:,:));
    num_SINR(k) = P*abs((h_BF(:,k)')*(eye(NM) - temp_A_K)'*h(:,k))^2/(norm((eye(NM) - temp_A_K)*h_BF(:,k))^2);

    den_SINR_temp = 0;
    for i = 1:K
        if i~=k
            den_SINR_temp = den_SINR_temp + P*abs((h_BF(:,k)')*(eye(NM) - temp_A_K)'*h(:,i))^2/(norm((eye(NM) - temp_A_K)*h_BF(:,k))^2);
        end
    end
    den_SINR(k) = den_SINR_temp;

    SINR(1,k) = num_SINR(k)/(den_SINR(k) + 1);
end