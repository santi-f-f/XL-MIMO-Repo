function A_K = Ak_ZeroForcing_Calculation(h_BF)
% Calculation of A_k for Zero Forcing

NM = length(h_BF(:,1));
K = length(h_BF(1,:));

A_K = zeros(K,NM,NM);
B_k = zeros(K,NM,K-1);
for k = 1:K
    h_i = [];    
    for i = 1:K
        if i ~= k
            h_i = [h_i h_BF(:,i)];            
        end
    end
    B_k(k,:,:) = h_i;   

    B_k_temp = squeeze(B_k(k,:,:));
    A_K(k,:,:) = B_k_temp*(B_k_temp'*B_k_temp)^(-1)*B_k_temp';    
end
