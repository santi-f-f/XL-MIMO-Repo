function C = C_MMSE_Calculation(h_BF, P)
% Calculation of C for MMSE 

NM = length(h_BF(:,1));
K = length(h_BF(1,:));

C = zeros(length(P),K,NM,NM);
for pp = 1:length(P)
    for k = 1:K
        C_temp = 0;        
        for i = 1:K
            if i ~= k
                C_temp = C_temp + P(pp)*h_BF(:,i)*h_BF(:,i)';                
            end
        end
        C(pp,k,:,:) = C_temp + eye(NM);        
    end    
end
