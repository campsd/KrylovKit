function [V, K, L] = EXQZIR(V,K,L,s)
    % filters the RKS triple (V,K,L) with the shifts in vector s
    for i=1:length(s)
       [Q,~] = qr(L - s(i)*K);
       g = Q(:,end)'*(conj(s(i))*L + K);
       Z = Z2(g);
       Q = Q(:,1:end-1);
       K = Q'*K*Z;
       L = Q'*L*Z;
       V = V*Q;
    end
end

function Z = Z1(g)
    % first option for Z (Karl)
    [Z,~] = qr(g');
    Z = Z(:,2:end);
end

function Z = Z2(g)
    % second option for Z (Paper De Samblanx)
    k = size(g,2);
    Z = zeros(k,k-1);
    %Z(1,1) = -g(2); Z(2,1) = g(1);
    for j=1:k-1
        [~,idx] = max(abs(g(1:j)));
        Z(idx,j) = -g(j+1); Z(j+1,j) = g(idx);
    end
    [Z,~] = qr(Z);
    Z = Z(:,1:end-1);
end