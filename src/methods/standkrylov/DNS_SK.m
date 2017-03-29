function [ V, H ] = DNS_SK( mv, V, H, m )
%DNS_SK -- dense standard Krylov algorithm
%   Arnoldi algorithm
% March 28, 2017

    tol = 1e-14;
    start_idx = size(V,2);
    % ensure startvec is normalized
    if start_idx == 1
        V(:,1) = V(:,1)/norm(V(:,1));
    end
    % zero padding arrays to appropriate size
    V = padarray( V, [0 m], 0, 'post');
    H = padarray( H, [m m], 0, 'post');

    % main loop
    for i=1:m
        % matvec
        w = mv(V(:,start_idx + i - 1)); 
        % MGS
        h = zeros(start_idx + i,1);
        for kk=1:2 % two passes to ensure orthogonal vectors
            hc = V(:,1:start_idx + i - 1)'*w;
            h(1:start_idx + i - 1) = h(1:start_idx + i - 1) + hc;
            w = w - V(:,1:start_idx + i - 1)*hc;
        end

        h(start_idx + i) = norm(w,2);
        % check for breakdown
        if h(start_idx + i) < tol
            warning('breakdown of Krylov process, halted early');
            H(:,start_idx + i-1:end) = [];
            H(start_idx + i:end,:) = [];
            V(:,start_idx + i:end) = [];
            break
        end
        % Store new base vec and update hessenberg
        V(:,start_idx + i) = w / h(start_idx + i);
        H(1:start_idx+i,start_idx+i-1) = h(1:start_idx+i);
    end
end

