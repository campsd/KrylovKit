function [ V, H ] = DNS_SK( mv, V, H, m )
%[ V, H ] = DNS_SK( mv, V, H, m ) 
%-- dense (non-factorized) standard Krylov algorithm
%
% INPUT
% mv	handle to function that computes the matrix vector product
% V	Krylov basis to start from  (N x k+1)
%	for initial run: V contains start vector of norm 1
% H	upper Hessenberg matrix (k+1 x k)
%	for initial run: HR is 1x0 empty array
% m	number of steps to take
%
% OUTPUT
% V	enlarged Krylov basis with m additional vectors (N x k+m+1)
% H	enlarged upper Hessenberg (k+m+1 x k+m)
%
% Throughout the algorithm the recursion
%	 A*V(:,1:end-1) = V*H
% holds
%
% daan.camps@cs.kuleuven.be
% last edit: April 10, 2017
%
% See also: DNS_SK_BLK, CT_SK
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

