function [ V, H ] = DNS_SK( mv, V, H, m )
%[ V, H ] = DNS_SK( mv, V, H, m ) 
%-- dense (non-factorized) standard Krylov algorithm
%
% This function implements the Arnoldi iteration with two passes of
% modified Gram-Schmidt in every step.
%
% The function can be used in two ways:
% OPTION 1: [ V, H ] = DNS_SK( mv, V, m ) -> m Arnoldi steps are
% computed with startvector V and function handle mv
%
% OPTION 2: [ V, H ] = DNS_SK( mv, V, H, m )  -> an existing
% Arnoldi iteration is expanded with m vectors.
%
% More details on input and output arguments:
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
% last edit: May 22, 2017
%
% See also: DNS_SK_BLK, CT_SK

    % input parsing
    start_idx = size(V,2);
    if nargin == 3 %we expect mv, Vstart, m
       m = H; % third argument is actually the number of steps
       assert(start_idx==1,'error: should be a single start vector');
       V(:,1) = V(:,1)/norm(V(:,1)); % ensure startvec is normalized
       V(end,m+1) = 0;
       H = zeros(m+1,m);
    elseif nargin == 4 %default case
       assert(size(H,1)==start_idx,'error: dimension mismatch');
       assert(size(H,2)==start_idx-1,'error: dimension mismatch');
       % zero padding arrays to appropriate size
       V(end,start_idx+m) = 0;
       H(start_idx+m,start_idx+m-1) = 0;
    else
       error('unsupported input specification');
    end
    
    % check if any steps need to be taken
    if m == 0
        return
    end

    tol = 1e-14; %breakdown tolerance
    % main loop
    for i=1:m
        % matvec
        w = mv(V(:,start_idx + i - 1)); 
        % MGS
        h = zeros(start_idx + i,1);
        for kk=1:2 % two passes to ensure orthogonal vectors
            for j=1:start_idx + i - 1
                hc = V(:,j)'*w;
                h(j) = h(j) + hc;
                w = w - V(:,j)*hc;
            end
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
        % Store new base vec and update Hessenberg
        V(:,start_idx + i) = w / h(start_idx + i);
        H(1:start_idx+i,start_idx+i-1) = h(1:start_idx+i);
    end
end