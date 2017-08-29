function [ V, Hrot, HR ] = CT_SK( mv, V, Hrot, HR, m )
%[ V, Hrot, HR ] = CT_SK( mv, V, Hrot, HR, m ) 
%-- Core transformed Standard Krylov iteration
% 
% This function implements the Arnoldi iteration with two passes of
% modified Gram-Schmidt in every step. The upper Hessenberg matrix H is
% computed in a QR-factorised format which is updated every step. This
% computation has a cost of O(k) per step.
%
% The function can be used in two ways:
% OPTION 1: [ V, Hrot, HR ] = CT_SK( mv, V, m ) -> m Arnoldi steps are
% computed with startvector V and function handle mv
%
% OPTION 2: [ V, Hrot, HR ] = CT_SK( mv, V, Hrot, HR, m )  -> an existing
% Arnoldi iteration is expanded with m vectors.
%
% More details on input and output arguments:
%
% INPUT
% mv    handle to function that computes the matrix vector product
% V     Krylov basis to start from (N x k+1)
%	for initial run: V contains start vector of norm 1
% Hrot  Core transformations for upper Hessenberg matrix (2xk)
%	for initial run: Hrot is 2x0 empty array or skipped
% HR    upper triangular for upper Hessenberg matrix (k+1 x k)
%	for initial run: HR is 1x0 empty array or skipped
% m     number of steps to take
%
% OUTPUT
% V     enlarged Krylov basis with m additional vectors (N x k+m+1)
% Hrot  rotations for upper Hessenberg matrix (2 x k+m)
% HR    upper triangular for upper Hessenberg matrix (k+m+1 x k+m)
%
% Throughout the algorithm the recursion
%	 A*V(:,1:end-1) = V*H
% holds, with H = mat(Hrot) * HR.
%
% daan.camps@cs.kuleuven.be
% last edit: May 22, 2017
%
% See also: CT_TO_DNS_SK_HESS, CT_SK_TO_EK_LEFT, CT_SK_TO_EK_RIGHT,
% CT_SK_IR_SS, CT_SK_IR_DS, DNS_SK
    
    % input parsing
    start_idx = size(V,2);
    if nargin == 3 %we expect mv, Vstart, m
       assert(start_idx==1,'error: should be a single start vector');
       V(:,1) = V(:,1)/norm(V(:,1)); % ensure startvec is normalized
       m = Hrot; % third argument is actually the number of steps
       % check if steps need to be taken
       if m == 0
        return
       end
       % zero padding arrays
       V(end,m+1) = 0;
       Hrot = zeros(2,m);
       HR = zeros(m+1,m);
    elseif nargin == 5 %default case
       % check if steps need to be taken
       if m == 0
            return
       end
       assert(size(HR,1)==start_idx,'error: dimension mismatch');
       assert(size(HR,2)==start_idx-1,'error: dimension mismatch');
       assert(size(HR,2)==size(Hrot,2),'error: dimension mismatch');
       % zero padding arrays to appropriate size
       V(end,start_idx+m) = 0;
       HR(start_idx+m,start_idx+m-1) = 0;
       Hrot(end,start_idx+m-1) = 0;    
    else
       error('unsupported input specification');
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
        % normalisation
        h(start_idx + i) = norm(w,2);
        % check for breakdown
        if h(start_idx + i) < tol
            warning('breakdown of Krylov process, halted early');
            Hrot(:,start_idx + i -1:end) = [];
            HR(:,start_idx + i-1:end) = [];
            HR(start_idx + i:end,:) = [];
            V(:,start_idx + i:end) = [];
            break
        end
        % Store new base vec
        V(:,start_idx + i) = w / h(start_idx + i);
        % Apply previous rotations to h
        for kk=1:start_idx+i-2
            h(kk:kk+1) = CT_TO_MAT(CT_H(Hrot(:,kk))) * h(kk:kk+1);
        end
        % Compute new rotation
        [c,s,r] = CT_GIV(h(start_idx+i-1),h(start_idx+i));
        h(start_idx+i-1) = r;
        HR(1:start_idx+i-1,start_idx+i-1) = h(1:start_idx+i-1);
        Hrot(:,start_idx+i-1) = [conj(c); -s];
    end
end
