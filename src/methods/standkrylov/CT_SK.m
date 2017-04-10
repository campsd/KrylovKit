function [ V, Hrot, HR ] = CT_SK( mv, V, Hrot, HR, m )
%[ V, Hrot, HR ] = CT_SK( mv, V, Hrot, HR, m ) 
%-- Core transformed Standard Krylov iteration
% 
% INPUT
% mv    handle to function that computes the matrix vector product
% V     Krylov basis to start from (N x k+1)
%	for initial run: V contains start vector of norm 1
% Hrot  Core transformations for upper Hessenberg matrix (2xk)
%	for initial run: Hrot is 2x0 empty array
% HR    upper triangular for upper Hessenberg matrix (k+1 x k)
%	for initial run: HR is 1x0 empty array
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
% last edit: April 10, 2017
%
% See also: CT_SK_HESS, CT_SK_TO_EK_LEFT, CT_SK_TO_EK_RIGHT,
% CT_SK_IR_SS
    tol = 1e-14; %breakdown tolerance
    start_idx = size(V,2);
    % ensure startvec is normalized
    if start_idx == 1
        V(:,1) = V(:,1)/norm(V(:,1));
    end
    
    % zero padding arrays to appropriate size
    V = padarray( V, [0 m], 0, 'post');
    HR = padarray( HR, [m m], 0, 'post');
    Hrot = padarray( Hrot, [0 m], -1, 'post');
    
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
