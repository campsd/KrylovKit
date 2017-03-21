function [ V, Hrot, HR ] = CT_SK( mv, V, Hrot, HR, m )
%CT_SK -- Core transformed Standard Krylov
% 
% INPUT
% mv    handle to function that computes the matrix vector product
% V     Krylov basis to start from
% Hrot  rotations from upper Hessenberg matrix
% HR    upper triangular from upper Hessenberg matrix
% m     number of steps to compute
%
% OUTPUT
% V     enlarged Krylov basis with m additional vectors
% Hrot  rotations from upper Hessenberg matrix
% HR    upper triangular
%
% Recursion: A*V(:,1:end-1) = V*H
%
% daan.camps@cs.kuleuven.be
% March 21, 2017
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
            h(kk:kk+1) = CT_TO_MAT(RotH(Hrot(:,kk))) * h(kk:kk+1);
        end
        % Compute new rotation
        [c,s,r] = CT_GIV(h(start_idx+i-1),h(start_idx+i));
        h(start_idx+i-1) = r;
        HR(1:start_idx+i-1,start_idx+i-1) = h(1:start_idx+i-1);
        Hrot(:,start_idx+i-1) = [conj(c); -s];
    end
end