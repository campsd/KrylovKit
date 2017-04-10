function [ V, Hrot, Hrow, HR ] = CT_SK_BLK( mv, V, Hrot, Hrow, HR, bs, m )
%[ V, Hrot, HR ] = CT_SK_BLK( mv, V, Hrot, HR, bs, m ) 
% -- Core transformed Standard block-Krylov
%
% INPUT
% mv    function that computes matvec product
% V     existing Krylov basis (N x (k+1)*bs)
%	for initial run: V should contain bs orthonormal startvectors
% Hrot  existing core transformations, stored per rhomb (2xbs^2xk)
%	for initial run: Hrot is 2xbs^2x0 empty array
% Hrow  keeps track of the row indices of the core transformations (1xbs^2xk)
%	for initial run: Hrow is 1xbs^2x0 empty array
% HR    existing upper triangular for banded upper Hessenberg ((k+1*bs x k*bs))
%	for initial run: HR is bsx0 empty array
% bs    blocksize
% m     number of blocks to add to subspace
%
% OUTPUT
% V     enlarged Krylov basis (N x (k+m+1)*bs)
% Hrot  updated rotations stored per rhomb (2 x bs^2 x k+m)
% Hrow  updated row indices (1 x bs^2 x k+m)
% HR    updated upper triangular ((k+m+1)*bs x (k+m)*bs)
%
% The recurrence that holds throughout the algorithm is
%	A * V(:,1:end-bs) = V * H
%
% daan.camps@cs.kuleuven.be
% last edit: April 10, 2017
%
% See also: CT_SK_HESS_BLK, CT_SK_TO_EK_RIGHT_BLK

% We store the rotations Hrot in a 2 x bs^2 x m array. Each slice stores a
% rhomb such that from front to back the descending pattern forms.
    tol = 1e-14; %breakdown tolerance
    start_idx = size(V,2);
    assert(mod(start_idx,bs)==0,'incorrect blocksize');
    start_blck = start_idx/bs;
    
    if start_idx == bs
        assert(norm(V'*V - eye(bs),'fro')<1e-14,'starting block must be unitary');
    end
    
    % zero padding arrays to appropriate size
    V = padarray( V, [0 bs*m], 0, 'post');
    HR = padarray( HR, [bs*m bs*m], 0, 'post');
    Hrot = padarray( Hrot, [0 0 m], 0, 'post');
    Hrow = padarray( Hrow, [0 0 m], 0, 'post');
    
    % main loop
    for i=1:m
       %matvecs
       W = mv(V(:, start_idx+1+(i-2)*bs:start_idx+(i-1)*bs));
       %MGS
       Ht = zeros(start_idx+i*bs,bs);
       for j=1:start_blck + i-1
           for kk=1:2 %two passes to ensure orthogonal vectors
               Hc = V(:,(j-1)*bs+1:j*bs)'*W; %(bs*m)xbs
               Ht((j-1)*bs+1:j*bs,:) = Ht((j-1)*bs+1:j*bs,:) + Hc;
               W = W - V(:,(j-1)*bs+1:j*bs)*Hc;
           end
       end
       % Update recurrence
       [V(:,start_idx+1+(i-1)*bs:start_idx+i*bs), Ht(end-bs+1:end,:)] = qr(W,0);
       
       % Check for breakdown
       if min(svd(Ht(end-bs+1:end,:))) < tol
           warning('breakdown of block-Krylov process, halted early');
           % handle arrays TODO
           break
       end
       
       % Apply previous rotations to Ht
       for kk=1:start_blck + i-2 % loop over blocks of rotations
        for jj=1:size(Hrot,2) 
           Ht(Hrow(1,jj,kk):Hrow(1,jj,kk)+1,:) = CT_TO_MAT(CT_H(Hrot(:,jj,kk))) * Ht(Hrow(1,jj,kk):Hrow(1,jj,kk)+1,:);
        end
       end
       
       % Compute new rotations
       RhombRot = zeros(2,bs^2);
       RhombRow = zeros(1,bs^2);
       start_row = (i + start_blck -1)*bs + 1;
       for kk=1:bs
           for jj=1:bs
               [c,s,r] = CT_GIV(Ht(start_row-kk+jj-1,jj),Ht(start_row-kk+jj,jj));
               Ht(start_row-kk+jj-1,jj) = r; Ht(start_row-kk+jj,jj) = 0;
               Ht(start_row-kk+jj-1:start_row-kk+jj,jj+1:end) = CT_TO_MAT([c;s]) * Ht(start_row-kk+jj-1:start_row-kk+jj,jj+1:end);
               RhombRot(:,(kk-1)*bs + jj) = [conj(c); -s];
               RhombRow((kk-1)*bs + jj) = start_row-kk+jj-1;
           end
       end
       HR(1:start_idx+i*bs,start_idx+1+(i-2)*bs:start_idx+(i-1)*bs) = Ht;
       
       Hrot(:,:,start_blck + i-1) = RhombRot;
       Hrow(:,:,start_blck + i-1) = RhombRow;
    end
end


    
