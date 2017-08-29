function [ V, Hrot, Hrow, HR ] = CT_SK_BLK( mv, V, Hrot, Hrow, HR, bs, m )
%[ V, Hrot, Hrow, HR ] = CT_SK_BLK( mv, V, Hrot, Hrow, HR, bs, m ) 
% -- Core transformed Standard block-Krylov
%
% This function implements the block Arnoldi iteration with two passes of
% modified Gram-Schmidt in every step. The banded upper Hessenberg matrix
% is conputed and updated in QR factorized format.
%
% The function can be called in two ways:
% OPTION 1: [ V, Hrot, Hrow, HR ] = CT_SK_BLK( mv, V, bs, m ) -> m block
% Arnoldi steps are computed with startblock V of blocksize bs, mv is a
% function handle for the matvec. If the starting block is not orthonormal,
% it is explicitly orthonormalized.
%
% OPTION 2: [ V, Hrot, Hrow, HR ] = CT_SK_BLK( mv, V, Hrot, Hrow, HR, bs, m )
% -> an existing block Arnoldi iteration is expanded with m blocks.
%
% More details on input and output arguments:
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
% last edit: May 22, 2017
%
% See also: CT_TO_DNS_SK_HESS_BLK, CT_SK_TO_EK_RIGHT_BLK, DNS_SK_BLK

% We store the rotations Hrot in a 2 x bs^2 x m array. Each slice stores a
% rhomb such that from front to back the descending pattern forms.
    
    % input parsing
    start_idx = size(V,2);
    if nargin == 4 %we expect mv, Vstart, bs, m
        bs = Hrot; m = Hrow;
        assert(start_idx==bs,'error: incorrect start block');
        start_blck = 1;
        if norm(V'*V - eye(bs),'fro')>1e-13
            warning('starting block must be unitary, will be orthogonalized');
            [V,~] = qr(V,0);
        end
        % check if steps need to be taken
        if m == 0
            return
        end
        % create arrays of appropriate size
        V(end,bs*(m+1)) = 0;
        HR = zeros(bs*(m+1), bs*m);
        Hrot = zeros(2, bs^2, m);
        Hrow = zeros(1, bs^2, m);
    elseif nargin == 7 %default case
        assert(rem(start_idx,bs)==0,'error: incorrect blocksize');
        start_blck = start_idx/bs;
        assert(size(HR,1)==start_idx,'error: dimension mismatch');
        assert(size(HR,2)==start_idx-bs,'error: dimension mismatch');
        assert(size(Hrot,2)==bs^2,'error: dimension mismatch');
        assert(size(Hrot,3)==start_blck-1,'error: dimension mismatch');
        assert(size(Hrow,2)==bs^2,'error: dimension mismatch');
        assert(size(Hrow,3)==start_blck-1,'error: dimension mismatch');
        % check if steps need to be taken
        if m == 0
            return
        end
        % zero padding arrays to appropriate size
        V(end,start_idx + bs*m) = 0;
        HR(start_idx + bs*m, start_idx + bs*(m-1)) = 0;
        Hrot(end,end,start_blck+m-1) = 0;
        Hrow(end,end,start_blck+m-1) = 0;
    else
        error('unsupported input specification');
    end

    tol = 1e-14; %breakdown tolerance  
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
           % handle arrays
           Hrot(:,:,start_blck+i-1:end) = [];
           Hrow(:,:,start_blck+i-1:end) = [];
           HR(:,start_idx + bs*(i-1):end) = [];
           HR(start_idx + bs*i:end,:) = [];
           V(:,start_idx + bs*i:end) = [];
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