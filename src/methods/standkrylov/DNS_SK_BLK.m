function [ V, H ] = DNS_SK_BLK( mv, V, H, bs, m )
% [ V, H ] = DNS_SK_BLK( mv, V, H, bs, m )
% -- Dense Standard block-Krylov
%
% This function implements the block Arnoldi iteration with two passes of
% modified Gram-Schmidt in every iteration.
%
% The function can be called in two ways:
% OPTION 1: [ V, H ] = DNS_SK_BLK( mv, V, bs, m ) -> m block Arnoldi steps
% are computed with startblock V of blocksize bs, mv is a function handle
% for the matvec. If the startingblock is not orthonormal, it is explicitly
% orthonormalized
%
% OPTION 2: [ V, H ] = DNS_SK_BLK( mv, V, H, bs, m ) -> an existing block
% Arnoldi iteration is expanded with m blocks
%
% More details on input and output arguments:
%
% INPUT
% mv    function that computes matvec product
% V     existing Krylov basis (N x (k+1*bs))
%	for initial run: V should contain bs orthonormal startvectors
% H     existing banded upper Hessenberg ((k+1)*bs x k*bs)
%	for initial run: H is 1x0 empty array
% bs    blocksize
% m     number of blocks to add to subspace
%
% OUTPUT
% V     enlarged Krylov basis (N x (k+m+1)*bs)
% H  	updated banded upper Hessenberg ((k+m+1)*bs x (k+m)*bs)
%
% daan.camps@cs.kuleuven.be
% last edit: May 23, 2017
%
% See also: CT_SK_BLK, DNS_SK, DNS_TO_CT_SK_HESS_BLK

    % input parsing
    start_idx = size(V,2);
    if nargin == 4 %we expect mv, Vstart, bs, m
         m = bs; bs = H;
         assert(start_idx==bs,'error: incorrect start block');
         start_blck = 1;
         if norm(V'*V - eye(bs),'fro') > 10 * start_idx * eps
             warning('starting block must be unitary, will be orthogonalized');
            [V,~] = qr(V,0);
         end
         % check if steps need to be taken
         if m == 0
            return
         end
         % create arrays of appropriate size
         V(end,bs*(m+1)) = 0;
         H = zeros(bs*(m+1), bs*m);
    elseif nargin == 5 %we expect mv, V, H, bs, m
        assert(rem(start_idx,bs)==0,'error: incorrect blocksize');
        start_blck = start_idx/bs;
        assert(size(H,1)==start_idx,'error: dimension mismatch');
        assert(size(H,2)==start_idx-bs,'error: dimension mismatch');
        % check if steps need to be taken
        if m == 0
            return
        end
        % zero padding arrays to appropriate size
        V(end,start_idx + bs*m) = 0;
        H(start_idx + bs*m, start_idx + bs*(m-1)) = 0;
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
       H(1:start_idx+i*bs,start_idx+1+(i-2)*bs:start_idx+(i-1)*bs) = Ht;
       
       % Check for breakdown
       if min(svd(Ht(end-bs+1:end,:))) < tol
           warning('breakdown of block-Krylov process, halted early');
           % handle arrays
           H(:,start_idx + bs*(i-1):end) = [];
           H(start_idx + bs*i:end,:) = [];
           V(:,start_idx + bs*i:end) = [];
           break
       end
       
    end
end