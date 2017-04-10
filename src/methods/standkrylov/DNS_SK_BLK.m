function [ V, H ] = DNS_SK_BLK( mv, V, H, bs, m )
% [ V, H ] = DNS_SK_BLK( mv, V, H, bs, m )
% -- Dense Standard block-Krylov
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
% last edit: April 10, 2017
%
% See also: CT_SK_BLK, DNS_SK

% TODO breakdown check
    tol = 1e-14; %breakdown tolerance
    start_idx = size(V,2);
    assert(mod(start_idx,bs)==0,'incorrect blocksize');
    start_blck = start_idx/bs;
    
    if start_idx == bs
        assert(norm(V'*V - eye(bs),'fro')<1e-14,'starting block must be unitary');
    end
    
    % zero padding arrays to appropriate size
    V = padarray( V, [0 bs*m], 0, 'post');
    H = padarray( H, [bs*m bs*m], 0, 'post');

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
    end
end

