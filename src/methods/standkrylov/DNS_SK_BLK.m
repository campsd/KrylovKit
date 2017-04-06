function [ V, H ] = DNS_SK_BLK( mv, V, H, bs, m )
% [ V, H ] = DNS_SK_BLK( mv, V, H, bs, m )
% -- Dense Standard block-Krylov
%
% INPUT
% mv    function that computes matvec product
% V     existing Krylov basis
% H     existing banded upper Hessenberg
% bs    blocksize
% m     number of blocks to add to subspace
%
% OUTPUT
% V     extended Krylov basis containing bs*m more vectors
% Hrot  updated banded upper Hessenberg
%
% daan.camps@cs.kuleuven.be
!
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

