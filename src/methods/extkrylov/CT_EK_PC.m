function [ PC ] = CT_EK_PC( KLrot, KLidx, KR, LR, PC )
% CT_EK_PC - Core Transformed Projected Counterpart
% Computes and updates the projected counterpart of the extended Krylov
% iteration in an iterative fashion.
%
% Assumption: KLidx(end) == 1, otherwise the PC is incorrect
%
% INPUT
% KLrot     Array with Core Transformations
% KLidx     Position of Core Transformations (0=K, 1=L)
% KR        Upper triangular matrix of K
% LR        Upper triangular matirx of L
% PC        Previous PC of size < KR
%
% OUTPUT
% PC        Updated PC
%
% November 16, 2016
% daan.camps@cs.kuleuven.be

assert(KLidx(end)==1,'Last CT not at L!');

% resize arrays
if size(KR,1) == size(KR,2)+1
    KR(end,:) = [];
end
if size(LR,1) == size(LR,2)+1
    LR(end-1:end,end) = CT_TO_MAT(KLrot(:,end)) * LR(end-1:end,end); % last rotation at L
    LR(end,:) = [];
    
end
if size(KLrot,2) == size(LR,1)
    KLrot(:,end) = [];
    KLidx(end) = [];
end

i = size(PC,1); % initial size of projected counterpart
n = size(KR,1); % final size of projected counterpart
j = n-i; % difference in size

if (i==n)
    % nothing to do
    return
end

PC12 = [];
if i > 0
    LKi = LR(1:i,1:i)/ KR(1:i,1:i);
    %PC12 = ((-LR(1:i,1:i)*( KR(1:i,1:i)\KR(1:i,i+1:n) )) + LR(1:i,i+1:n)) / KR(i+1:n,i+1:n); 
    PC12 = ((-LKi * KR(1:i,i+1:n) ) + LR(1:i,i+1:n)) / KR(i+1:n,i+1:n); 
end
PC22 = LR(i+1:n,i+1:n) / KR(i+1:n,i+1:n);
PCX2 = [PC12; PC22];

for k=n-1:-1:1
    if(KLidx(k) == 1) % L rotation, apply left
         PCX2(k:k+1,:) = CT_TO_MAT(KLrot(:,k)) * PCX2(k:k+1,:);
    elseif (KLidx(k) == 0) && (k>i)% K rotation, apply H right
         PCX2(:,k-i:k-i+1) = PCX2(:,k-i:k-i+1) * CT_TO_MAT(RotH(KLrot(:,k)));
    end
end

PC21 = zeros(j,size(PC,2));
if i > 0
    G = CT_TO_MAT(KLrot(:,i)) * [LKi(end,end);0];
    PC21(1,end) = G(2);
    if(KLidx(i-1) ~= 1)
        k = i-1;
        while (KLidx(k) == 0)
            PC21(1,k:k+1) = PC21(1,k:k+1) * CT_TO_MAT(RotH(KLrot(:,k)));
            k = k - 1; 
        end
    end
end

PC = [[PC; PC21] PCX2];

end

