function [ PC ] = CTPC( KLrot, KLidx, KR, LR, PC )
%CTPC - Core Transformed Projected Counterpart
% Computes and updates the projected counterpart of the extended Krylov
% iteration.
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
% September 30, 2016
% daan.camps@cs.kuleuven.be

assert(KLidx(end)==1,'Last CT not at L!');

if size(KR,1) == size(KR,2)+1,
    KR(end,:) = [];
end
if size(LR,1) == size(LR,2)+1,
    LR(end-1:end,end) = CreateRotMat(KLrot(:,end)) * LR(end-1:end,end);
    LR(end,:) = [];
    
end
if size(KLrot,2) == size(LR,1),
    KLrot(:,end) = [];
    KLidx(end) = [];
end

i = size(PC,1);
n = size(KR,1);
j = n-i;

% settle ith rotation back
if i > 0,
    PC(end-1:end,end) = CreateRotMat(RotH(KLrot(:,i))) * PC(end-1:end,end);
end

PC12 = [];
if i > 0,
    PC12 = ((-LR(1:i,1:i)*( KR(1:i,1:i)\KR(1:i,i+1:n) )) + LR(1:i,i+1:n)) / KR(i+1:n,i+1:n);
    %for k=i-1:-1:1
    for k=n-1:-1:1
        if(KLidx(k) == 1) && (k<i)% L rotation
            PC12(k:k+1,:) = CreateRotMat(KLrot(:,k)) * PC12(k:k+1,:);
        elseif (KLidx(k) == 1) && (k>i)% K rotation
            PC12(:,k-i:k-i+1) = PC12(:,k-i:k-i+1) * CreateRotMat(RotH(KLrot(:,k)));
        end
    end
end

PC22 = LR(i+1:n,i+1:n) / KR(i+1:n,i+1:n);
%if i == 0, i = 1; end
for k=n-1:-1:i+1
    if(KLidx(k) == 1) % L rotation, apply left
         PC22(k-i:k-i+1,:) = CreateRotMat(KLrot(:,k)) * PC22(k-i:k-i+1,:);
    else % K rotation, apply H right
         PC22(:,k-i:k-i+1) = PC22(:,k-i:k-i+1) * CreateRotMat(RotH(KLrot(:,k)));
    end
end

PC21 = zeros(j,size(PC,2));

PC = [PC, PC12; PC21, PC22];

end

