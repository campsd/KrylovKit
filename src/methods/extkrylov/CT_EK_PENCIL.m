function [KR,LR] = CT_EK_PENCIL(KLrot, KLidx, KR, LR)
% [K,L] = CT_EK_PENCIL(KLrot, KLidx, KR, LR) 
% -- creates the upper Hessenberg pencil L,K from the factorized 
% representation
%
% INPUT
% KLrot	Core transformations for L,K pencil (2xm)
% KLidx	Indicates side of core transformation
%	(0 = K // 1 = L) (m)
% KR	Upper triangular K (m+1 x m)
% LR	Upper triangular L (m+1 x m)
%
% OUTPUT
% K	upper Hessenberg (m+1 x m)
% L	upper Hessenberg (m+1 x m)
%
% daan.camps@cs.kuleuven.be
% last edit: April 10, 2017
%
% See also: CT_EK, CT_EK_PC
    for i=size(KLrot,2):-1:1
       if(KLidx(i) == 0) %K
           KR(i:i+1,:) = CT_TO_MAT(KLrot(:,i)) * KR(i:i+1,:);
       else %L
           LR(i:i+1,:) = CT_TO_MAT(KLrot(:,i)) * LR(i:i+1,:);
       end
    end
end
