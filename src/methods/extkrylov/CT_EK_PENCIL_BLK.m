function [ KR, LR ] = CT_EK_PENCIL_BLK( KLrot, KLrow, KLidx, KR, LR )
%[ K, L ] = CT_EK_PENCIL_BLK( KLrot, KLrow, KLidx, KR, LR )
% -- creates the banded upper Hessenberg pencil L,K from the factorized 
% representation. Valid for block-Krylov
%
% INPUT
% KLrot	core transformations banded Hessenberg pencil, stored per rhomb (2x bs^2 x m-1)
% KLrow	row indices core transformations (1 x bs^2 x m-1)
% KLidx	indicates side on which rhomb of core transformations acts
%	(0 = K // 1 = L) (m-1)
% KR	upper triangular Hessenberg K (m*bs x m*bs)
% LR	upper triangular Hessenberg L (m*bs x m*bs)
%
% OUTPUT
% K	banded Hessenberg K (m*bs x m*bs)
% L	banded Hessenberg L (m*bs x m*bs)
%
% daan.camps@cs.kuleuven.be
% last edit: April 10, 2017
%
% See also: CT_SK_BLK, CT_SK_HESS_BLK
for i=size(KLrot,3):-1:1
    for j=size(KLrot,2):-1:1
        if KLrow(1,j,i) ~= -1
           if KLidx(i) == 0 %K
              KR(KLrow(1,j,i):KLrow(1,j,i)+1,:) = CT_TO_MAT(KLrot(:,j,i)) * KR(KLrow(1,j,i):KLrow(1,j,i)+1,:); 
           else % L
              LR(KLrow(1,j,i):KLrow(1,j,i)+1,:) = CT_TO_MAT(KLrot(:,j,i)) * LR(KLrow(1,j,i):KLrow(1,j,i)+1,:);  
           end
        end
    end
end
end

