function [ HR ] = CT_SK_HESS_BLK( Hrot, Hrow, HR )
%[ HR ] = CT_SK_HESS_BLK( Hrot, Hrow, HR )
% -- Full banded Hessenberg matrix from block-Krylov
%
% This computes the full banded Hessenberg matrix from the
% factorized result of CT_SK_BLK.
%
% INPUT
% Hrot  rotations from block Krylov method, ordered per rhomb (2 x bs^2 x m)
% Hrow  row indices from block Krylov, ordered per rhomb (1 x bs^2 x m)
% HR    upper triangular from block Krylov ((m+1)*bs x m*bs)
%
% OUTPUT
% HR    banded upper Hessenberg ((m+1)*bs x bs)
%
% daan.camps@cs.kuleuven.be
% last edit: April 10, 2017
%
% See also: CT_SK_BLK
for i=size(Hrot,3):-1:1
    for j=size(Hrot,2):-1:1
       HR(Hrow(1,j,i):Hrow(1,j,i)+1,:) = CT_TO_MAT(Hrot(:,j,i))* HR(Hrow(1,j,i):Hrow(1,j,i)+1,:);
    end
end
end

