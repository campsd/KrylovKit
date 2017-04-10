function [ HR ] = CT_SK_HESS( Hrot, HR )
% [ HR ] = CT_SK_HESS( Hrot, HR )
% -- Full upper Hessenberg matrix from factorized format
%
% This computes the dense, unreduced Hessenberg matrix from the
% factorized version from CT_SK
%
% INPUT
% Hrot  rotations for upper Hessenberg (2xm)
% HR    upper triangular for upper Hessenberg (m+1 x m)
%
% OUTPUT
% HR    upper Hessenberg (m+1 x m)
%
% daan.camps@cs.kuleuven.be
% last edit: April 10, 2017
%
% See also: CT_SK
for i=size(Hrot,2):-1:1
    HR(i:i+1,:) = CT_TO_MAT(Hrot(:,i)) * HR(i:i+1,:);
end
end

