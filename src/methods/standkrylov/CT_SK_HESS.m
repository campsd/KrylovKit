function [ HR ] = CT_SK_HESS( Hrot, HR )
%CT_SK_HESS -- Dense Hessenberg matrix from CT_SK
%
% This computes the dense, unreduced Hessenberg matrix from the CT
% factorized version
%
% IN
% Hrot  rotations
% HR    upper triangular
%
% OUT
% HR    upper Hessenberg

for i=size(Hrot,2):-1:1
    HR(i:i+1,:) = CT_TO_MAT(Hrot(:,i)) * HR(i:i+1,:);
end

end

