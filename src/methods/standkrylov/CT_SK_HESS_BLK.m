function [ HR ] = CT_SK_HESS_BLK( Hrot, Hrow, HR )
%CT_SK_HESS_BLK -- Dense banded Hessenberg matrix from block-Krylov
%algorithm

for i=size(Hrot,3):-1:1
    for j=size(Hrot,2):-1:1
       HR(Hrow(1,j,i):Hrow(1,j,i)+1,:) = CT_TO_MAT(Hrot(:,j,i))* HR(Hrow(1,j,i):Hrow(1,j,i)+1,:);
    end
end

end

