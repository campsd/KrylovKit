function [ KR, LR ] = CT_EK_PENCIL_BLK( KLrot, KLrow, KLidx, KR, LR )
%CT_EK_PENCIL_BLK Summary of this function goes here
%   Detailed explanation goes here

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
