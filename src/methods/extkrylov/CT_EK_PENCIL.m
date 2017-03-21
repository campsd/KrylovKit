function [KR,LR] = CT_EK_PENCIL(KLrot, KLidx, KR, LR)
% CT_EK_PENCIL -- creates the unreduced upper Hessenberg pencil from the
% factorized representation
    for i=size(KLrot,2):-1:1
       if(KLidx(i) == 0) %K
           KR(i:i+1,:) = CT_TO_MAT(KLrot(:,i)) * KR(i:i+1,:);
       else %L
           LR(i:i+1,:) = CT_TO_MAT(KLrot(:,i)) * LR(i:i+1,:);
       end
    end
end