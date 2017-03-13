function [KR,LR] = CONS_CTEK_PENCIL(KLrot, KLidx, KR, LR)
    for i=length(KLrot):-1:1
       if(KLidx(i) == 0) %K
           KR(i:i+1,:) = CreateRotMat(KLrot(:,i)) * KR(i:i+1,:);
       else %L
           LR(i:i+1,:) = CreateRotMat(KLrot(:,i)) * LR(i:i+1,:);
       end
    end
end