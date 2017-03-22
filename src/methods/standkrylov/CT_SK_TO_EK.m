function [ V, KLrot, KLidx, KR, LR ] = CT_SK_TO_EK( V, Hrot, HR, s )
%CT_SK_TO_EK  -- converts a standard Krylov recurrence to an extended
%Krylov recurrence
%
% The method is similar to the results from Mach et al. (2013) but doesn't
% influence the residual.
%
% IN
% V     Krylov basis (N x m+1)
% Hrot  Core transformations upper Hessenberg (m)
% HR    Upper triangular (m+1 x m)
% s     selection vector to which the recurrence should be transformed
%       (0:to K - 1:stay at L) (m)
% Recurrence that holds: A * V(:,1:end-1) = V * mat(Hrot) * HR
%
% OUT
% V     adjusted Krylov basis
% KLrot adjusted Core transformations pencil
% KLidx indicates side of core transformations
% KR    upper triangular K
% LR    upper triangular L
% Recurrence that holds: A * V * mat(KLrot(KLidx==0)) * KR = V * mat(KLrot(KLidx==1)) * LR 
%
% March 22, 2017
% daan.camps@cs.kuleuven.be

    m = length(Hrot);
    KR = eye(m+1,m);
    CTSV = zeros(3,0);
    KLidx = zeros(m,1);
    KR = eye(m+1,m);
    
    for i=m:-1:1
       if s(i) == 1 % stay at L
           KLidx(i) = 1;
       else % move to K
           CTSV = [CTSV [Hrot(:,1:i); 1:i]];
           Hrot(:,i) = CT_H(Hrot(:,i));
           for j=1:i-1
               HR(1:j+1,j:j+1) = HR(1:j+1,j:j+1) * CT_TO_MAT(Hrot(:,j));
               [cos,sin,~]=CT_GIV(HR(j,j),HR(j+1,j));
               Grot = [cos,sin];
	           HR(j:j+1,j:end) = CT_TO_MAT(Grot)*HR(j:j+1,j:end);
               Hrot(:,j) = CT_H(Grot);
           end
       end
    end
    LR = HR;
    KLrot = Hrot;
    
    % Apply to basis
    Q = eye(m+1);
    for i=size(CTSV,2):-1:1
        Q(CTSV(3,i):CTSV(3,i)+1,:) = CT_TO_MAT(CTSV(1:2,i)) * Q(CTSV(3,i):CTSV(3,i)+1,:);
    end
    V = V*Q;
end

