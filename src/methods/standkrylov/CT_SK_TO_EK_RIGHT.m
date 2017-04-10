function [ V, KLrot, KLidx, KR, LR, r, em ] = CT_SK_TO_EK_RIGHT( V, Hrot, HR, s )
% [ V, KLrot, KLidx, KR, LR ] = CT_SK_TO_EK_RIGHT( V, Hrot, HR, s )
%   -- converts a standard Krylov recurrence to an extended Krylov
%   recurrence via an initial removal from the right.
%
% This is the method from Mach et al. (2013). It transforms a standard
% Krylov subspace K_m(A,v) to an approximate extended Krylov subspace 
% K_ext(A,v) with the same start vector, but the residual is affected 
% which makes that the result approximate. The error can be checked by
% inspecting the residual.
%
% INPUT
% V		Krylov basis (N x m+1)	
% Hrot		Core transformations upper Hessenberg (m)
% HR		Upper triangular from upper Hessenberg (m+1 x m)
% s		selection vector to which the recurrence should be
%		transformed 0:move to K - 1:stay at L) (m)
% Recurrence that holds: 
%	A * V(:,1:end-1) = V * mat(Hrot) * HR
%
% OUTPUT
% V		adjusted extended Krylov basis (N x m)
% KLrot		adjusted Core transformations L,K pencil (m-1)
% KLidx		indicates side on which core transformation acts 
%		(0 = K // 1 = L) (m-1)
% KR		upper triangular Hessenberg K (m x m)
% LR		upper triangular Hessenberg L (m x m)
% r		residual vector (N x 1)
% em		adjusted mth canonical vector (1 x m)
% Recurrence that holds: 
% 	A * V * mat(KLrot(KLidx==0)) * KR = 
%	V * mat(KLrot(KLidx==1)) * LR + r * em
%
% Before the transformation, em is the mth canonical basis vector.
% During the transformation, core transformations are applied to em which
% pushes the energy to the front. The approximate extended Krylov will be 
% accurate when the recurrence is cut-off at a point where the elements of 
% em are still small. For more details, consult Mach et al. (2013)
%
% daan.camps@cs.kuleuven.be
% last edit: April 10,2017
%
% See also: CT_SK, CT_EK, CT_SK_TO_EK_LEFT, CT_EK_PENCIL
    m = size(Hrot,2);
    CTSV = zeros(3,0);
    CTSW = zeros(3,0);
    KLidx = ones(1,m);
    KR = eye(m,m);
    
    % First we explicitly extract the residual from the recurrence
    HR(m:m+1,m) = CT_TO_MAT(Hrot(:,m)) * HR(m:m+1,m);
    hmp1m = HR(m+1,m); vmp1 = V(:,m+1); em = transpose(flipud(eye(m,1)));
    HR(m+1,:) = []; V(:,m+1) = []; Hrot(:,m) = [];
    
    % Now we modify the square part of H
    for i=1:m-1
        if s(i) == 1 % stay at L
            %KLidx(i) = 1;
        else % move to K
            KLidx(i) = 0;
            for j= m-1:-1:i
                HR(j:j+1,j:m) = CT_TO_MAT(Hrot(:,j)) * HR(j:j+1,j:m);
                [cos,sin,~]=CT_GIV(HR(j+1,j+1),HR(j+1,j));
                Grot = [cos; sin];
                HR(1:j+1,j:j+1) = HR(1:j+1,j:j+1) * CT_TO_MAT(Grot);
                Hrot(:,j) = Grot; % 
            end
            CTSWc = [Hrot(:,m-1:-1:i); m-1:-1:i];
            % we will need to apply these from left to right
            CTSW = [CTSW CTSWc]; 
            CTSV = [CTSV CTSWc(:,1:end-1)];
            Hrot(:,i+1:m-1) = CT_H(Hrot(:,i+1:m-1));
        end
    end
   
    % Apply to V
    Q = eye(m);
    for i = size(CTSV,2):-1:1
        Q(CTSV(3,i):CTSV(3,i)+1,:) = CT_TO_MAT(CTSV(1:2,i)) * Q(CTSV(3,i):CTSV(3,i)+1,:);
    end
    V = V*Q;
    
    % Apply to residual (W)
    for i = 1:size(CTSW,2)
        em(CTSW(3,i):CTSW(3,i)+1) = em(CTSW(3,i):CTSW(3,i)+1) * CT_TO_MAT(CTSW(1:2,i));
    end
    
    LR = HR;
    KLrot = Hrot;
    r = hmp1m * vmp1;
end
