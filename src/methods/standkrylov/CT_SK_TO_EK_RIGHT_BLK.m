function [ V, KLrot, KLrow, KLidx, KR, LR, F, E ] = CT_SK_TO_EK_RIGHT_BLK( V, Hrot, Hrow, HR, bs, s )
%[ V, KLrot, KLrow, KR, LR, F, E ] = CT_SK_TO_EK_RIGHT_BLK( V, Hrot, Hrow, HR, bs, s )
% -- converts a standard block-Krylov recurrence to an extended block- Krylov 
% recurrence via an initial removal from the right.
%
% This algorithm is based on the work form Mach et al. (2014) on
% approximate rational block-Krylov methods.
%
% INPUT
% V	standard block-Krylov basis (N x (m+1)*bs)
% Hrot	core transformations banded Hessenberg, stored per rhomb (2 x bs^2 x m)
% Hrow  row indices core transformations (1 x bs^2 x m)
% HR 	upper triangular banded upper Hessenberg
% bs	blocksize of block-Krylov subspace
% s	selection vector to which the subspace needs to be transformed
% Recurrence that holds:
%	A * V(:,1:end-bs) = V * mat(Hrot) * HR
%
% OUTPUT
% V	approximate extended block-Krylov basis (N x m*bs)
% KLrot	core transformations banded Hessenberg pencil, stored per rhomb (2x bs^2 x m-1)
% KLrow row indices core transformations (1 x bs^2 x m-1)
% KLidx indicates side on which rhomb of core transformations acts
%	(0 = K // 1 = L) (m-1)
% KR	upper triangular Hessenberg K (m*bs x m*bs)
% LR	upper triangular Hessenberg L (m*bs x m*bs)
% F	residual vectors (N x bs)
% E	adjusted mth canonical vectors (bs x m)
% Recurrence that holds:
%	A * V * mat(KLrot(KLidx==0)) * KR = 
%	V * mat(KLrot(KLidx==1)) * LR + F * E
%
% Before the transformation, E contains canonical basis vectors.
% During the transformation, core transformations are applied to E, which
% push the energy in the vectors to the front. The approximate extended
% block-Krylov subspace will be accurate if it is extracted at a point where
% all elements of E are still small.
%
% daan.camps@cs.kuleuven.be
% last edit: April 10, 2017
%
% See also: CT_SK_BLK, CT_EK_PENCIL_BLK, CT_TO_DNS_SK_HESS_BLK
m = size(Hrot,3); %nb of blocks
CTSV = zeros(3,0);
CTSW = zeros(3,0);
KLidx = ones(1,m);
KR = eye(m*bs,m*bs);
    
% First, we extract the residual from the recurrence
% (a) apply last rhomb to HR
for j=size(Hrot,2):-1:1
    HR(Hrow(1,j,m):Hrow(1,j,m)+1,:) = CT_TO_MAT(Hrot(:,j,m))* HR(Hrow(1,j,m):Hrow(1,j,m)+1,:);
end
% (b) residual
F = V(:,bs*m+1:bs*(m+1)) * HR(bs*m+1:bs*(m+1),bs*(m-1)+1:bs*m);
E = eye(bs*m,bs*m); E = rot90(E(1:bs,:),2);
HR(bs*m+1:bs*(m+1),:) = []; V(:,bs*m+1:bs*(m+1)) = [];
% (c) last half rhomb
Hrot(:,:,m) = -1; Hrow(:,:,m) = -1;
cnt = 1;
for j=1:bs-1
   for k=bs-j:-1:1
       [cos,sin,r] = CT_GIV(HR(j+k+(m-1)*bs-1,j+(m-1)*bs),HR(j+k+(m-1)*bs,j+(m-1)*bs));
       HR(j+k+(m-1)*bs-1,j+(m-1)*bs) = r; HR(j+k+(m-1)*bs,j+(m-1)*bs) = 0;
       HR(j+k+(m-1)*bs-1:j+k+(m-1)*bs,j+(m-1)*bs+1:m*bs) = CT_TO_MAT([cos;sin]) * HR(j+k+(m-1)*bs-1:j+k+(m-1)*bs,j+(m-1)*bs+1:m*bs);
       Hrot(:,cnt,m) = [conj(cos); -sin];
       Hrow(1,cnt,m) = j+k+(m-1)*bs-1;
       cnt = cnt + 1;
   end
end

% Now modify the top square part of H
for i=1:m-1
    if s(i) ==1 %stay at L
        %nothing needs to be done
    else %move the ith rhomb to K
        KLidx(i) = 0;
        % move the half rhomb to the other side
        for j=bs*(bs-1)/2:-1:1
            HR(Hrow(1,j,m):Hrow(1,j,m)+1,:) = CT_TO_MAT(Hrot(:,j,m))* HR(Hrow(1,j,m):Hrow(1,j,m)+1,:);
            [cos,sin,~]=CT_GIV(HR(Hrow(1,j,m)+1,Hrow(1,j,m)+1),HR(Hrow(1,j,m)+1,Hrow(1,j,m)));
            Grot = [cos; sin];
            HR(1:Hrow(1,j,m)+1,Hrow(1,j,m):Hrow(1,j,m)+1) = HR(1:Hrow(1,j,m)+1,Hrow(1,j,m):Hrow(1,j,m)+1) * CT_TO_MAT(Grot);
            Hrot(:,j,m) = Grot;
        end
        for j=m-1:-1:i % move all other rhombs to the other side
            for k=size(Hrot,2):-1:1
                HR(Hrow(1,k,j):Hrow(1,k,j)+1,:) = CT_TO_MAT(Hrot(:,k,j))* HR(Hrow(1,k,j):Hrow(1,k,j)+1,:);
                [cos,sin,~]=CT_GIV(HR(Hrow(1,k,j)+1,Hrow(1,k,j)+1),HR(Hrow(1,k,j)+1,Hrow(1,k,j)));
                Grot = [cos; sin];
                HR(1:Hrow(1,k,j)+1,Hrow(1,k,j):Hrow(1,k,j)+1) = HR(1:Hrow(1,k,j)+1,Hrow(1,k,j):Hrow(1,k,j)+1) * CT_TO_MAT(Grot);
                Hrot(:,k,j) = Grot;
            end
        end
        CTSWc = [Hrot(:,bs*(bs-1)/2:-1:1,m); Hrow(1,bs*(bs-1)/2:-1:1,m)];
        CTSWc = [CTSWc [reshape(Hrot(:,bs^2:-1:1,m-1:-1:i),2,bs^2*(m-i)); reshape(Hrow(1,bs^2:-1:1,m-1:-1:i),1,bs^2*(m-i))]];
        CTSW = [CTSW CTSWc];
        CTSV = [CTSV CTSWc(:,1:end-bs^2)];
        for j=i+1:m
            Hrot(:,:,j) = CT_H(Hrot(:,:,j));
        end
        % reverse order in i block
        Hrot(:,:,i) = fliplr(Hrot(:,:,i));
        Hrow(:,:,i) = fliplr(Hrow(:,:,i));
    end
end

% Apply to V
Q = eye(m*bs,m*bs);
for i = size(CTSV,2):-1:1
    Q(CTSV(3,i):CTSV(3,i)+1,:) = CT_TO_MAT(CTSV(1:2,i)) * Q(CTSV(3,i):CTSV(3,i)+1,:);
end
V = V*Q;

% Apply to residual
for i = 1:size(CTSW,2)
    E(:,CTSW(3,i):CTSW(3,i)+1) = E(:,CTSW(3,i):CTSW(3,i)+1) * CT_TO_MAT(CTSW(1:2,i));
end

LR = HR;
KLrot = Hrot;
KLrow = Hrow;

end

