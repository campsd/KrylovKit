function [ V, KLrot, KLidx, KR, LR, r, em ] = CT_EK_TO_EK_RIGHT( V, KLrot, KLidx, KR, LR, s )
%[ V, KLrot, KLidx, KR, LR, r, em ] = CT_EK_TO_EK_RIGHT( V, KLrot, KLidx, KR, LR, s )
%-- converts the selection vector of an EK recurrence to another selection
%vector by transforming the pattern of core transformations from the right
%
% INPUT
% V         extended Krylov basis (N x m+1)
% KLrot     core transformations L,K pencil (2xm)
% KLidx     indicates the side of the core transformation (0 = K // 1 = L)
%           (m)
% KR        upper triangular Hessenberg K (m+1 x m)
% LR        upper triangular Hessenberg L (m+1 x m)
% s         desired selection vector (m)
% Recurrence that holds
% 	A * V * mat(KLrot(KLidx==0)) * KR = 
%	V * mat(KLrot(KLidx==1)) * LR
%
% OUTPUT
% V         adjusted extended Krylov basis (N x m)
% KLrot     adjusted core transformations L,K pencil (m-1)
% KLidx     indicates the side on which the core transformation acts (0 = K
%           // 1 = L) This satisfies the input s (m-1)
% KR        upper triangular Hessenberg K (m x m)
% LR        upper triangular Hessenberg L (m x m)
% r         residual vector (N x 1)
% en        adjusted mth canonical vector (1 x m)
% Recurrence that holds:
% 	A * V * mat(KLrot(KLidx==0)) * KR = 
%	V * mat(KLrot(KLidx==1)) * LR + (A*)r * em
% The multiplication with A is required when initially KLidx(m)==0
%
% daan.camps@cs.kuleuven.be
% last edit: April 14, 2017
%
% See also: CT_SK, CT_EK, CT_SK_TO_EK_RIGHT, CT_EK_TO_EK_LEFT

m = size(KLrot,2);
CTSV = zeros(3,0);
CTSW = zeros(3,0);

% Explicitly extract the residual from the recurrence
if KLidx(m) == 1 %residual at L
   LR(m:m+1,m) = CT_TO_MAT(KLrot(:,m)) *  LR(m:m+1,m);
   lmp1m = LR(m+1,m); vmp1 = V(:,m+1); em = transpose(flipud(eye(m,1)));
   r = lmp1m*vmp1;
else %residual at K
   KR(m:m+1,m) = CT_TO_MAT(KLrot(:,m)) *  KR(m:m+1,m);
   kmp1m = KR(m+1,m); vmp1 = V(:,m+1); em = transpose(flipud(eye(m,1)));
   r = -kmp1m * vmp1;
end
KR(m+1,:) = []; LR(m+1,:) = []; V(:,m+1) = [];
KLrot(:,m) = []; KLidx(m) = [];

% Modify the top square part of the pencil to satisfy the new selection
% vector
for i=1:m-1
    if ~((s(i) == 1 && KLidx(i) == 1) || (s(i) == -1 && KLidx(i) == 0))
       % the core transformation needs to be brought to the other side       
       if s(i) == 1 %ith needs to move from K to L
           % move all L CTs from i+1 to m-1 to K via V
           CTSV = [CTSV [KLrot(:,logical([zeros(1,i) (KLidx(i+1:m-1)==1)])); find(KLidx(i+1:m-1)==1)+i]];
           KLrot(:,logical([zeros(1,i) (KLidx(i+1:m-1)==1)])) = CT_H(KLrot(:,logical([zeros(1,i) (KLidx(i+1:m-1)==1)])));
           % shift all K CTs from i to m-1 to L via pencil
           for j=m-1:-1:i
              if KLidx(j) == 0
                 Grot = KLrot(:,j);
                 ShiftRotLeftKToRightK(j);
                 CTSW = [CTSW [Grot;j]];
                 ShiftRotRightLToLeftL(j);
                 KLrot(:,j) = CT_H(Grot);
              end
           end
           % the ith is at the correct place now
           KLidx(i) = 1;
           % move the other CTs back to K via V
           CTSVc = [KLrot(:,logical([zeros(1,i) (KLidx(i+1:m-1)==0)])); find(KLidx(i+1:m-1)==0)+i];
           CTSV = [CTSV fliplr(CTSVc) ];
           KLrot(:,logical([zeros(1,i) (KLidx(i+1:m-1)==0)])) = CT_H(KLrot(:,logical([zeros(1,i) (KLidx(i+1:m-1)==0)])));
           % and back to L via pencil
           for j=i+1:m-1
              if KLidx(j) == 1
                 Grot = KLrot(:,j);
                 ShiftRotLeftKToRightK(j);
                 CTSW = [CTSW [Grot;j]];
                 ShiftRotRightLToLeftL(j);
                 KLrot(:,j) = CT_H(Grot);
              end
           end
       else %ith needs to move from L to K
           % move all K CTs from i+1 to m-1 to L via V
           CTSV = [CTSV [KLrot(:,logical([zeros(1,i) (KLidx(i+1:m-1)==0)])); find(KLidx(i+1:m-1)==0)+i]];
           KLrot(:,logical([zeros(1,i) (KLidx(i+1:m-1)==0)])) = CT_H(KLrot(:,logical([zeros(1,i) (KLidx(i+1:m-1)==0)])));
           % shift all L CTs from i to m-1 to K via pencil
           for j=m-1:-1:i
              if KLidx(j) == 1
                 Grot = KLrot(:,j);
                 ShiftRotLeftLToRightL(j); % now,CT_H(Grot) is at the right of L
                 CTSW = [CTSW [Grot;j]]; % remove from the right
                 ShiftRotRightKToLeftK(j); % bring to left of K
                 KLrot(:,j) = CT_H(Grot); % now, CT_H(Grot) is at K
              end
           end
           % the ith is at the correct place now
           KLidx(i) = 0; %KLrot(:,i) = CT_H(KLrot(:,i));
           % move the other CTs back to L via V
           CTSVc = [KLrot(:,logical([zeros(1,i) (KLidx(i+1:m-1)==1)])); find(KLidx(i+1:m-1)==1)+i];
           CTSV = [CTSV fliplr(CTSVc) ];
           KLrot(:,logical([zeros(1,i) (KLidx(i+1:m-1)==1)])) = CT_H(KLrot(:,logical([zeros(1,i) (KLidx(i+1:m-1)==1)])));
           % and back to K via pencil
           for j=i+1:m-1
              if KLidx(j) == 0
                 Grot = KLrot(:,j);
                 ShiftRotLeftLToRightL(j);
                 CTSW = [CTSW [Grot;j]];
                 ShiftRotRightKToLeftK(j);
                 KLrot(:,j) = CT_H(Grot); %?
              end
           end
       end
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

% Nested functions
% -------------------------------------------------------------------------
function ShiftRotLeftLToRightL(i)
% 	% Shifts a rotation from left to right through the upper LR triangular
	LR(i:i+1,i:m) = CT_TO_MAT(Grot)*LR(i:i+1,i:m);
	[cos,sin,~]=CT_GIV(LR(i+1,i+1),LR(i+1,i));
	Grot = [cos;sin];
	LR(1:i+1,i:i+1) = LR(1:i+1,i:i+1)*CT_TO_MAT(Grot);
end

function ShiftRotLeftKToRightK(i)
	% Shifts a rotation from left to right through the upper LR triangular
	KR(i:i+1,i:m) = CT_TO_MAT(Grot)*KR(i:i+1,i:m);
	[cos,sin,~]=CT_GIV(KR(i+1,i+1),KR(i+1,i));
	Grot = [cos;sin];
	KR(1:i+1,i:i+1) = KR(1:i+1,i:i+1)*CT_TO_MAT(Grot);
end

function ShiftRotRightKToLeftK(i)
	% Shifts a rotation from right to left through the upper triangular
	KR(1:i+1,i:i+1) = KR(1:i+1,i:i+1)*CT_TO_MAT(Grot);
	[cos,sin,~]=CT_GIV(KR(i,i),KR(i+1,i));
	Grot = [cos;sin];
	KR(i:i+1,i:m) = CT_TO_MAT(Grot)*KR(i:i+1,i:m);
end

function ShiftRotRightLToLeftL(i)
	% Shifts a rotation from right to left through the upper triangular
	LR(1:i+1,i:i+1) = LR(1:i+1,i:i+1)*CT_TO_MAT(Grot);
	[cos,sin,~]=CT_GIV(LR(i,i),LR(i+1,i));
	Grot = [cos;sin];
	LR(i:i+1,i:m) = CT_TO_MAT(Grot)*LR(i:i+1,i:m);
end

end

