function [ V, KR, Lrot, LR ] = CTIR( V, KR, Lrot, LR, mu )
%CTIR - Core transformations implicit restart for rational Krylov
%
% Executes single shift QZ steps for the implicit restart of the modified
% RKS sequence from the CTRKS algorithm. 
%
% INPUT
% V     Krylov basis
% KR    upper triangular at K side
% Lrot  n-1 rotations at L side
% LR    upper triangular at L side
% mu    shifts (applied as single shift steps)
%
% OUTPUT
% V     Adjusted Krylov basis
% KR    upper triangular of adjusted K
% Lrot  rotations of adjusted L
% LR    upper triangular of adjusted L
%
% The recursion that hold before and after restart is:
%       A*V*KR = B*V*L
%
% daan.camps@cs.kuleuven.be
% January 12, 2017
    
    for kk=1:length(mu)
       n = size(KR,1);
       % -----------------------------------------------------------------
       % Determine initial perturbing rotation of shifted L - mu*K
       % -----------------------------------------------------------------
       R = CreateRotMat(Lrot(:,1)) * LR(1:2,1) - [mu(kk)*KR(1,1);0];
       [c,s,~] = RotGIV(R(1),R(2));
       ChaRot = [c;s];
       
       % -----------------------------------------------------------------
       % Step 0
       % -----------------------------------------------------------------
       V = ApplyRotToV(RotH(ChaRot), V, 1);
       Lrot(:,1) = RotFUS(ChaRot, Lrot(:,1));
       [ChaRot, KR] = ShiftRotLeftToRight(ChaRot, KR, 1, n-1);
	   [ChaRot, LR] = ShiftRotRightToLeft(ChaRot, LR, 1, n-1);
       ChaRot = RotH(ChaRot);
       
       % -----------------------------------------------------------------
       % Main loop
       % -----------------------------------------------------------------
       for i=1:n-3
           [ChaRot, Lrot(:,i),Lrot(:,i+1)] = RotST(Lrot(:,i),Lrot(:,i+1), ChaRot) ;
           V = ApplyRotToV(ChaRot, V, i+1);
           [ChaRot, KR] = ShiftRotLeftToRight(RotH(ChaRot), KR, i+1, n-1);
           [ChaRot, LR] = ShiftRotRightToLeft(ChaRot, LR, i+1, n-1);
           ChaRot = RotH(ChaRot) ;
       end
       
       % -----------------------------------------------------------------
       % Final rotation
       % -----------------------------------------------------------------
       i = n-2;
	   [ChaRot, Lrot(:,i),Lrot(:,i+1)] = RotST(Lrot(:,i),Lrot(:,i+1), ChaRot) ;
	   V = ApplyRotToV(ChaRot, V, i+1);
       
       % -----------------------------------------------------------------
       % Shrink the Krylov recurrence
       % -----------------------------------------------------------------
	   V(:,end) = [];
	   KR(:,end) = []; KR(end,:) = [];
	   LR(:,end) = []; LR(end,:) = [];
	   Lrot(:,end) = [];
    end

end

function [V] = ApplyRotToV(Grot, V, i)
	% Applies a rotation to the Krylov basis
	V(:,i:i+1) = V(:,i:i+1)*CreateRotMat(Grot);
end

function [Grot, R] = ShiftRotLeftToRight(Grot, R, i, k)
	% Shifts a rotation from left to right through the upper triangular
	R(i:i+1,i:k) = CreateRotMat(Grot)*R(i:i+1,i:k);
	[c,s,~]=RotGIV(R(i+1,i+1),R(i+1,i));
	Grot = [c,s];
	R(1:i+1,i:i+1) = R(1:i+1,i:i+1)*CreateRotMat(Grot);
end

function [Grot, R] = ShiftRotRightToLeft(Grot, R, i, k)
	% Shifts a rotation from right to left through the upper triangular
	R(1:i+1,i:i+1) = R(1:i+1,i:i+1)*CreateRotMat(Grot);
	[c,s,~]=RotGIV(R(i,i),R(i+1,i));
	Grot = [c,s];
	R(i:i+1,i:k) = CreateRotMat(Grot)*R(i:i+1,i:k);
end

