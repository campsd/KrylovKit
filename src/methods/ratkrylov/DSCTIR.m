function [ V, KR, Lrot, LR ] = DSCTIR( V, KR, Lrot, LR, mu )
%DSCTIR - Double Shift Core-Transformations Implicit Restart for rational
%Krylov
%
% INPUT
% V     rational Krylov basis
% KR    upper triangular at K side
% Lrot  n-1 rotations at L side
% LR    upper triangular at L side
% mu    2 by k array applied as double shifts
%
% OUTPUT
% V     Adjusted rational Krylov basis
% KR    upper triangular of adjusted K
% Lrot  rotations of adjusted L
% LR    upper triangular of adjusted L

    for kk=1:size(mu,2)
       n=size(KR,1);
       DescBranch = zeros(2,3);
       Limb = zeros(2,1);
       % -----------------------------------------------------------------
       % Determine initial perturbing rotations of shifted L - mu*K
       % -----------------------------------------------------------------
       R = LR(1:3,1:3);
       R(2:3,1) = CreateRotMat(Lrot(:,2)) * R(2:3,1);
       R(1:2,1) = CreateRotMat(Lrot(:,1)) * R(1:2,1);
       R = (R - mu(2,kk)*KR(1:3,1:3))*(R - mu(1,kk)*KR(1:3,1:3));
       R = R(:,1);
       if max(abs(imag(R)) < 100*eps)
            R = real(R);
        end
       [c,s,r] = RotGIV(R(2),R(3));
       DescBranch(:,2) = [c;s];
       [c,s,~] = RotGIV(R(1),r);
       DescBranch(:,1) = [c;s];
       AscBranch = DescBranch(:,1:2);
       
        % At this moment AscBranch is equal to DescBranch, so at this
        % moment the AscBranch is actually descending (order:[1,2])
        % -----------------------------------------------------------------
        % Step 0
        % -----------------------------------------------------------------
		V = ApplyRotToV(RotH(DescBranch(:,2)), V, 2);
        V = ApplyRotToV(RotH(DescBranch(:,1)), V, 1);

        % LLL
        % Shift the AscBranch through (KR,LR) - this makes the
        % AscBranch truly ascending (order:[2,1])
        [AscBranch, KR, LR] = ShiftBranchThroughPencil(AscBranch, KR, LR, 1, n-1);
        [Limb, DescBranch(:,1),DescBranch(:,2)] = RotST(DescBranch(:,1),DescBranch(:,2), Lrot(:,1));
        DescBranch(:,2) = RotFUS(DescBranch(:,2), Lrot(:,2));
        % ready to augment
        
        % -----------------------------------------------------------------
        % Main loop
        % -----------------------------------------------------------------	
        for i=1:n-4
            
            % Augment descending branch
            DescBranch(:,3) = Lrot(:,i+2);
            % Translate Ascending branch through Descending branch
            [AscBranch(:,2), DescBranch(:,2),DescBranch(:,3)] = RotST(DescBranch(:,2),DescBranch(:,3), AscBranch(:,2));
            [AscBranch(:,1), DescBranch(:,1),DescBranch(:,2)] = RotST(DescBranch(:,1),DescBranch(:,2), AscBranch(:,1));
            % Move limb through Ascending branch to the middle
            [AscBranch(:,2), AscBranch(:,1),Limb] = RotST(Limb,AscBranch(:,2), AscBranch(:,1));
            % Fix CT
            Lrot(:,i) = DescBranch(:,1);
            % Remove Ascending branch from the left
            V = ApplyRotToV(AscBranch(:,2), V, i+2);
            V = ApplyRotToV(AscBranch(:,1), V, i+1);
            % The AscBranch becomes a descending branch (order:[1,2]) at
            % the other side
            
            % We always stay working at L side
            % At this moment the Ascending branch is again
            % descending at the other side
            AscBranch(:,1) = RotH(AscBranch(:,1));
            AscBranch(:,2) = RotH(AscBranch(:,2));
            [AscBranch, KR, LR] = ShiftBranchThroughPencil(AscBranch, KR, LR, i+1, n-1);
            % Do not forget to shift both branches !!
            DescBranch(:,1:2) = DescBranch(:,2:3);
            % ready to augment
               
        end
        
        % -----------------------------------------------------------------
        % Final rotation
        % -----------------------------------------------------------------
		i = n-3;
        % Augment descending branch
        DescBranch(:,3) = Lrot(:,i+2);
        % Translate Ascending branch through Descending branch
        [AscBranch(:,2), DescBranch(:,2),DescBranch(:,3)] = RotST(DescBranch(:,2),DescBranch(:,3), AscBranch(:,2));
        [AscBranch(:,1), DescBranch(:,1),DescBranch(:,2)] = RotST(DescBranch(:,1),DescBranch(:,2), AscBranch(:,1));
        % Move limb through Ascending branch to the middle
        [AscBranch(:,2), AscBranch(:,1),Limb] = RotST(Limb,AscBranch(:,2), AscBranch(:,1));
        % Fix CT
        Lrot(:,i) = DescBranch(:,1);
        % Remove Ascending branch from the left
        V = ApplyRotToV(AscBranch(:,2), V, i+2);
        V = ApplyRotToV(AscBranch(:,1), V, i+1);
        
        % -----------------------------------------------------------------
        % Shrink the Krylov recurrence
        % -----------------------------------------------------------------
        % 1st
		V(:,end) = [];
		KR(:,end) = []; KR(end,:) = [];
		LR(:,end) = []; LR(end,:) = [];
		Lrot(:,end) = [];
        % 2nd
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

function [DBranch, Rstart, Rend] = ShiftBranchThroughPencil(DBranch, Rstart, Rend, i, k)
    % Shifts a Descending branch of two core transformations through the
    % pencil from the Rstart side to the Rend side. The core
    % transformations act on rows i,i+1,i+2. The result is an ascending
    % branch at the Rend side
    [DBranch(:,2), Rstart] = ShiftRotLeftToRight(DBranch(:,2), Rstart, i+1, k);
    [DBranch(:,2), Rend] = ShiftRotRightToLeft(DBranch(:,2), Rend, i+1, k);
    DBranch(:,2) = RotH(DBranch(:,2)) ;
    [DBranch(:,1), Rstart] = ShiftRotLeftToRight(DBranch(:,1), Rstart, i, k);
    [DBranch(:,1), Rend] = ShiftRotRightToLeft(DBranch(:,1), Rend, i, k);
    DBranch(:,1) = RotH(DBranch(:,1)) ;
end

