function [V,KLrot,KLidx,KR,LR] = CT_EK_IR_DS(V,KLrot,KLidx,KR,LR,mu)
% CT_EK_IR_DS -- Double Shift Core-Transformations Implicit Restart for 
% Extended Krylov
%
% INPUT
% V      Krylov Basis
% KLrot  n-1 c,s Givens rotations K and L Hessenbergs
% KLidx  side of Givens rotations 0 = K, 1 = L
% KR     Upper triangular matrix K Hessenberg
% LR     Upper triangular matrix L Hessenberg
% mu     2 by k array applied as double shifts
%
% OUTPUT
% V     Adjusted Krylov basis of dim-len(mu)
% KLrot  c,s Givens rotations adjusted K Hessenberg
% KLidx  order Givens rotations adjusted K Hessenberg
% KR     Upper triangular matrix adjusted K Hessenberg
% LR     Upper triangular matrix adjusted L Hessenberg
%
% The recursion relation that holds before and after the restart is:
%       A*V*K = V*L
%
%
% daan.camps@cs.kuleuven.be
% November 22, 2016

% -------------------------------------------------------------------------
% Main function
% -------------------------------------------------------------------------

	for kk = 1:size(mu,2)
		n = size(KR,1);
        DescBranch = zeros(2,3);
        Limb = zeros(2,1);
        % -----------------------------------------------------------------
        % Determine initial perturbing rotations of shifted (L -
        % mu(1,kk)*K)*(L - mu(2,kk)*K)
    	% -----------------------------------------------------------------

        % This should always be valid
        [tmpK,tmpL] = CT_EK_PENCIL(KLrot, KLidx, KR, LR);
        tmpK(end,:) = []; tmpL(end,:) = [];
        R = eye(size(tmpK,1),1);
        for k = 1:2
            if (KLidx(k) == 1) % Hessenberg
                R = (tmpL*(tmpK\R)) - mu(k,kk)*R;
            else % Inverse Hessenberg
                R = R - mu(k,kk)*tmpK*(tmpL\R); % This can cause issues with R becoming singular (cavity)
            end
        end
        if max(abs(imag(R)) < 100*eps)
            R = real(R);
        end
        [c,s,r] = CT_GIV(R(2),R(3));
		DescBranch(:,2) = [c;s];
        [c,s,~] = CT_GIV(R(1),r);
        DescBranch(:,1) = [c;s];
        AscBranch = DescBranch(:,1:2);
        % At this moment AscBranch is equal to DescBranch, so at this
        % moment the AscBranch is actually descending (order:[1,2])
        % -----------------------------------------------------------------
        % Step 0
        % -----------------------------------------------------------------
		V = ApplyRotToV(CT_H(DescBranch(:,2)), V, 2);
        V = ApplyRotToV(CT_H(DescBranch(:,1)), V, 1);

        if (KLidx(1) == 0) % Kxx
            if (KLidx(2) == 0) % KKx
                if (KLidx(3) == 0) % KKK
                    % Shift the AscBranch through (LR,KR) - this makes the
                    % AscBranch truly ascending (order:[2,1])
                    [AscBranch, LR, KR] = ShiftBranchThroughPencil(AscBranch, LR, KR, 1, n-1);
                    [Limb, DescBranch(:,1),DescBranch(:,2)] = CT_TURNOVER(DescBranch(:,1),DescBranch(:,2), KLrot(:,1));
                    LimbSide = 0;
                    DescBranch(:,2) = CT_FUSE(DescBranch(:,2), KLrot(:,2));
                    % ready to augment
                else % KKL
                    [KLrot(:,1:2), KR, LR] = ShiftBranchThroughPencil(KLrot(:,1:2), KR, LR, 1, n-1);
                    [AscBranch, KR, LR] = ShiftBranchThroughPencil(AscBranch, KR, LR, 1, n-1);
                    [KLrot(:,1), KLrot(:,2),AscBranch(:,2)] = CT_TURNOVER(KLrot(:,2),KLrot(:,1), AscBranch(:,2));
                    AscBranch(:,1) = CT_FUSE(AscBranch(:,2), AscBranch(:,1));
                    AscBranch(:,2) = KLrot(:,2);
                    [Limb, DescBranch(:,1),DescBranch(:,2)] = CT_TURNOVER(DescBranch(:,1),DescBranch(:,2), KLrot(:,1));
                    LimbSide = 1;
                    % ready to augment
                end
            else % KLx
                if (KLidx(3) == 0) % KLK
                    AscBranch(:,2) = CT_FUSE(AscBranch(:,2), KLrot(:,2));
                    [AscBranch, LR, KR] = ShiftBranchThroughPencil(AscBranch, LR, KR, 1, n-1);
                    [Limb, DescBranch(:,1),DescBranch(:,2)] = CT_TURNOVER(DescBranch(:,1),DescBranch(:,2), KLrot(:,1));
                    LimbSide = 0;
                    % ready to augment
                else % KLL
                    [Limb, DescBranch(:,1),DescBranch(:,2)] = CT_TURNOVER(DescBranch(:,1),DescBranch(:,2), KLrot(:,1));
                    V = ApplyRotToV(Limb, V, 2);
                    Limb = CT_H(Limb);
                    [DescBranch(:,1:2), KR, LR] = ShiftBranchThroughPencil(DescBranch(:,1:2), KR, LR, 1, n-1);
                    AscBranch(:,2) = CT_FUSE(AscBranch(:,2), KLrot(:,2));
                    tmpBranch = AscBranch; AscBranch = DescBranch(:,1:2); DescBranch(:,1:2) = tmpBranch;
                    LimbSide = 1;
                    % ready to augment
                end
            end
            %
		elseif (KLidx(1) == 1) % Lxx
            if (KLidx(2) == 0) % LKx
                if (KLidx(3) == 0) % LKK
                    [Limb, DescBranch(:,1),DescBranch(:,2)] = CT_TURNOVER(DescBranch(:,1),DescBranch(:,2), KLrot(:,1));
                    V = ApplyRotToV(Limb, V, 2);
                    Limb = CT_H(Limb);
                    [DescBranch(:,1:2), LR, KR] = ShiftBranchThroughPencil(DescBranch(:,1:2), LR, KR, 1, n-1);
                    AscBranch(:,2) = CT_FUSE(AscBranch(:,2), KLrot(:,2));
                    tmpBranch = AscBranch; AscBranch = DescBranch(:,1:2); DescBranch(:,1:2) = tmpBranch;
                    LimbSide = 0;
                    % ready to augment
                else % LKL
                    AscBranch(:,2) = CT_FUSE(AscBranch(:,2), KLrot(:,2));
                    [AscBranch, KR, LR] = ShiftBranchThroughPencil(AscBranch, KR, LR, 1, n-1);
                    [Limb, DescBranch(:,1),DescBranch(:,2)] = CT_TURNOVER(DescBranch(:,1),DescBranch(:,2), KLrot(:,1));
                    LimbSide = 1;
                    % ready to augment
                end
            else % LLx
                if (KLidx(3) == 0) % LLK
                    [KLrot(:,1:2), LR, KR] = ShiftBranchThroughPencil(KLrot(:,1:2), LR, KR, 1, n-1);
                    [AscBranch, LR, KR] = ShiftBranchThroughPencil(AscBranch, LR, KR, 1, n-1);
                    [KLrot(:,1), KLrot(:,2),AscBranch(:,2)] = CT_TURNOVER(KLrot(:,2),KLrot(:,1), AscBranch(:,2));
                    AscBranch(:,1) = CT_FUSE(AscBranch(:,2), AscBranch(:,1));
                    AscBranch(:,2) = KLrot(:,2);
                    [Limb, DescBranch(:,1),DescBranch(:,2)] = CT_TURNOVER(DescBranch(:,1),DescBranch(:,2), KLrot(:,1));
                    LimbSide = 0;
                    % ready to augment
                else % LLL
                    % Shift the AscBranch through (KR,LR) - this makes the
                    % AscBranch truly ascending (order:[2,1])
                    [AscBranch, KR, LR] = ShiftBranchThroughPencil(AscBranch, KR, LR, 1, n-1);
                    [Limb, DescBranch(:,1),DescBranch(:,2)] = CT_TURNOVER(DescBranch(:,1),DescBranch(:,2), KLrot(:,1));
                    LimbSide = 1;
                    DescBranch(:,2) = CT_FUSE(DescBranch(:,2), KLrot(:,2));
                    % ready to augment
                end
            end
            %
        end
        
        % -----------------------------------------------------------------
        % Main loop
        % -----------------------------------------------------------------	
        for i=1:n-4
            
            % Augment descending branch
            DescBranch(:,3) = KLrot(:,i+2);
            % Translate Ascending branch through Descending branch
            [AscBranch(:,2), DescBranch(:,2),DescBranch(:,3)] = CT_TURNOVER(DescBranch(:,2),DescBranch(:,3), AscBranch(:,2));
            [AscBranch(:,1), DescBranch(:,1),DescBranch(:,2)] = CT_TURNOVER(DescBranch(:,1),DescBranch(:,2), AscBranch(:,1));
            % Move limb through Ascending branch to the middle
            [AscBranch(:,2), AscBranch(:,1),Limb] = CT_TURNOVER(Limb,AscBranch(:,2), AscBranch(:,1));
            % Fix CT
            KLrot(:,i) = DescBranch(:,1);
            KLidx(i) = LimbSide;
            % Remove Ascending branch from the left
            V = ApplyRotToV(AscBranch(:,2), V, i+2);
            V = ApplyRotToV(AscBranch(:,1), V, i+1);
            % The AscBranch becomes a descending branch (order:[1,2]) at
            % the other side
            
            if (LimbSide == 0) % Currently working at K side
               if (KLidx(i+3) == 0) % We stay working at K side
                   AscBranch(:,1) = CT_H(AscBranch(:,1));
                   AscBranch(:,2) = CT_H(AscBranch(:,2));
                   [AscBranch, LR, KR] = ShiftBranchThroughPencil(AscBranch, LR, KR, i+1, n-1);
                   DescBranch(:,1:2) = DescBranch(:,2:3);
                    % ready to augment
               else % We move to L side
                   AscBranch(:,1) = CT_H(AscBranch(:,1));
                   AscBranch(:,2) = CT_H(AscBranch(:,2));
                   tmpBranch = AscBranch;
                   [DescBranch(:,2:3), KR, LR] = ShiftBranchThroughPencil(DescBranch(:,2:3), KR, LR, i+1, n-1);
                   AscBranch = DescBranch(:,2:3);
                   DescBranch(:,1:2) = tmpBranch;
                   V = ApplyRotToV(Limb, V, i+2);
                   Limb = CT_H(Limb);
                   LimbSide = 1;
                   % ready to augment
               end
            else % Currently working at L side
                if (KLidx(i+3) == 0) % We move to K side
                   AscBranch(:,1) = CT_H(AscBranch(:,1));
                   AscBranch(:,2) = CT_H(AscBranch(:,2));
                   tmpBranch = AscBranch;
                   [DescBranch(:,2:3), LR, KR] = ShiftBranchThroughPencil(DescBranch(:,2:3), LR, KR, i+1, n-1);
                   AscBranch = DescBranch(:,2:3);
                   DescBranch(:,1:2) = tmpBranch;
                   V = ApplyRotToV(Limb, V, i+2);
                   Limb = CT_H(Limb);
                   LimbSide = 0;
                   % ready to augment
                else % We stay working at L side
                   % At this moment the Ascending branch is again
                   % descending at the other side
                   AscBranch(:,1) = CT_H(AscBranch(:,1));
                   AscBranch(:,2) = CT_H(AscBranch(:,2));
                   [AscBranch, KR, LR] = ShiftBranchThroughPencil(AscBranch, KR, LR, i+1, n-1);
                   % Do not forget to shift both branches !!
                   DescBranch(:,1:2) = DescBranch(:,2:3);
                   % ready to augment 
                end
            end
        end
        %break
        
        % -----------------------------------------------------------------
        % Final rotation
        % -----------------------------------------------------------------
		i = n-3;
        % Augment descending branch
        DescBranch(:,3) = KLrot(:,i+2);
        % Translate Ascending branch through Descending branch
        [AscBranch(:,2), DescBranch(:,2),DescBranch(:,3)] = CT_TURNOVER(DescBranch(:,2),DescBranch(:,3), AscBranch(:,2));
        [AscBranch(:,1), DescBranch(:,1),DescBranch(:,2)] = CT_TURNOVER(DescBranch(:,1),DescBranch(:,2), AscBranch(:,1));
        % Move limb through Ascending branch to the middle
        [AscBranch(:,2), AscBranch(:,1),Limb] = CT_TURNOVER(Limb,AscBranch(:,2), AscBranch(:,1));
        % Fix CT
        KLrot(:,i) = DescBranch(:,1);
        KLidx(i) = LimbSide;
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
		KLrot(:,end) = []; KLidx(end) = [];
        % 2nd
        V(:,end) = [];
		KR(:,end) = []; KR(end,:) = [];
		LR(:,end) = []; LR(end,:) = [];
		KLrot(:,end) = []; KLidx(end) = [];
        end    
	end
%end

function [V] = ApplyRotToV(Grot, V, i)
	% Applies a rotation to the Krylov basis
	V(:,i:i+1) = V(:,i:i+1)*CT_TO_MAT(Grot);
end

function [Grot, R] = ShiftRotLeftToRight(Grot, R, i, k)
	% Shifts a rotation from left to right through the upper triangular
	R(i:i+1,i:k) = CT_TO_MAT(Grot)*R(i:i+1,i:k);
	[c,s,~]=CT_GIV(R(i+1,i+1),R(i+1,i));
	Grot = [c,s];
	R(1:i+1,i:i+1) = R(1:i+1,i:i+1)*CT_TO_MAT(Grot);
end

function [Grot, R] = ShiftRotRightToLeft(Grot, R, i, k)
	% Shifts a rotation from right to left through the upper triangular
	R(1:i+1,i:i+1) = R(1:i+1,i:i+1)*CT_TO_MAT(Grot);
	[c,s,~]=CT_GIV(R(i,i),R(i+1,i));
	Grot = [c,s];
	R(i:i+1,i:k) = CT_TO_MAT(Grot)*R(i:i+1,i:k);
end

function [DBranch, Rstart, Rend] = ShiftBranchThroughPencil(DBranch, Rstart, Rend, i, k)
    % Shifts a Descending branch of two core transformations through the
    % pencil from the Rstart side to the Rend side. The core
    % transformations act on rows i,i+1,i+2. The result is an ascending
    % branch at the Rend side
    [DBranch(:,2), Rstart] = ShiftRotLeftToRight(DBranch(:,2), Rstart, i+1, k);
    [DBranch(:,2), Rend] = ShiftRotRightToLeft(DBranch(:,2), Rend, i+1, k);
    DBranch(:,2) = CT_H(DBranch(:,2)) ;
    [DBranch(:,1), Rstart] = ShiftRotLeftToRight(DBranch(:,1), Rstart, i, k);
    [DBranch(:,1), Rend] = ShiftRotRightToLeft(DBranch(:,1), Rend, i, k);
    DBranch(:,1) = CT_H(DBranch(:,1)) ;
end
