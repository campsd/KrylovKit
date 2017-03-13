% CTIR - Core-Transformations Implicit Restart for Extended Krylov
%
% INPUT
% V      Krylov Basis
% KLrot  n-1 c,s Givens rotations K and L Hessenbergs
% KLidx  side of Givens rotations 0 = K, 1 = L
% KR     Upper triangular matrix K Hessenberg
% LR     Upper triangular matrix L Hessenberg
% mu     shifts (applied as single shift steps)
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
% September 30, 2016

% -------------------------------------------------------------------------
% Main function
% -------------------------------------------------------------------------
function [V,KLrot,KLidx,KR,LR] = CTIR(V,KLrot,KLidx,KR,LR,mu)
	
	
	for kk = 1:length(mu)
		n = size(KR,1);
        % -----------------------------------------------------------------
        % Determine initial perturbing rotation of shifted L - mu*K
    	% -----------------------------------------------------------------
		if (KLidx(1) == 0) % First rotation at K (system solve)
			R = LR(1:2,1) - mu(kk)*CreateRotMat(KLrot(:,1))*KR(1:2,1);
		elseif (KLidx(1) == 1) % First rotation at L (multiplication)
			R = CreateRotMat(KLrot(:,1))*LR(1:2,1) - [mu(kk)*KR(1,1);0];
		end
		[c,s,~] = RotGIV(R(1),R(2));
		ChaRot = [c;s];

        % -----------------------------------------------------------------
        % Step 0
        % -----------------------------------------------------------------
		V = ApplyRotToV(RotH(ChaRot), V, 1);
		KLrot(:,1) = RotFUS(ChaRot, KLrot(:,1));
        if (KLidx(1) == 0) % First rotation at K
            if (KLidx(2) == 0) % Second rotation also at K
                %       K       L          K     L
                %fused->|  cha->|    -->   | |
                %        |                  |	
                [ChaRot, LR] = ShiftRotLeftToRight(ChaRot, LR, 1, n-1);
                [ChaRot, KR] = ShiftRotRightToLeft(ChaRot, KR, 1, n-1);
                ChaRot = RotH(ChaRot) ;
            else % Second rotation swapped to L
                %       K       L          K      L
                %fused->|  cha->|    -->		  | |
                %                |			       |
                tmp = ChaRot;
                ChaRot = KLrot(:,1);
                KLrot(:,1) = tmp;
                KLidx(1) = 1;

                [ChaRot, KR] = ShiftRotLeftToRight(ChaRot, KR, 1, n-1); % RotH(ChaRot) lives at right of K, annihilated by ChaRot
                [ChaRot, LR] = ShiftRotRightToLeft(ChaRot, LR, 1, n-1); 
                ChaRot = RotH(ChaRot) ;
            end
		elseif (KLidx(1) == 1) % First rotation at L
			if (KLidx(2) == 1) % Second rotation also at L
				%       K     L           K     L
				%  cha->|fus->|   -->           | |
				%              |                 |
				[ChaRot, KR] = ShiftRotLeftToRight(ChaRot, KR, 1, n-1);
				[ChaRot, LR] = ShiftRotRightToLeft(ChaRot, LR, 1, n-1);
                ChaRot = RotH(ChaRot) ;
			else % Second rotation swapped to K
				%        K     L          K       L
				%   fus->|cha->|	-->	  | |
				%         |                |
				tmp = ChaRot;
				ChaRot = KLrot(:,1);
				KLrot(:,1) = tmp;
				KLidx(1) = 0;

				[ChaRot, LR] = ShiftRotLeftToRight(ChaRot, LR, 1, n-1);
				[ChaRot, KR] = ShiftRotRightToLeft(ChaRot, KR, 1, n-1);
                ChaRot = RotH(ChaRot) ;
			end
        end
        
        % -----------------------------------------------------------------
        % Main loop
        % -----------------------------------------------------------------	
        for i=1:n-3
            if ( KLidx(i) == 0 ) % ith rotation at K side
                if ( KLidx(i+1) == 0 ) % i+1st rotation also at K side
                    [ChaRot, KLrot(:,i),KLrot(:,i+1)] = RotST(KLrot(:,i),KLrot(:,i+1), ChaRot) ;
                    V = ApplyRotToV(ChaRot, V, i+1);
                    if ( KLidx(i+2) == 0 ) % i+2nd rotation also at K Side
                        %	K	L
                        %	|
                        %	 |
                        %	  |
                        [ChaRot, LR] = ShiftRotLeftToRight(RotH(ChaRot), LR, i+1, n-1);
                        [ChaRot, KR] = ShiftRotRightToLeft(ChaRot, KR, i+1, n-1);
                        ChaRot = RotH(ChaRot) ;
                    else % i+2nd rotation swapped to L side				
                        %	K	L
                        %	|
                        %	 |
                        %		|
                        tmp = RotH(ChaRot);
                        [ChaRot, KR] = ShiftRotLeftToRight(KLrot(:,i+1), KR, i+1, n-1);
                        [ChaRot, LR] = ShiftRotRightToLeft(ChaRot, LR, i+1, n-1);
                        ChaRot = RotH(ChaRot) ;
                        KLrot(:,i+1) = tmp;
                        KLidx(i+1) = 1;
                    end
                end
            elseif ( KLidx(i) == 1 ) % ith rotation at L side
                if ( KLidx(i+1) == 1 ) % i+1st rotation also at L side
                    [ChaRot, KLrot(:,i),KLrot(:,i+1)] = RotST(KLrot(:,i),KLrot(:,i+1), ChaRot) ;
                    V = ApplyRotToV(ChaRot, V, i+1);
                    if ( KLidx(i+2) == 1 ) % i+2nd rotation also at L side
                        %	K	L
                        %		|
                        %		 |
                        %		  |		
                        [ChaRot, KR] = ShiftRotLeftToRight(RotH(ChaRot), KR, i+1, n-1);
                        [ChaRot, LR] = ShiftRotRightToLeft(ChaRot, LR, i+1, n-1);
                        ChaRot = RotH(ChaRot) ;
                    else % i+2ndrotation swapped to K side
                        %	K	L
                        %		|
                        %		 |
                        %	|
                        tmp = RotH(ChaRot);
                        [ChaRot, LR] = ShiftRotLeftToRight(KLrot(:,i+1), LR, i+1, n-1);
                        [ChaRot, KR] = ShiftRotRightToLeft(ChaRot, KR, i+1, n-1);
                        ChaRot = RotH(ChaRot) ;
                        KLrot(:,i+1) = tmp;
                        KLidx(i+1) = 0;
                    end
                end
            end
        end
        
        % -----------------------------------------------------------------
        % Final rotation
        % -----------------------------------------------------------------
		i = n-2;
		[ChaRot, KLrot(:,i),KLrot(:,i+1)] = RotST(KLrot(:,i),KLrot(:,i+1), ChaRot) ;
		V = ApplyRotToV(ChaRot, V, i+1);
        
        % -----------------------------------------------------------------
        % Shrink the Krylov recurrence
        % -----------------------------------------------------------------
		V(:,end) = [];
		KR(:,end) = []; KR(end,:) = [];
		LR(:,end) = []; LR(end,:) = [];
		KLrot(:,end) = []; KLidx(end) = [];
        
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
