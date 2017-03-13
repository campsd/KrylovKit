% CTIR - Core-Transformations Implicit Filter for Extended Krylov
%
% INPUT
% V      Krylov Basis
% KLrot  n-1 c,s Givens rotations K and L Hessenbergs
% KLidx  side of Givens rotations 0 = K, 1 = L
% KR     Upper triangular matrix K Hessenberg
% LR     Upper triangular matrix L Hessenberg
% mu     shifts (applied as single shift steps)
% s      determines the side at which the last rotation is appended
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
% October 4, 2016
% TODO: How to complete the chase??

% -------------------------------------------------------------------------
% Main function
% -------------------------------------------------------------------------
function [V,KLrot,KLidx,KR,LR] = CTIF(A,V,KLrot,KLidx,KR,LR,mu,ss)
	% -----------------------------------------------------------------
    % Extend the pencil
    % -----------------------------------------------------------------
    k = size(KR,1);
    N = null([V A*V]);
    tl = N(1:k,1);
    tk = (A*V)\(V*tl);
    % test
    norm(A*V*tk - V*tl,'fro')
    for i=length(KLrot):-1:1
       if(KLidx(i) == 0) %K
           tk(i:i+1) = CreateRotMat(RotH(KLrot(:,i))) * tk(i:i+1);
       else %L
           tl(i:i+1) = CreateRotMat(RotH(KLrot(:,i))) * tl(i:i+1);
       end
    end
    
%     % still apply the rotations !
%     for i=1:nbk
%         tk(KGi(i,1):KGi(i,2)) = CreateRotMat(RotH(KGr(:,i))) * tk(KGi(i,1):KGi(i,2));
%     end
%     for i=1:nbl
%         tl(LGi(i,1):LGi(i,2)) = CreateRotMat(RotH(LGr(:,i))) * tl(LGi(i,1):LGi(i,2));
%     end
    KR = [KR, tk];    % add a column to make KR square
    LR = [LR, tl];    % add a column to make LR square
	
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
        
        if ss(kk) == KLidx(n-1), %we keep it at the same side
            if KLidx(n-1) == 0, % that is the K side
                [ChaRot, LR] = ShiftRotLeftToRight(RotH(ChaRot), LR, n-1, n-1);
                [ChaRot, KR] = ShiftRotRightToLeft(ChaRot, KR, n-1, n-1);
            else % that side is the L side
                [ChaRot, KR] = ShiftRotLeftToRight(RotH(ChaRot), KR, n-1, n-1);
                [ChaRot, LR] = ShiftRotRightToLeft(ChaRot, LR, n-1, n-1);
            end
                KLrot(:,n-1) = RotFUS(KLrot(:,n-1),RotH(ChaRot));
        else %we keep it at the other side
            if KLidx(n-1) == 0, % that is the L side
                [TmpRot, KR] = ShiftRotLeftToRight(KLrot(:,n-1), KR, n-1, n-1);
                [TmpRot, LR] = ShiftRotRightToLeft(TmpRot, LR, n-1, n-1);
            else % that is the K side
                [TmpRot, LR] = ShiftRotLeftToRight(KLrot(:,n-1), LR, n-1, n-1);
                [TmpRot, KR] = ShiftRotRightToLeft(TmpRot, KR, n-1, n-1);
            end
            KLidx(n-1) = s(kk);
            KLrot(:,n-1) = RotFUS(RotH(ChaRot),RotH(TmpRot));
        end
        
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
