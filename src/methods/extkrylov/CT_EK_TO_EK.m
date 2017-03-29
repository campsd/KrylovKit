function [ V, KLrot, KLidx, KR, LR ] = CT_EK_TO_EK( V, KLrot, KLidx, KR, LR, s )
%CT_EK_TO_EK -- converts the selection vector of an EK recurrence to
%another selection vector
%
% March 28, 2017

% work in progress
    m = length(KLrot);
    CTSV = zeros(3,0);
    
    for i=m:-1:1
        if ~((s(i) == 1 && KLidx(i) == 1) || (s(i) == -1 && KLidx(i) == 0))
            % the core transformation needs to be brought to the other side
            
            if s(i) == 1 %the ith core transformation needs to be brought from K to L
                CTSV = [CTSV [KLrot(:,(KLidx(1:i)==0)); find(KLidx(1:i)==0)']];
                KLrot(:,i) = CT_H(KLrot(:,i));
                KLidx(i) = 1;
                % bring the original rotations at L side through L,K such 
                % that we get the correct order again
                for j=i-1:-1:1
                    if KLidx(j) == 1
                        Grot = KLrot(:,j);
                        ShiftRotLeftLToRightL(j)
                        ShiftRotRightKToLeftK(j)
                        CTSV = [CTSV [CT_H(Grot); j]]; % CT_H maybe vice versa
                        KLrot(:,j) = Grot;
                    end
                end
                % bring the rotations that have moved to the L side back to
                % the K side
                for j=1:i-1
                    if KLidx(j) == 0
                       Grot = KLrot(:,j);
                       ShiftRotLeftLToRightL(j)
                       ShiftRotRightKToLeftK(j)
                       KLrot(:,j) = CT_H(Grot);
                    end
                end
            else %the ith core transformation needs to be brought from L to K
                CTSV = [CTSV [KLrot(:,(KLidx(1:i)==1)); find(KLidx(1:i)==1)']];
                KLrot(:,i) = CT_H(KLrot(:,i));
                KLidx(i) = 0;
                % bring the original rotations at K side through K,L such 
                % that we get the correct order again
                for j=i-1:-1:1
                    if KLidx(j) == 0
                        Grot = KLrot(:,j);
                        ShiftRotLeftKToRightK(j)
                        ShiftRotRightLToLeftL(j)
                        CTSV = [CTSV [CT_H(Grot); j]]; % CT_H maybe vice versa
                        KLrot(:,j) = Grot;
                    end
                end
                % bring the rotations that have moved to the K side back to
                % the L side
                for j=1:i-1
                    if KLidx(j) == 1
                       Grot = KLrot(:,j);
                       ShiftRotLeftKToRightK(j)
                       ShiftRotRightLToLeftL(j)
                       KLrot(:,j) = CT_H(Grot);
                    end
                end
            end
        end
    end
    
    % Apply the rotations to V
    Q = eye(size(KR,1));
    for i=size(CTSV,2):-1:1
        Q(CTSV(3,i):CTSV(3,i)+1,:) = CT_TO_MAT(CTSV(1:2,i)) * Q(CTSV(3,i):CTSV(3,i)+1,:);
    end
    V = V*Q;

% Nested functions
% -------------------------------------------------------------------------
function ShiftRotLeftLToRightL(i)
% 	% Shifts a rotation from left to right through the upper LR triangular
	LR(i:i+1,i:m) = CT_TO_MAT(Grot)*LR(i:i+1,i:m);
	[cos,sin,~]=CT_GIV(LR(i+1,i+1),LR(i+1,i));
	Grot = [cos,sin];
	LR(1:i+1,i:i+1) = LR(1:i+1,i:i+1)*CT_TO_MAT(Grot);
end

function ShiftRotLeftKToRightK(i)
	% Shifts a rotation from left to right through the upper LR triangular
	KR(i:i+1,i:m) = CT_TO_MAT(Grot)*KR(i:i+1,i:m);
	[cos,sin,~]=CT_GIV(KR(i+1,i+1),KR(i+1,i));
	Grot = [cos,sin];
	KR(1:i+1,i:i+1) = KR(1:i+1,i:i+1)*CT_TO_MAT(Grot);
end

function ShiftRotRightKToLeftK(i)
	% Shifts a rotation from right to left through the upper triangular
	KR(1:i+1,i:i+1) = KR(1:i+1,i:i+1)*CT_TO_MAT(Grot);
	[cos,sin,~]=CT_GIV(KR(i,i),KR(i+1,i));
	Grot = [cos,sin];
	KR(i:i+1,i:m) = CT_TO_MAT(Grot)*KR(i:i+1,i:m);
end

function ShiftRotRightLToLeftL(i)
	% Shifts a rotation from right to left through the upper triangular
	LR(1:i+1,i:i+1) = LR(1:i+1,i:i+1)*CT_TO_MAT(Grot);
	[cos,sin,~]=CT_GIV(LR(i,i),LR(i+1,i));
	Grot = [cos,sin];
	LR(i:i+1,i:m) = CT_TO_MAT(Grot)*LR(i:i+1,i:m);
end

end

