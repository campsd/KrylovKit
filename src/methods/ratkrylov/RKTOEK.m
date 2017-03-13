function [ V, KLrot, KLidx, KR, LR ] = RKTOEK( V, Krot, KR, Lrot, LR, s, stop )
%RKTOEK - rational krylov to extended krylov transformation
%   This function transforms a rational Krylov pencil (V,K,L) to an
%   extended Krylov pencil (V,K,L) by shifting the poles either to infinity
%   (s_i = 1) or zero (s_i = -1).
%
%   INPUT
%   V       rational Krylov basis matrix
%   Krot    CTs for K
%   KR      upper triangular for K
%   Lrot    CTs for L
%   LR      upper triangular for L
%   s       selection vector of the final EK subspace
%   stop    boolean. If true, the last CTs are left as they are (Ritz
%           values remain unaffected). If false, the transformation is
%           carried out until the end
%
%   OUTPUT
%   V       extended Krylov basis matrix with same span as input V
%   KLrot   condensed CTS for extended Krylov recursion
%   KLidx   indicates the side of the EK CTs (0 = K, 1 = L)
%   KR      upper triangular for K
%   LR      upper triangular for L
%
%   NOTES
%   
%
%   daan.camps@cs.kuleuven.be
%   March 10, 2017

% Check if dimensions match
assert(size(KR) == size(LR),'Dimension mismatch 1');
assert(size(V,2) == size(KR,1),'Dimension mismatch 2');
assert(size(Lrot) == size(Krot),'Dimension mismatch 3');
assert(size(Lrot,2) == size(KR,2),'Dimension mismatch 4');
assert(length(s) == size(Lrot,2),'Dimension mismatch 5');

% Initialize some variables
n = size(Lrot,2);   % Hence Hessenberg matrices are n+1xn, V is Nxn+1 and we have n CTs 
CTSV = zeros(3,0);  % Collects the CTs that need to be applied to V in order to sustian the recursion
PyRot = zeros(3,0); % This will contain the "pyramid" of CTs that will be condensed

% To start we bring all CTs from K to L
for i=1:n
    PyRot(1:2,i) = RotH(Krot(:,n-i+1));
    PyRot(3,i) = n-i+1;
    CTSV(1:2,i) = Krot(:,i);
    CTSV(3,i) = i;
end

% Fuse the top one
PyRot(1:2,n) = RotFUS(PyRot(:,n), Lrot(:,1));
% Add remaining L rotations
PyRot(1:2,n+1:2*n-1) = Lrot(:,2:n);
PyRot(3,n+1:2*n-1) = 2:n;

dirswitch = 1; %keeps the idx of the row with last change in direction
% Main loop
for i=2:n
    curr_i = i;
    % Pick the CT that needs to be removed
    if s(i-1) > 0 % the leftmost needs to be removed
       if (i>2) &&  s(i-2) < 0, dirswitch = i-1; end
       col = find(PyRot(3,:)==i,1,'first');
       ChaseRot = PyRot(1:2,col);
       while curr_i > dirswitch
           [PyRot(1:2,curr_i+col-1),PyRot(1:2,curr_i+col), ChaseRot] = RotST(ChaseRot,PyRot(1:2,curr_i+col-1),PyRot(1:2,curr_i+col));
           curr_i = curr_i -1;
           % update i and col in while loop
           ShiftRotLeftToRight(curr_i);
           ShiftRotRightToLeft(curr_i);
           CTSV(1:2,end+1) = ChaseRot;
           CTSV(3,end) = curr_i;
           ChaseRot = RotH(ChaseRot); % !!
       end
       % Fuse CTS
       PyRot(1:2,curr_i+col-1) = RotFUS(ChaseRot, PyRot(1:2,curr_i+col-1));
       % delete the initial ChaseRot from the PyRot pattern
       PyRot(:,col) = [];
    else % the rightmost needs to be removed
       if (i>2) &&  s(i-2) < 0, dirswitch = i-1; end
       col = find(PyRot(3,:)==i,1,'last');
       ChaseRot = PyRot(1:2,col);
       while curr_i > dirswitch
           [ChaseRot, PyRot(1:2,col-curr_i-1),PyRot(1:2,col-curr_i)] = RotST(PyRot(1:2,col-curr_i-1),PyRot(1:2,col-curr_i),ChaseRot);
           curr_i = curr_i - 1;
           CTSV(1:2,end+1) = ChaseRot;
           CTSV(3,end) = curr_i;
           ChaseRot = RotH(ChaseRot)
           ShiftRotRightToLeft(curr_i);
           ShiftRotLeftToRight(curr_i);
           ChaseRot = RotH(ChaseRot);
       end
       % Fuse CTS
       PyRot(1:2,col-curr_i-1) = RotFUS(PyRot(1:2,col-curr_i-1), ChaseRot);
       % delete the initial ChaseRot from the PyRot pattern
       PyRot(:,col) = [];
    end
    
end

function ShiftRotLeftToRight(i)
	% Shifts a rotation from left to right through the upper LR triangular
	LR(i:i+1,i:k) = CreateRotMat(ChaseRot)*LR(i:i+1,i:n);
	[cos,sin,~]=RotGIV(LR(i+1,i+1),LR(i+1,i));
	ChaseRot = [cos,sin];
	LR(1:i+1,i:i+1) = LR(1:i+1,i:i+1)*CreateRotMat(ChaseRot);
end

function [Grot, R] = ShiftRotRightToLeft(Grot, R, i, k)
	% Shifts a rotation from right to left through the upper triangular
	KR(1:i+1,i:i+1) = KR(1:i+1,i:i+1)*CreateRotMat(ChaseRot);
	[cos,sin,~]=RotGIV(KR(i,i),KR(i+1,i));
	ChaseRot = [cos,sin];
	KR(i:i+1,i:n) = CreateRotMat(ChaseRot)*KR(i:i+1,i:k);
end

end
