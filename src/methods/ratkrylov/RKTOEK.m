function [ V, KLrot, KLidx, KR, LR ] = RKTOEK( A,B,V, Krot, KR, Lrot, LR, s )
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
%   This version first moves all rotations to one side, modifies the
%   pattern and shifts them back. The last rotation is always at the L side
%   in this version.
%   
%
%   daan.camps@cs.kuleuven.be
%   March 10, 2017

% Check if dimensions match
assert(sum(size(KR) == size(LR))==2,'Dimension mismatch 1');
assert(size(V,2) == size(KR,1),'Dimension mismatch 2');
assert(sum(size(Lrot) == size(Krot))==2,'Dimension mismatch 3');
assert(size(Lrot,2) == size(KR,2),'Dimension mismatch 4');
assert(length(s) == size(Lrot,2),'Dimension mismatch 5');

% Initialize some variables
n = size(Lrot,2);   % Hence Hessenberg matrices are n+1xn, V is Nxn+1 and we have n CTs 
CTSV = zeros(3,0);  % Collects the CTs that need to be applied to V in order to sustian the recursion
PyRot = zeros(2,0); % This will contain the "pyramid" of CTs that will be condensed
PyInd = zeros(2,0); % This contains the row indices and order of the rotations
% To start we bring all CTs from K to L
for i=1:n
    PyRot(1:2,i) = RotH(Krot(:,n-i+1));
    PyInd(1,i) = n-i+1; %row idx
    PyInd(2,i) = 2*n - i+1; %order
    CTSV(1:2,i) = Krot(:,i);
    CTSV(3,i) = i;
end

% Fuse the top one
PyRot(1:2,n) = RotFUS(PyRot(:,n), Lrot(:,1));
% Add remaining L rotations
PyRot(1:2,n+1:2*n-1) = Lrot(:,2:n);
PyInd(1,n+1:2*n-1) = 2:n; % row
PyInd(2,n+1:2*n-1) = n:-1:2; % order

% At this moment we have all rotations at the L side in an empty
% pyramid shape: /\

% Test for correctness:
% -------------------------------------------------------------------------
% Apply the rotations to V
Q = eye(size(KR,1));
for i=size(CTSV,2):-1:1
    Q(CTSV(3,i):CTSV(3,i)+1,:) = CreateRotMat(CTSV(1:2,i)) * Q(CTSV(3,i):CTSV(3,i)+1,:);
end
Vt = V*Q;

L = LR;
[~,I] = sort(PyInd(2,:));
for i=1:size(PyRot,2)
    L(PyInd(1,I(i)):PyInd(1,I(i))+1,:) = CreateRotMat(PyRot(1:2,I(i)))*L(PyInd(1,I(i)):PyInd(1,I(i))+1,:);
end
norm(A*Vt*KR-B*Vt*L,'fro')/norm(A*Vt*KR,'fro')
% -------------------------------------------------------------------------

% Now we'll modify the py ramid shape to a condensed pattern satisfing the
% selection vector s. The pattern is built up from bottom to top: while
% chasing core transformations up, the pattern gradually moves downwards.

KLrot = zeros(2,n); % Will contain final rotation pattern
KLidx = zeros(2,n);
KLrot(:,1) = PyRot(:,n);
KLidx(:,1) = [1; n];
PyRot(:,n) = [];
PyInd(:,n) = [];  
% Main loop
for i=2:n
    % We do the initial turnover decide which one the ChaseRot will be and
    % chase it upwards. At the top we select on of both based on the entry
    % of s.
    col_i = find(PyInd(1,:)==i); % two results
        
    if ((i==2) && (s(n-1) > 0)) || ((i>2) && (KLidx(2,i-2) > KLidx(2,i-1))) % The i-2 rotation comes after the i-1 (right will be free after turnover)
        % Hessenberg
        %[KLrot(:,i-1),KLrot(:,i),ChaseRot] = RotST( PyRot(:,col_i(1)),KLrot(:,i-1),PyRot(:,col_i(2)));
        ChaseRot = PyRot(:,col_i(1));
        KLrot(:,i) = PyRot(:,col_i(2));
        KLidx(:,i) = [i; KLidx(2,i-1)-1];
    else % Left will be free after turnover
        % inv Hessenberg
        %[ChaseRot,KLrot(:,i),KLrot(:,i-1)] = RotST( PyRot(:,col_i(1)),KLrot(:,i-1),PyRot(:,col_i(2)));
        ChaseRot = PyRot(:,col_i(2));
        KLrot(:,i) = PyRot(:,col_i(1));
        KLidx(:,i) = [i; KLidx(2,i-1)+1];
    end
    % remove the rotation at row i from PyRot
    PyRot(:,col_i) = [];
    PyInd(:,col_i) = [];    
    
    %ChaseRot is currently situated at row i-1, need to chase it up to row 1
    for k=i:-1:3        
        if KLidx(2,k) < KLidx(2,k-1) % Hessenberg shape
            if KLidx(2,k-1) > KLidx(2,k-2) % inverse Hessenberg shape (DIR!)
               % Turnover
               [ChaseRot,KLrot(:,k),KLrot(:,k-1)] = RotST( ChaseRot,KLrot(:,k-1),KLrot(:,k));
               % The order of rotations has changed
               %KLidx(2,k) = KLidx(2,k) + 2;
               KLidx(2,1:k-1) = KLidx(2,1:k-1) - 2;
               % Left removal
               CTSV(1:2,end+1) = ChaseRot;
               CTSV(3,end) = k-1;
               ChaseRot = RotH(ChaseRot);
               % Double transfer
               ShiftRotLeftKToRightK(k-1);
               ShiftRotRightLToLeftL(k-1);
               ChaseRot = RotH(ChaseRot);
            else % also Hessenberg shape
               % Turnover
               [KLrot(:,k-1),KLrot(:,k), ChaseRot] = RotST( ChaseRot,KLrot(:,k-1),KLrot(:,k));
               % Double transfer
               ShiftRotLeftLToRightL(k-1);
               ShiftRotRightKToLeftK(k-1);
               % Left removal
               CTSV(1:2,end+1) = RotH(ChaseRot);
               CTSV(3,end) = k-1;
            end
        else % inverse Hessenberg shape
            if (KLidx(2,k-1) < KLidx(2,k-2)) % Hessenberg shape (DIR!)
               % Turnover
               [KLrot(:,k-1),KLrot(:,k),ChaseRot] = RotST(KLrot(:,k),KLrot(:,k-1),ChaseRot);
               % The order of rotations has changed
               KLidx(2,k:i) = KLidx(2,k:i) - 2;
               % Double transfer
               ShiftRotLeftLToRightL(k-1);
               ShiftRotRightKToLeftK(k-1);
               % Left removal
               CTSV(1:2,end+1) = RotH(ChaseRot);
               CTSV(3,end) = k-1;
            else % also inverse Hessenberg shape
               % Turnover
               [ChaseRot, KLrot(:,k),KLrot(:,k-1)] = RotST(KLrot(:,k),KLrot(:,k-1),ChaseRot);
               % Left removal
               CTSV(1:2,end+1) = ChaseRot;
               CTSV(3,end) = k-1;
               ChaseRot = RotH(ChaseRot);
               ShiftRotLeftKToRightK(k-1);
               ShiftRotRightLToLeftL(k-1);
               ChaseRot = RotH(ChaseRot);
            end
        end
    end
    
    % Final operation - the ChaseRot has now been chased to the second row of
    % KLrot. We perform the final turnover. Based on the entry s(n-i+1), 
    % we keep one of both rotations at the first row.
    
    if KLidx(2,2) < KLidx(2,1) % Hessenberg shape
        % This means that the ChaseRot is at the left side on row 2
        if s(n-i+1) >0 % we want a Hessenberg order
           [KLrot(:,1),KLrot(:,2), ChaseRot] = RotST( ChaseRot,KLrot(:,1),KLrot(:,2));
            % Double transfer
           ShiftRotLeftLToRightL(1);
           ShiftRotRightKToLeftK(1);
           % Left removal
           CTSV(1:2,end+1) = RotH(ChaseRot);
           CTSV(3,end) = 1;
           % Fuse
           KLrot(:,1) = RotFUS(ChaseRot, KLrot(:,1));   
        else % we want an inverse Hessenberg order
           % Turnover
           [ChaseRot,KLrot(:,2),KLrot(:,1)] = RotST( ChaseRot,KLrot(:,1),KLrot(:,2));
           % The order of rotations has changed
           %KLidx(2,k) = KLidx(2,k) + 2;
           KLidx(2,1) = KLidx(2,1) - 2;
           % Left removal
           CTSV(1:2,end+1) = ChaseRot;
           CTSV(3,end) = 1;
           ChaseRot = RotH(ChaseRot);
           % Double transfer
           ShiftRotLeftKToRightK(1);
           ShiftRotRightLToLeftL(1);
           ChaseRot = RotH(ChaseRot);
           % Fuse
           KLrot(:,1) = RotFUS(KLrot(:,1), ChaseRot);
        end
    else % inv Hessenberg shape
        % This means that the ChaseRot is at the right side on row 2
        if s(n-i+1) > 0 % we want a Hessenberg order
            % Turnover
           [KLrot(:,1),KLrot(:,2),ChaseRot] = RotST(KLrot(:,2),KLrot(:,1),ChaseRot); 
           % The order of rotations has changed
           KLidx(2,2:i) = KLidx(2,2:i) - 2;
           % Double transfer
           ShiftRotLeftLToRightL(1);
           ShiftRotRightKToLeftK(1);
           % Left removal
           CTSV(1:2,end+1) = RotH(ChaseRot);
           CTSV(3,end) = 1;
           % Fuse
           KLrot(:,1) = RotFUS(ChaseRot, KLrot(:,1));
        else % we want an inverse Hessenberg order
            % Turnover
           [ChaseRot, KLrot(:,2),KLrot(:,1)] = RotST(KLrot(:,2),KLrot(:,1),ChaseRot);
           % Left removal
           CTSV(1:2,end+1) = ChaseRot;
           CTSV(3,end) = 1;
           ChaseRot = RotH(ChaseRot);
           ShiftRotLeftKToRightK(1);
           ShiftRotRightLToLeftL(1);
           ChaseRot = RotH(ChaseRot);
           % Fuse
           KLrot(:,1) = RotFUS(KLrot(:,1), ChaseRot);
        end
        
    end

end
           
% Test for correctness:
% -------------------------------------------------------------------------
% Apply the rotations to V
Q = eye(size(KR,1));
for i=size(CTSV,2):-1:1
    Q(CTSV(3,i):CTSV(3,i)+1,:) = CreateRotMat(CTSV(1:2,i)) * Q(CTSV(3,i):CTSV(3,i)+1,:);
end
V = V*Q;

L = LR;
[~,I] = sort(KLidx(2,:));
for i=1:size(KLrot,2)
    L(KLidx(1,I(i)):KLidx(1,I(i))+1,:) = CreateRotMat(KLrot(1:2,I(i)))*L(KLidx(1,I(i)):KLidx(1,I(i))+1,:);
end
norm(A*V*KR-B*V*L,'fro')/norm(A*V*KR,'fro')
% -------------------------------------------------------------------------

% Now we need to bring the correct rotations to the K side
KLidx = zeros(n,1);
for i=1:n-1
    if s(i) > 0 % CT at correct side
        %KLrot(:,i) = KLrot(:,i);
        KLidx(i) = 1;
    else % Transfer to the other side
        ChaseRot = KLrot(:,i);
        ShiftRotLeftLToRightL(i);
        ShiftRotRightKToLeftK(i);
        KLrot(:,i) = RotH(ChaseRot);
        KLidx(i) = 0;
    end
end
KLidx(n) = 1;

% Test again
[K,L] = CONS_CTEK_PENCIL(KLrot,KLidx,KR,LR);
norm(A*V*K-B*V*L,'fro')/norm(A*V*K,'fro')

function ShiftRotLeftLToRightL(i)
	% Shifts a rotation from left to right through the upper LR triangular
	LR(i:i+1,i:n) = CreateRotMat(ChaseRot)*LR(i:i+1,i:n);
	[cos,sin,~]=RotGIV(LR(i+1,i+1),LR(i+1,i));
	ChaseRot = [cos,sin];
	LR(1:i+1,i:i+1) = LR(1:i+1,i:i+1)*CreateRotMat(ChaseRot);
end

function ShiftRotLeftKToRightK(i)
	% Shifts a rotation from left to right through the upper LR triangular
	KR(i:i+1,i:n) = CreateRotMat(ChaseRot)*KR(i:i+1,i:n);
	[cos,sin,~]=RotGIV(KR(i+1,i+1),KR(i+1,i));
	ChaseRot = [cos,sin];
	KR(1:i+1,i:i+1) = KR(1:i+1,i:i+1)*CreateRotMat(ChaseRot);
end

function ShiftRotRightKToLeftK(i)
	% Shifts a rotation from right to left through the upper triangular
	KR(1:i+1,i:i+1) = KR(1:i+1,i:i+1)*CreateRotMat(ChaseRot);
	[cos,sin,~]=RotGIV(KR(i,i),KR(i+1,i));
	ChaseRot = [cos,sin];
	KR(i:i+1,i:n) = CreateRotMat(ChaseRot)*KR(i:i+1,i:n);
end

function ShiftRotRightLToLeftL(i)
	% Shifts a rotation from right to left through the upper triangular
	LR(1:i+1,i:i+1) = LR(1:i+1,i:i+1)*CreateRotMat(ChaseRot);
	[cos,sin,~]=RotGIV(LR(i,i),LR(i+1,i));
	ChaseRot = [cos,sin];
	LR(i:i+1,i:n) = CreateRotMat(ChaseRot)*LR(i:i+1,i:n);
end

end
