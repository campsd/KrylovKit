function [ V, KR, Lrot, LR ] = CTRKTOSK( A, B, funpos, funneg, V, KR, Lrot, LR, Xi,T)
%CTRKS - Core-Transformations factorized rational Krylov sequence
%   This algorithm creates and updates a rational Krylov recursion with the
%   recurrence matrices in QR factorised format. The core transformations
%   are manipulated such that they all appear at the L side. This
%   essentially shifts all the poles to infinity and transforms the RK
%   relation to an Arnoldi relation.
%
%   INPUT
%   A,B     matrices of the generalised eigenvalue problem
%   funpos  function handle that executes the matvec product AB^-1
%   funneg  function handle that solves the system BA^-1
%   V       matrix with Krylov basis vectors
%   Krot    core transformation acting on last two rows of K
%   KR      upper triangular matrix for K
%   Lrot    core transformations of L
%   LR      upper triangular matrix for L
%   Xi      cell array with poles. Options: 0->funneg, Inf->funpos, or
%           [alpha, beta, gamma, delta]
%   T       matrix with continuation vectors
%
%   OUTPUT
%   V
%   Krot
%   KR
%   Lrot
%   LR

    k = length(Xi);
    start_idx = size(V,2);
    
    % zero padding arrays to appropriate size
    V = padarray( V , [0 k], 0, 'post');
    KR = padarray( KR , [k k], 0, 'post');
    LR = padarray( LR , [k k], 0, 'post');
    Lrot = padarray( Lrot, [0 k], -1, 'post');
    
    % main loop
    for i=1:k
        curr_ip1 = start_idx + i;
        curr_i = curr_ip1 - 1;
        
        % (1) Execute the operation
        if Xi{i} == Inf %AB^-1v
            w = funpos(V(:,1:curr_i)*T(1:curr_i,i));
        elseif Xi{i} == 0 %BA^-1v
            w = funneg(V(:,1:curr_i)*T(1:curr_i,i));
        else %[alpha, beta, gamma, delta]
            pole = Xi{i};
            w = (pole(1)*A + pole(2)*B) \ ((pole(3)*A + pole(4)*B)*(V(:,1:curr_i)*T(1:curr_i,i)));
        end
        
        % (2) Orthogonalise the new vector
        h = zeros(curr_ip1,1);
        for kk=1:2
            for j=1:curr_i
                    hc = V(:,j)'*w;
                    h(j) = h(j) + hc;
                    w = w - hc*V(:,j);
            end
        end
        h(curr_ip1) = norm(w,2);
        
        % (3) Update the recursion
        V(:,curr_ip1) = w/h(curr_ip1);
        if Xi{i} == Inf %AB-1v
            % In this case there only appears a rotation at the L side 
            % Vectors
            ki = -T(1:curr_ip1,i);
            li = -h;
            % Apply the previous rotations 
            [li] = ApplyLrot(li,Lrot,curr_i);
            % Compute new rotation at the L side
            [c,s,r] = RotGIV(li(curr_i),li(curr_ip1));
            li(curr_i:curr_ip1) = [r,0];
            % Update upper triangulars
            KR(1:curr_ip1,curr_i) = ki;
            LR(1:curr_ip1,curr_i) = li;
            % New rotations
            Lrot(:,curr_i) = [conj(c); -s]; 
        elseif Xi{i} == 0 %BA-1v
            % In this case there only appears a rotation at the K side
            % Vectors
            ki = h;
            li = T(1:curr_ip1,i);
            % Apply the previous rotations 
            [li] = ApplyLrot(li,Lrot,curr_i);
            % Compute the new rotation at the K side
            [c,s,r] = RotGIV(ki(curr_i),ki(curr_ip1));
            ki(curr_i:curr_ip1) = [r,0];
            % Update upper triangulars
            KR(1:curr_ip1,curr_i) = ki;
            LR(1:curr_ip1,curr_i) = li;
            % New rotations
            Krot = [conj(c); -s];
            Lrot(:,curr_i) = [1; 0];
            % Chase Krot up
            [V,Lrot,KR,LR,T] = ChaseKrotUp(V,Lrot,KR,LR,Krot,curr_i,T);
        else %[alpha, beta, gamma, delta]
            % In this case there appear rotations at both sides
            % Vectors
            ki = pole(1)*h - pole(3)*T(1:curr_ip1,i);
            li = -pole(2)*h + pole(4)*T(1:curr_ip1,i);
            % Apply the previous rotations 
            [li] = ApplyLrot(li,Lrot,curr_i);
            % Compute the new rotations at both sides
            [ck,sk,rk] = RotGIV(ki(curr_i),ki(curr_ip1));
            ki(curr_i:curr_ip1) = [rk,0];
            [cl,sl,rl] = RotGIV(li(curr_i),li(curr_ip1));
            li(curr_i:curr_ip1) = [rl,0];
            % Update upper triangulars
            KR(1:curr_ip1,curr_i) = ki;
            LR(1:curr_ip1,curr_i) = li;
            % New rotations
            Lrot(:,curr_i) = [conj(cl); -sl]; 
            Krot = [conj(ck); -sk];
            % Chase the previous Krot up
            [V,Lrot,KR,LR,T] = ChaseKrotUp(V,Lrot,KR,LR,Krot,curr_i,T);
        end
    end
end

function [li] = ApplyLrot(li,Lrot,curr_i)
    % Applies the rotations in Lrot to li
    for i=1:curr_i-1
        li(i:i+1) = CreateRotMat(RotH(Lrot(:,i))) * li(i:i+1);
    end
end

function [V,Lrot,KR,LR,T] = ChaseKrotUp(V,Lrot,KR,LR,Krot,curr_i,T)
    % Chases the Krot upwards and fuses it with the first rotation
    for i=curr_i:-1:2
       V(:,i:i+1) = V(:,i:i+1) * CreateRotMat(Krot);
       %T(i:i+1,:) = CreateRotMat(RotH(Krot)) * T(i:i+1,:);
       [Lrot(:,i-1), Lrot(:,i), Krot] = RotST(RotH(Krot), Lrot(:,i-1), Lrot(:,i));
       [Krot, LR] = ShiftRotLeftToRight(Krot,LR,i-1,curr_i);
       [Krot, KR] = ShiftRotRightToLeft(Krot,KR,i-1,curr_i);
       Krot = RotH(Krot);
    end
    V(:,1:2) = V(:,1:2) * CreateRotMat(Krot);
    %T(1:2,:) = CreateRotMat(RotH(Krot)) * T(1:2,:); % Issues in complex
    %arithmetic
    Lrot(:,1) = RotFUS(RotH(Krot),Lrot(:,1));
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