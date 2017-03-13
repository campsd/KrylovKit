function [V,KR,Lrot,LR] = CTRKS1(A,B,funcPole,poleFac,Xi,V,KR,Lrot,LR,T)
%CTRKS1 - Core-Transformations factorized rational Krylov sequence
% This algorithm executes one rational Krylov step

start_idx = size(V,2);

% zero padding arrays
V = padarray( V, [0 1], 0, 'post');
KR = padarray( KR, [1 1], 0, 'post');
LR = padarray( LR, [1 1], 0,'post');
Lrot = padarray(Lrot, [0 1], 'post');

curr_ip1 = start_idx + 1;
curr_i = curr_ip1 - 1;

% (1) Execute the operation
w = funcPole(1/poleFac(1),(Xi(3)*A + Xi(4)*B)*(V(:,1:curr_i)*T(1:curr_i)));
%w = (Xi(1)*A + Xi(2)*B) \ ((Xi(3)*A + Xi(4)*B)*(V(:,1:curr_i)*T(1:curr_i)));
% (2) Orthogonalise
h = zeros(curr_ip1,1);
for kk=1:2
    for j=1:curr_i
            hc = V(:,j)'*w;
            h(j) = h(j) + hc;
            w = w - hc*V(:,j);
    end
end
h(curr_ip1) = norm(w,2);

% (3) Update recursion 
V(:,curr_ip1) = w/h(curr_ip1);
ki = Xi(1)*h - Xi(3)*T;
li = -Xi(2)*h + Xi(4)*T;
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
[V,Lrot,KR,LR] = ChaseKrotUp(V,Lrot,KR,LR,Krot,curr_i);

end

function [li] = ApplyLrot(li,Lrot,curr_i)
    % Applies the rotations in Lrot to li
    for i=1:curr_i-1
        li(i:i+1) = CreateRotMat(RotH(Lrot(:,i))) * li(i:i+1);
    end
end

function [V,Lrot,KR,LR] = ChaseKrotUp(V,Lrot,KR,LR,Krot,curr_i)
    % Chases the Krot upwards and fuses it with the first rotation
    for i=curr_i:-1:2
       V(:,i:i+1) = V(:,i:i+1) * CreateRotMat(Krot);
       [Lrot(:,i-1), Lrot(:,i), Krot] = RotST(RotH(Krot), Lrot(:,i-1), Lrot(:,i));
       [Krot, LR] = ShiftRotLeftToRight(Krot,LR,i-1,curr_i);
       [Krot, KR] = ShiftRotRightToLeft(Krot,KR,i-1,curr_i);
       Krot = RotH(Krot);
    end
    V(:,1:2) = V(:,1:2) * CreateRotMat(Krot);
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
