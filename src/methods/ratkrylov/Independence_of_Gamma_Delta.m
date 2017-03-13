n = 10;
A = randn(n,n); B = randn(n,n); v =  randn(n,1);
%A = randn(n,n) + 1i * randn(n,n); B = randn(n,n) + 1i * randn(n,n); v =  randn(n,1) + 1i * randn(n,1);
[U,S,V] = svd(B);  
S = zeros(n,n); S(1,1) =1; S(2,2) = 0.5;
B = U*S*V';
v = V(:,2);
tol = 1e-10;
%% First vector
a = 2;
b = 1;

c1 = 200; c2 = 0;
d1 = 2; d2 = 4;
v=v/norm(v,2);

% compute
v11 = (a*A + b*B)\((c1*A + d1*B)*v);
v12 = (a*A + b*B)\((c2*A + d2*B)*v);
v11./v12
%orthogonalize
v11o = v11 - v*v'*v11;
v12o = v12 - v*v'*v12;
if (norm(v11o) < tol || norm(v12o) < tol), warning('RatArn Halted'), return,end
v11o./v12o
%normalize
v11o = v11o/norm(v11o);
v12o = v12o/norm(v12o);
v11o./v12o
%orthogonalize v2
PVo = eye(n,n) - v*v'; % Projector
PVo = v*v';
PVo = v';
v13 = PVo * ((a*A + b*B)\(A*v));
v14 = PVo * ((a*A + b*B)\(B*v));
v13./v14
% matrices
M1 = PVo * ((a*A + b*B)\(A));
M2 = PVo * ((a*A + b*B)\(B));
Mi = PVo / (a*A + b*B);
(M1*v)./(M2*v)
(Mi*A*v)./(Mi*B*v)
% t =randn(n,1);
% (M1*t)./(M2*t)
%C = B\A;
%OP = (a*eye(n) + b*inv(C))\(a*C + b*eye(n));
%OP2 = (a*C + b*eye(n))/(a*eye(n) + b*inv(C));
%% Second vector
vc = v11o;

a = 3;
b = 2;

c1 = 3; c2 = 4.3;
d1 = 5.1; d2 = 2.1;

v21 = (a*A + b*B)\((c1*A + d1*B)*vc);
v22 = (a*A + b*B)\((c2*A + d2*B)*vc);
v21./v22
v21o = v21 - v*v'*v21 - v11o*v11o'*v21;
v22o = v22 - v*v'*v22 - v12o*v12o'*v22;
if (norm(v21o) < tol || norm(v22o) < tol), warning('RatArn Halted'), return,end
v21o./v22o
%normalize
v21o = v21o/norm(v21o);
v22o = v22o/norm(v22o);
v21o./v22o
%% Span of 3dim space
V1 = [v v11o v21o];
V2 = [v v12o v22o];
V1 - V2*V2'*V1
%% Projectors 3dim space
%orthogonalize v2
Vo = [v v11o];
PVo2 = eye(n,n) - Vo*Vo'; % Projector
%PVo2 = Vo*Vo';
%PVo2 = Vo';
v13 = PVo2 * ((a*A + b*B)\(A*vc));
v14 = PVo2 * ((a*A + b*B)\(B*vc));
v13./v14
% matrices
M1 = PVo2 * ((a*A + b*B)\(A));
M2 = PVo2 * ((a*A + b*B)\(B));
Mi = PVo2 / (a*A + b*B);
(M1*v11o)./(M2*v11o)
(Mi*A*vc)./(Mi*B*vc)
