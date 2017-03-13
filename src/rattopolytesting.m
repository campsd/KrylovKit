clear all; clc;
global A LfacA UfacA pfacA B LfacB UfacB pfacB;
%% 1 POLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 200;
A = randn(n,n); B = eye(n,n);
v = randn(n,1);
[LfacA,UfacA,pfacA] = lu(A,'vector');
[LfacB,UfacB,pfacB] = lu(B,'vector');
T = [1; 0];
Xi = {[1 2 4 4]}; pole = Xi{1};
Vr = zeros(n,1); Vr(:,1) = v/norm(v,2);
K = zeros(1,0); L = zeros(1,0);
[ Vr, K, L ] = RKS( A, B, @funcpos_generalized, @funcneg_generalized, Vr, K, L, Xi,T);
[Q,~,~] = arnoldiRef(A,v,1);
V =  (pole(1) * A + pole(2) * B) * Vr;
%V = randn(n,2);
R = Q' * V;
norm(V - Q*R,'fro') % This must be very small for the spaces to span the same space
%% 2 POLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 20;
A = randn(n,n); B = randn(n,n);
v = randn(n,1);
[LfacA,UfacA,pfacA] = lu(A,'vector');
[LfacB,UfacB,pfacB] = lu(B,'vector');
T = [1 1; 0 0; 0 0];
Xi = {[1 3 2 1], [4 1 1 2.2]}; pole1 = Xi{1}; pole2 = Xi{2};
Vr = zeros(n,1); Vr(:,1) = v/norm(v,2);
K = zeros(1,0); L = zeros(1,0);
[ Vr, K, L ] = RKS( A, B, @funcpos_generalized, @funcneg_generalized, Vr, K, L, Xi,T);
[Q,~,~] = arnoldiRef(B\A,v,2);
V =  (B\(pole2(1) * A + pole2(2) * B)) * (B\(pole1(1) * A + pole1(2) * B)) * Vr;
%V = randn(n,2);
R = Q' * V;
norm(V - Q*R,'fro') % This must be very small for the spaces to span the same space
%%
pole3 = [2.5 1 1 2];
[Qs,~] = qr(pole3(1)*L + pole3(2)*K);
q = Qs(:,end)
sigma = -pole3(2)/pole3(1);
s3 = [sigma^2;sigma; 1]
t3 = R*s3
Xi = {pole3};
T = [q; 0];
[ Vr, K, L ] = RKS( A, B, @funcpos_generalized, @funcneg_generalized, Vr, K, L, Xi,T);
%%
