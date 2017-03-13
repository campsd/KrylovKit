%% (*) TEST CTRKS
n = 50;
global A LfacA UfacA pfacA B LfacB UfacB pfacB;
%Xi = {0,[0.12 0.17 0.09 0.42],0,0,Inf,Inf,[2 3 3 2], Inf,[0.1 0.5 0.6 0.3],Inf};
%Xi ={[2 3 5 6]};
%Xi = {Inf,Inf,Inf,Inf,Inf,[0.12 0.17 0.09 0.42],Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf};
Xi =
A = randn(n,n); B = eye(n,n);
[LfacA,UfacA,pfacA] = lu(A,'vector');
[LfacB,UfacB,pfacB] = lu(B,'vector');
T = [1 1 1; 0 0 1; 0 0 1; 0 0 0];
T = ones(11,10); T = triu(T);
T = eye(21,20);
v = randn(n,1);
eigAB = eig(A,B);
%% RKS variables
V = zeros(n,1); V2 = zeros(n,1);
V(:,1) = v/norm(v,2); V2(:,1) = v/norm(v,2);
Lrot = zeros(2,0);
KR = zeros(1,0); LR = zeros(1,0);
K2 = zeros(1,0); L2 = zeros(1,0);
%% CTRKS
[ V, KR, Lrot, LR ] = CTRKS( A, B, @funcpos_generalized, @funcneg_generalized, V, KR, Lrot, LR, Xi,T);
%
K = KR;
L = LR;
for i=length(Lrot):-1:1
    L(i:i+1,:) = CreateRotMat(Lrot(:,i)) * L(i:i+1,:);
end
norm(A*V*K-B*V*L,'fro')/norm(A*V*K,'fro')
A*V*K-B*V*L;
rval = eig(K(1:end-1,:),L(1:end-1,:));
%% CTIR restart
mu = -4.984244035477887;
[ Vr, KRr, Lrotr, LRr ] = CTIR( V, KR, Lrot, LR, mu );
%
Kr = KRr;
Lr = LRr;
for i=length(Lrotr):-1:1
    Lr(i:i+1,:) = CreateRotMat(Lrotr(:,i)) * Lr(i:i+1,:);
end
norm(A*Vr*Kr-B*Vr*Lr,'fro')/norm(A*Vr*Kr,'fro')
rvalrestart = eig(Kr(1:end-1,:),Lr(1:end-1,:));
%% RKS reference
[ V2, K2, L2 ] = RKS( A, B, @funcpos_generalized, @funcneg_generalized, V2, K2, L2, Xi,T);
norm(A*V2*K2-B*V2*L2,'fro')/norm(A*V2*K2,'fro')
rval2 = eig(K2(1:end-1,:),L2(1:end-1,:));
%% EXQZIR restart
[V2r, K2r, L2r] = EXQZIR(V2,K2,L2,mu);
norm(A*V2r*K2r-B*V2r*L2r,'fro')/norm(A*V2r*K2r,'fro')
rvalrestart2 = eig(K2r(1:end-1,:),L2r(1:end-1,:));
%% Compare
norm(sort(rval)-sort(rval2))
norm(sort(rvalrestart) - sort(rvalrestart2))
%% Plot
figure;
plot(real(eigAB),imag(eigAB),'.');
hold on
plot(real(rval),imag(rval),'o');
plot(real(rval2),imag(rval2),'s');
plot(real(rvalrestart),imag(rvalrestart),'d');
plot(real(rvalrestart2),imag(rvalrestart2),'x');
%% (*) Shift through testing (obsolete)
C = eye(3);
C(2:3,:) = CreateRotMat(Lrot(:,i)) * C(2:3,:);
C(1:2,:) = CreateRotMat(Lrot(:,i-1)) * C(1:2,:);
C(2:3,:) = CreateRotMat(RotH(Krot)) * C(2:3,:)
D = eye(3);
[CT1,CT2,CT3] = RotST(RotH(Krot), Lrot(:,i-1), Lrot(:,i));
D(1:2,:) = CreateRotMat(CT3) * D(1:2,:);
D(2:3,:) = CreateRotMat(CT2) * D(2:3,:);
D(1:2,:) = CreateRotMat(CT1) * D(1:2,:)
%%
C = eye(3);
C(1:2,:) = CreateRotMat(Lrot(:,i)) * C(1:2,:);
C(2:3,:) = CreateRotMat(Lrot(:,i-1)) * C(2:3,:);
C(1:2,:) = CreateRotMat(RotH(Krot)) * C(1:2,:)
D = eye(3);
[CT1,CT2,CT3] = RotST(RotH(Krot), Lrot(:,i-1), Lrot(:,i));
D(2:3,:) = CreateRotMat(CT1) * D(2:3,:);
D(1:2,:) = CreateRotMat(CT2) * D(1:2,:);
D(2:3,:) = CreateRotMat(CT3) * D(2:3,:)
%% (*) Comparison between CTRKS, RKS and Arnoldi
% We will study if the CTRKS relation is in agreement with a certain
% Arnoldi recursion and if so, which one. Also the agreement between CTRKS
% and RKS is investigated.
global A LfacA UfacA pfacA B LfacB UfacB pfacB;
n = 50;
%A = diag(1:1:n);
%A = A + 1i*A;
A = randn(n,n);
%B = eye(size(A));
%B = randn(n,n);
B = eye(n);
%[U,S,V] = svd(B); S(end,end) = 0; B = U*S*V'; %make B singular
[LfacA,UfacA,pfacA] = lu(A,'vector');
[LfacB,UfacB,pfacB] = lu(B,'vector');
v = ones(size(A,1),1);
%v = v + 1i*v;
%v = randn(size(A,1),1);
v=v/norm(v);
%Xi = {[1 -5.1 1 1],[1 -3.1+1i 1 1],[5 -1.1 1 1],[1 1.1 1 1],[1 3.1 1 1],[1 -50.1 1 1],[1 -170.1 1 1],[1 -190.1 1 1],[1 -195.1 1 1],[1 59.1 1 1],[1 59.1 1 1],[1 59.1 1 1],[1 59.1 1 1]};
Xi = {[1 -1 2 1],[1 -1.5 3 1],[1 -1.25 4 1],[1 -2.5 5 1],[1 -3.5 6 1],[1 -4.5 5.5 1],[1 -4.5 4.5 1],[1 -4.5 3.5 1],[1 -4.5 2.23 1],[1 -4.5 3.42 1],[1 -4.5 1.1 1]};
%Xi = {[1 -1 1 1],[1 -1.5 1 1]};
T = eye(n+1,n);
%T = randn(n+1,n); T = triu(T);
%% Run CTRKS
V = zeros(size(A,1),1);
V(:,1) = v;
Lrot = zeros(2,0);
KR = zeros(1,0); LR = zeros(1,0);
[ V, KR, Lrot, LR ] = CTRKS( A, B, @funcpos_generalized, @funcneg_generalized, V, KR, Lrot, LR, Xi,T);
L = LR;
for i=size(Lrot,2):-1:1
    L(i:i+1,:) = CreateRotMat(Lrot(:,i)) * L(i:i+1,:);
end
norm(A*V*KR-B*V*L,'fro')/norm(A*V*KR,'fro')
H = L/KR(1:end-1,:);
norm(A*V(:,1:end-1)-B*V*H,'fro')/norm(A*V(:,1),'fro')
rval = eig(KR(1:end-1,:),L(1:end-1,:));
%% Run Arnoldi with modified start vector
[U,Ha,flag] = arnoldiRef(B\A,V(:,1),length(Xi));
%% Run RKS
Vr = zeros(size(A,1),1);
Vr(:,1) = v;
Kr = zeros(1,0); Lr = zeros(1,0);
[ Vr, Kr, Lr ] = RKS( A, B, @funcpos_generalized, @funcneg_generalized, Vr, Kr, Lr, Xi,T);
norm(A*Vr*Kr-B*Vr*Lr,'fro')/norm(A*Vr*Kr,'fro')
[Lh,Kh,Q,Z] = hess(Lr(1:end-1,:),Kr(1:end-1,:));
%% Confirm that the span of V and Vr is the same
nrmsVVr = zeros(size(V,2),1);
for i=1:size(V,2)
    nrmsVVr(i) = norm(Vr(:,i)-V*V'*Vr(:,i),2);
end
% Span U and Vr ( this is much worse since it is numerically unstable )
nrmsUVr = zeros(size(V,2),1);
for i=1:size(V,2)
    nrmsUVr(i) = norm(Vr(:,i)-U*U'*Vr(:,i),2);
end
% Span U and V
nrmsUV = zeros(size(V,2),1);
for i=1:size(V,2)
    nrmsUV(i) = norm(V(:,i)-U*U'*V(:,i),2);
end
figure;
semilogy(nrmsVVr);
hold on
semilogy(nrmsUVr);
semilogy(nrmsUV);
legend('V vs Vr','U vs Vr', 'U vs V');
%% Regular plot
m = ceil(sqrt(size(V,2)));
n = m;
figure;
subplot(m,n,1)
hold on
plot(abs(V(:,1)));
vs = v;
for i=1:length(Xi)
    pole = Xi{i};
    vs = ((pole(1)*A+pole(2)*B)\(B)) * vs;
end
%
vs = v;
for i=1:length(Xi)
    pole = Xi{i};
    vs = (I - (pole(2)/pole(1)) * (A\B)) * vs;
end
%
plot(abs(vs/norm(vs)));
plot(abs(U(:,1)));
for i = 2:size(V,2)
    subplot(m,n,i);
    plot(abs(V(:,i)));
    hold on;
    plot(abs(U(:,i)));
    plot(abs(Vr(:,i)));
end
%% Semilogy plot
figure;
subplot(m,n,1)
semilogy(abs(V(:,1)));
hold on
semilogy(abs(vs/norm(vs)));
semilogy(abs(U(:,1)));
for i = 2:size(V,2)
    subplot(m,n,i);
    semilogy(abs(V(:,i)));
    hold on;
    semilogy(abs(U(:,i)));
    %semilogy(abs(Vr(:,i)));
end

%%
figure;
semilogy(abs(real(V(:,1))));
hold on
vsn = vs/norm(vs);
semilogy(abs(real(vsn)));
figure;
plot(abs(V(:,1))./abs(vsn));