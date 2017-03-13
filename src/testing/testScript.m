clear all; close all; clc;
format short e;
m = 1210; B = 5.45; C = 2; Du = 0.004; Dv = 0.008; L = 1;
A = brusselator(m,C,B,Dv,Du,L);
A = sparse(A);
P = amd(A); A2 = A(P,P);
global Lfac Ufac pfac;
%%
tic;
[Lfac,Ufac,pfac] = lu(A,'vector');
toc;
%%
tic;
[Lfac,Ufac,pfac] = lu(A2,'vector');
toc;
%%
v=ones(size(A,1),1); v=v/norm(v);
n = size(A,1);
V = zeros(n,1);
V(:,1) = v/norm(v,2);
KLrot = zeros(2,0); KLidx = zeros(1,0);
KR = zeros(1,0); LR = zeros(1,0);
global Lfac Ufac pfac;
[Lfac,Ufac,pfac] = lu(A,'vector');

pattern = [-1,-1,-1,1];
selection = repmat(pattern,[1,5]);
%new
[V,KLrot,KLidx,KR,LR] = CTEK(@funcpos,@funcneg,V,KLrot,KLidx,KR,LR,selection);
mu = 1;
%[V,KLrot,KLidx,KR,LR] = CTIF(A,V,KLrot,KLidx,KR,LR,mu,0);

[V,KLrot,KLidx,KR,LR] = CTIR(V,KLrot,KLidx,KR,LR,mu);
%[V,KLrot,KLidx,KR,LR] = CTEK(@funcpos,@funcneg,V,KLrot,KLidx,KR,LR,selection);
[K,L] = CONS_CTEK_PENCIL(KLrot,KLidx,KR,LR);
norm(A*V*K-V*L,'fro')

%% CTPC testing
PC = [];
[V,KLrot,KLidx,KR,LR] = CTEK(@funcpos,@funcneg,V,KLrot,KLidx,KR,LR,selection);
[ PC ] = CTPC( KLrot, KLidx, KR, LR, PC );
[V,KLrot,KLidx,KR,LR] = CTEK(@funcpos,@funcneg,V,KLrot,KLidx,KR,LR,selection);
[ PC ] = CTPC( KLrot, KLidx, KR, LR, PC );
% Check correctness
diff = PC-V(:,1:end-1)'*A*V(:,1:end-1);
norm(PC-V(:,1:end-1)'*A*V(:,1:end-1),'fro')/norm(PC,'fro')

figure, colormap('summer')
imagesc(log10(abs(diff)))
colorbar, set(gca,'CLim',[-15,0]); axis square

%% figure;plot(eig(PC),'o');
figure;
plot(eig(PC),'o');

%% old
[VO,KGO,KRO, LGO, LRO, res] = RIREK_REFERENCE(A, v, 5,pattern);

figure;
spy(K);
figure;
spy(L);