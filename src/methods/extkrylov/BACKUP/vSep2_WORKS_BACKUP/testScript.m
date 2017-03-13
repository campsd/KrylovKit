clear all; close all; clc;
format short e;
m = 484; B = 5.45; C = 2; Du = 0.004; Dv = 0.008; L = 1;
A = brusselator(m,C,B,Dv,Du,L);
A = sprand(1500,1500,0.01);
%v=ones(size(A,1),1); v=v/norm(v);
v = randn(size(A,1),1); v=v/norm(v);
A\v;
n = size(A,1);
V = zeros(n,1);
V(:,1) = v/norm(v,2);
KLrot = zeros(2,0); KLidx = zeros(1,0);
K = zeros(1,0); L = zeros(1,0);      
KR = zeros(1,0); LR = zeros(1,0);
[Lfac,Ufac,pfac] = lu(A,'vector');

pattern = [+1,-1,-1,-1,+1];
selection = repmat(pattern,[1,2]);
%% new
[V,KLrot,KLidx,KR,LR] = CTEK(@funcpos,@funcneg,V,KLrot,KLidx,KR,LR,selection);

[K,L] = CONS_CTEK_PENCIL(KLrot,KLidx,KR,LR);
norm(A*V*K-V*L,'fro')
%%
norm(L(1:end-1,:)/K(1:end-1,:) - V(:,1:end-1)'*A*V(:,1:end-1),'fro')

rank(K(1:end-1,:))
%%
[V,K,L] = EK(@funcpos,@funcneg,V,K,L,selection);
%% residual in case of s(k) = -1
I = eye(size(K,2),size(K,2));
RES = K(end,end) * V(:,1:end-1)'*A*V(:,end)*I(end,:)/K(1:end-1,:);
norm(RES,'fro')
norm(L(1:end-1,:)/K(1:end-1,:) - RES - V(:,1:end-1)'*A*V(:,1:end-1),'fro')

figure, colormap('summer')
imagesc(log10(abs(RES)))
colorbar, set(gca,'CLim',[-15,0]); axis square

figure, colormap('summer')
imagesc(log10(abs(L(1:end-1,:)/K(1:end-1,:) + RES)))
colorbar, set(gca,'CLim',[-15,0]); axis square

figure, colormap('summer')
imagesc(log10(abs(V(:,1:end-1)'*A*V(:,1:end-1))))
colorbar, set(gca,'CLim',[-15,0]); axis square

figure, colormap('summer')
imagesc(log10(abs(L(1:end-1,:)/K(1:end-1,:) + RES - V(:,1:end-1)'*A*V(:,1:end-1))))
colorbar, set(gca,'CLim',[-15,0]); axis square
%%
mu = ones(10,1);
[V,KLrot,KLidx,KR,LR] = CTIR(V,KLrot,KLidx,KR,LR,mu);
[V,KLrot,KLidx,KR,LR] = CTEK(@funcpos,@funcneg,V,KLrot,KLidx,KR,LR,selection);
[K,L] = CONS_CTEK_PENCIL(KLrot,KLidx,KR,LR);
norm(A*V*K-V*L,'fro')
%% old
[VO,KGO,KRO, LGO, LRO, res] = RIREK_REFERENCE(A, v, 5,pattern);
%%
tol = 1e-15;
for i=1:size(K,2)
    if abs(K(i,i)) < tol
        K(i,i) = 0;
    end
    if abs(L(i,i)) < tol
        L(i,i) = 0;
    end
end
%%
figure, colormap('summer')
imagesc(log10(abs(K)))
colorbar, set(gca,'CLim',[-15,0]); axis square
figure, colormap('summer')
imagesc(log10(abs(L)))
colorbar, set(gca,'CLim',[-15,0]); axis square

figure;
plot(eig(full(A)),'o');
hold on;
plot(eig(L(1:end-1,:),K(1:end-1,:)),'x');
T1 = L(1:end-1,:)/K(1:end-1,:);
plot(eig(T1),'s');
T2 = V(:,1:end-1)'*A*V(:,1:end-1);
plot(eig(T2),'d');

figure;
spy(abs(T1)>1e-12);

[Z,D] = eig(L(1:end-1,:),K(1:end-1,:));
norm(A*V(:,1:end-1)*K(1:end-1,:)*Z(:,1) - D(1,1)*V(:,1:end-1)*K(1:end-1,:)*Z(:,1))