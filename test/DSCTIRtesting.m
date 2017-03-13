clear all; close all; clc;
format short e;
%%
m = 484; B = 5.45; C = 2; Du = 0.004; Dv = 0.008; L = 1;
A = brusselator(m,C,B,Dv,Du,L);
A = sparse(A);
P = amd(A); A2 = A(P,P);
global Lfac Ufac pfac;
[Lfac,Ufac,pfac] = lu(A,'vector');
v=ones(size(A,1),1); v=v/norm(v);
n = size(A,1);
eigA = eig(full(A));
%%
V = zeros(n,1);
V(:,1) = v/norm(v,2);
KLrot = zeros(2,0); KLidx = zeros(1,0);
KR = zeros(1,0); LR = zeros(1,0);
global Lfac Ufac pfac;
[Lfac,Ufac,pfac] = lu(A,'vector');

pattern = [1,-1,-1,-1,1];
selection = repmat(pattern,[1,2]);
%new
[V,KLrot,KLidx,KR,LR] = CTEK(@funcpos,@funcneg,V,KLrot,KLidx,KR,LR,selection);
[K,L] = CONS_CTEK_PENCIL(KLrot,KLidx,KR,LR);
norm(A*V*K-V*L,'fro')
%%
mu = zeros(2,1);
%mu(1) = randn+randn*1i;
mu(1) = -1.509000422661994 + 2.297124747517389*1i;
mu(2) = conj(mu(1));
%mu(1) = -1.4524e+03;
%mu(2) = -2.4850e+03;
[V,KLrot,KLidx,KR,LR] = DSCTIR(A,V,KLrot,KLidx,KR,LR,mu);
[K,L] = CONS_CTEK_PENCIL(KLrot,KLidx,KR,LR);
norm(A*V*K-V*L,'fro')
[ PC ] = CTPC( KLrot, KLidx, KR, LR, [] );
figure, colormap('summer')
imagesc(log10(abs(V'*A*V)))
colorbar, set(gca,'CLim',[-15,0]); axis square

figure, colormap('summer')
imagesc(log10(abs(PC)))
colorbar, set(gca,'CLim',[-15,0]); axis square

figure,
plot(eigA,'.');
hold on
eigPC = eig(PC);
plot(real(eigPC),imag(eigPC),'o');
%%
[K,L] = CONS_CTEK_PENCIL(KLrot(:,1),KLidx(1),KR(1:2,1),LR(1:2,1));
norm(A*V(:,1:2)*K-V(:,1:2)*L,'fro')
%% OBSOLETE CODE
% perturbing CT
% 		if (KLidx(1) == 0) % First rotation at K (system solve)
%             tmp = KR(1:3,1:2);
%             tmp(2:3,:) = CreateRotMat(KLrot(:,2))*tmp(2:3,:);
%             tmp(1:2,:) = CreateRotMat(KLrot(:,1))*tmp(1:2,:);
%             L2(1) = LR(1,1)^2;
%             LK(1) = LR(1,1)*tmp(1,1); LK(2) = LR(1,1)*tmp(2,1);
%             K2(1) = tmp(1,1)^2 + tmp(2,1)*tmp(1,2);
%             K2(2) = tmp(1,1)*tmp(2,1) + tmp(2,1)*tmp(2,2);
%             K2(3) = tmp(3,2)*tmp(2,1);
% 		elseif (KLidx(1) == 1) % First rotation at L (multiplication)
%             tmp = LR(1:3,1:2);
%             tmp(2:3,:) = CreateRotMat(KLrot(:,2))*tmp(2:3,:);
%             tmp(1:2,:) = CreateRotMat(KLrot(:,1))*tmp(1:2,:);
% 			K2(1) = KR(1,1)^2;
%             LK(1) = tmp(1,1)*KR(1,1); LK(2) = KR(1,1)*tmp(2,1);
%             L2(1) = tmp(1,1)^2 + tmp(2,1)*tmp(1,2);
%             L2(2) = tmp(1,1)*tmp(2,1) + tmp(2,1)*tmp(2,2);
%             L2(3) = tmp(3,2)*tmp(2,1);
%         end

% This may not always be valid? correct!
%         [tmpK,tmpL] = CONS_CTEK_PENCIL(KLrot(:,1:2), KLidx(1:2), KR(1:3,1:3), LR(1:3,1:3));
%         L2 = tmpL*tmpL(:,1);
%         LK = tmpL*tmpK(:,1);
%         K2 = tmpK*tmpK(:,1);
%         R = L2 - (mu(1,kk) + mu(2,kk))*LK + (mu(1,kk)*mu(2,kk))*K2;