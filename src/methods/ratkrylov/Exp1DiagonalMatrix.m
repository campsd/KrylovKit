%% start with an empty workspace
clear all; close all; clc;
%% define the problem:
% We are going to define a diagonal matrix with n eigenvalues
% spread out over a couple of clusters. The idea would be to 
% determine the "edges" of the clusters with the implicitly 
% restart algorithm
global A LfacA UfacA pfacA B LfacB UfacB pfacB;
n = 200; %problem size
clusters = [1 2; 3 4; 5 5.5]; % define cluster edges
w = [0.35 0.55 0.1]; %define cluster weights (should add up to one)
% create the diagonal
dg = [];
for i=1:size(clusters,1)
    dg = [dg; (clusters(i,2)-clusters(i,1)).*rand(round(w(i)*n),1) + clusters(i,1)];
end
eigAB = dg;
% randomly permute the diagonal
dg = dg(randperm(length(dg)));
% Create the matrix
A = diag(dg);
B = eye(n,n);
[LfacA,UfacA,pfacA] = lu(A,'vector');
[LfacB,UfacB,pfacB] = lu(B,'vector');
v = ones(n,1); v = v/norm(v,2);
% RKS variables
V = zeros(n,1);
V(:,1) = v;
Lrot = zeros(2,0);
KR = zeros(1,0); LR = zeros(1,0);
% Poles
npoles = 10;
wpoles = [0.5 0.4 0.1]; % weight of pole distribution
pl = []; Xi = {};
for i=1:size(clusters,1)
    pl = [pl; (clusters(i,2)-clusters(i,1)).*rand(round(wpoles(i)*npoles),1) + clusters(i,1)];
end
pl = pl(randperm(length(pl)));
for i =1:npoles
    Xi{i} = [1 pl(i) 0 1];
    %Xi{i} = Inf;
end
%Xi{npoles} = Inf;
T = eye(npoles+1,npoles);
%% Compute the RKS recursion
res = {}; theta = {};
for i=1:1
    [ V, KR, Lrot, LR ] = CTRKS( A, B, @funcpos_generalized, @funcneg_generalized, V, KR, Lrot, LR, {Xi{i}},T(1:i+1,i));
    K = KR;
    L = LR;
    for j=size(Lrot,2):-1:1
        L(j:j+1,:) = CreateRotMat(Lrot(:,j)) * L(j:j+1,:);
    end
    theta{i} = eig(L(1:end-1,:),K(1:end-1,:));
    res{i} = RitzResiduals(eigAB,theta{i});
    theta{i} = real(theta{i});
end
RitzPlot(1:npoles,res,theta);
%%
K = KR;
L = LR;
for i=length(Lrot):-1:1
    L(i:i+1,:) = CreateRotMat(Lrot(:,i)) * L(i:i+1,:);
end
norm(A*V*K-B*V*L,'fro')/norm(A*V*K,'fro')
%rval = eig(L(1:end-1,:)/K(1:end-1,:));
rval = eig(L(1:end-1,:),K(1:end-1,:));
%%
figure;
plot(real(eigAB),imag(eigAB),'.');
hold on
%%plot(real(rval),imag(rval),'o');
%plot(real(rval2),imag(rval2),'x');
%%
m = 30;
l = 12;
p = 16;
tol = 1e-5;
info = true;
[ rval, res, theta ] = CTIRRKS( A, B, v, Xi, m, p, l, @selectRightmost, eigAB, tol, info );
RitzPlot(1:length(res),res,theta);
[ rval, res, theta ] = CTIRRKS( A, B, v, Xi, m, p, l, @selectLeftmost, eigAB, tol, info );
RitzPlot(1:length(res),res,theta);