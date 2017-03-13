clc; clear all; close all;
% This problem is related to a fluid flow problem. It is also studied in:
% An implicitly restarted rational Krylov strategy for Lyapunov inverse
% iteration - K. Meerbergen (2015)

%%
global A B;
A = mmread('howard_A7800c.mtx') ;
B = mmread('howard_Bc.mtx') ;
eigAB = load('eigAB_cavity.dat');
%% eigs
tic;
eigsAB = eigs(A,B,6,'LR');
toc;
%% Factorize the matrices
global LfacA UfacA pfacA;
tic;
[LfacA,UfacA,pfacA] = lu(A,'vector');
toc;
global LfacB UfacB pfacB;
tic;
[LfacB,UfacB,pfacB] = lu(B,'vector');
toc;
%% Initialize Krylov
tol = sqrt(eps);
n = size(A,1);
v=ones(n,1); v=v/norm(v);
V = zeros(n,1);
V(:,1) = v/norm(v,2);
KLrot = zeros(2,0); KLidx = zeros(1,0);
KR = zeros(1,0); LR = zeros(1,0);
%% Choose selection vector
pattern = [];
for i=1:60
    if mod(i,2) == 0
        pattern(i) = 1;
    else
        pattern(i) = -1;
    end
end
%%
pattern = [-1 1 -1 -1];
[ rval ] = CTIREK( @funcpos_generalized,@funcneg_generalized, pattern, 140, 10,100, 30);
%[ rval ] = CTIREK_CD( @funcpos_generalized,@funcneg_generalized, pattern, 40, 6,10, 24);
fprintf('Real part of rightmost eigenvalue: \n');
fprintf('True value: %12.6f \t Found value: %12.6f\n',max(eigAB(:,1)),max(real(rval)));
%%
tic;
[V,KLrot,KLidx,KR,LR] = CTEK(@funcpos_generalized,@funcneg_generalized,V,KLrot,KLidx,KR,LR,pattern);
toc;
[ PC ] = CTPC( KLrot, KLidx, KR, LR, [] );
[eigvecPC, eigPC] = eig(PC);
eigPC = diag(eigPC);

%%
eigAB = evalin('base','eigAB');
figure;
plot(eigAB(:,1),eigAB(:,2),'.');
hold on;
plot(real(rval),imag(rval),'o');
xlim([-2.5 0.5])
ylim([-5 5])

%%

figure;
plot(eigAB(:,1),eigAB(:,2),'.');
hold on;
plot(real(rval(ind_desired)),imag(rval(ind_desired)),'o');
xlim([-2.5 0.5])
ylim([-5 5])

%% save part of the spectrum for paper
eigABs = [];
for i = 1:size(eigAB,1)
    if (abs(eigAB(i,1)) < 2.5) && (abs(eigAB(i,2)) < 4)
        eigABs(end+1,:) = eigAB(i,:);
    end
end
fileID = fopen('/home/daanc/phd/Copy/Writing/ExtendedKrylov/data/numexp3/eigAB.dat','w');
fprintf(fileID','%12.6f %12.6f\n',[eigABs(:,1) eigABs(:,2)]');
fclose(fileID);