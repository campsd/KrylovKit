clear all; close all; clc;
%% option 1: http://www.cise.ufl.edu/research/sparse/matrices/GHS_indef/spmsrtls.html
load('/home/daanc/phd/Copy/Software/matlab/UFget/mat/GHS_indef/spmsrtls.mat');
eigA = load('eigA_spmsrtls.dat');
%% option 2: http://www.cise.ufl.edu/research/sparse/matrices/Dehghani/light_in_tissue.html
load('/home/daanc/phd/Copy/Software/matlab/UFget/mat/Dehghani/light_in_tissue.mat')
eigA = load('eigA_light_in_tissue.dat');
%% option 3: http://www.cise.ufl.edu/research/sparse/matrices/DNVS/thread.html
load('/home/daanc/phd/Copy/Software/matlab/UFget/mat/DNVS/thread.mat')
eigA = load('eigA_thread.dat');
%% option 4: http://www.cise.ufl.edu/research/sparse/matrices/Rommes/bips07_1998.html
load('/home/daanc/phd/Copy/Software/matlab/UFget/mat/Rommes/bips07_1998.mat');
eigA = load('eigA_bips07_1998.dat');
%%
A = Problem.A;
figure;spy(A)
v=ones(size(A,1),1); v=v/norm(v);
tic;
A*v;
toc;
%%
global Pamd  Lfacamd Ufacamd pfacamd;
tic;
Pamd = amd(A);
toc;
%%
tic;
[Lfacamd,Ufacamd,pfacamd] = lu(A(Pamd,Pamd),'vector');
toc;
%%
tic;
sol2(Pamd) = mldivide(Ufacamd,mldivide(Lfacamd,v(Pamd(pfacamd))));
toc;
%%
global Lfac Ufac pfac;
tic;
[Lfac,Ufac,pfac] = lu(A,'vector');
toc;
%%
n = size(A,1);
V = zeros(n,1);
V(:,1) = v/norm(v,2);
KLrot = zeros(2,0); KLidx = zeros(1,0);
KR = zeros(1,0); LR = zeros(1,0);
pattern = [+1,-1,+1,-1,-1,+1];
pattern = repmat(pattern,1,10);
[V,KLrot,KLidx,KR,LR] = CTEK(@funcpos,@funcneg,V,KLrot,KLidx,KR,LR,pattern);
[ PC ] = CTPC( KLrot, KLidx, KR, LR, [] ); 
eigPC = eig(PC);
%%
figure;
plot(eigA(:,1),eigA(:,2),'.')
hold on
plot(real(eigPC),imag(eigPC),'o');