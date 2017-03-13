clc; clear all; close all;
% This problem is related to a fluid flow problem. It is also studied in:
% An implicitly restarted rational Krylov strategy for Lyapunov inverse
% iteration - K. Meerbergen (2015)

%%
global A B;
A = mmread('howard_A7800c.mtx') ;
B = mmread('howard_Bc.mtx') ;
eigAB = load('eigAB_cavity.dat');
lambda = [-0.005135085169115 + 2.6984472335173851i; -0.005135085169115 - 2.6984472335173851i];
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

%% Set up and run extended Krylov method
s0 = [1 -1 -1 -1 -1]; %80
s1 = [1 -1 -1 -1]; %75
s2 = [1 -1 -1]; %66%
s = s2;
m = 100;
p = 50;
l = 12;
v = ones(size(A,1),1);
%v = randn(size(A,1),1);
info = true;
GEP = true;
% convergence tolerance
tol = 1.5e-8; %s0
%tol = 5e-5; %s1
%tol = 4e-4; %s2
[ rval, res_init, err_init, res_rest, err_rest, nbRest, last_sz ] = CTIREK( @funcpos_generalized, @funcneg_generalized, v, s, m, p, l, lambda, tol, info, GEP );
%[ rval, conv_init, conv_restart, ls ] = CTIREK( @funcpos,@funcneg, s, m, k, l, p);

%% Plot spectrum and ritz values
figure;
plot(eigAB(:,1),eigAB(:,2),'.');
hold on;
plot(real(rval),imag(rval),'o');
xlim([-2.5 0.5])
ylim([-5 5])

%% Convergence plots of eigenvalues of interest

for i=1:2
    ri = res_init{i};
    ei = err_init{i};
    figure;
    subplot(121)
    hRi = semilogy(ri(:,1),ri(:,2)); 
    
    % format hRi
       set(hRi                        , ...
      'LineWidth'       , 1.8           , ...
      'LineStyle'       , '-'        , ...
      'Color'       , [0 0 0]           , ...
      'Marker'          , 'o'         , ...
      'MarkerSize'      , 1           , ...
      'MarkerEdgeColor' , 'none'      , ...
      'MarkerFaceColor' , [.5 .5 .5] );
  
    hold on
    for j=1:nbRest
       rr = res_rest{i,j};
       hRr = semilogy(rr(:,1),rr(:,2));
       
       % format hRr
       set(hRr                        , ...
      'LineWidth'       , 0.8         , ...
      'LineStyle'       , '-.'        , ...
      'Color'       , [.5 .5 .5]           , ...
      'Marker'          , 'o'         , ...
      'MarkerSize'      , 1           , ...
      'MarkerEdgeColor' , 'none'      , ...
      'MarkerFaceColor' , [.5 .5 .5] );
    end
    subplot(122)
    hEi = semilogy(ei(:,1),ei(:,2)); 
    hold on
    
   % format hEi
   set(hEi                        , ...
  'LineWidth'       , 1.8           , ...
  'LineStyle'       , '-'        , ...
  'Color'       , [0 0 0]           , ...
  'Marker'          , 'o'         , ...
  'MarkerSize'      , 1           , ...
  'MarkerEdgeColor' , 'none'      , ...
  'MarkerFaceColor' , [.5 .5 .5] );
  
  
    for j=1:nbRest
       er = err_rest{i,j};
       hEr = semilogy(er(:,1),er(:,2));
       
       % format hEr
       set(hEr                        , ...
      'LineWidth'       , 0.8           , ...
      'LineStyle'       , '-.'        , ...
      'Color'       , [.5 .5 .5]           , ...
      'Marker'          , 'o'         , ...
      'MarkerSize'      , 1           , ...
      'MarkerEdgeColor' , 'none'      , ...
      'MarkerFaceColor' , [.5 .5 .5] );
  
    end
end
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
%% save even smaller part
figure;
plot(eigAB(:,1),eigAB(:,2),'.');
xlim([-0.05 -0.004])
ylim([-5 5])

eigABs = [];
for i = 1:size(eigAB,1)
    if (abs(eigAB(i,1)) < 0.05) && (abs(eigAB(i,2)) < 4)
        eigABs(end+1,:) = eigAB(i,:);
    end
end
fileID = fopen('/home/daanc/phd/Copy/Writing/ExtendedKrylov/data/numexp3/eigABr.dat','w');
fprintf(fileID','%12.6f %12.6f\n',[eigABs(:,1) eigABs(:,2)]');
fclose(fileID);

%% OBSOLETE
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
