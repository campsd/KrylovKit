clc; clear all; close all;
format short e;
%%
n = 5000;
global A
A = olmstead( 0.6, n );
%% Factorize the matrix
global Lfac Ufac pfac;
tic;
[Lfac,Ufac,pfac] = lu(A,'vector');
toc;
% eigenvalues
% tic;
% eigA = eig(full(A));
% toc;
load('eigA_olmstead10000.mat');
lambda = [-0.880980204206560 + 1.776059795576901i; -0.880980204206560 - 1.776059795576901i];
%% Set up and run extended Krylov method
s1 = [-1]; % 100%
s2 = [1 -1 -1 -1]; %75
s3 = [1 -1 -1]; %66%
s4 = [1 -1]; %50%
s5 = [1 1 -1]; %33%
s6 = [1 1 1 -1]; %25%
s7 = [1 1 1 1 1 -1]; %20%
s8 = [1 1 1 1 1 1 -1]; %17%
%s7 = [1]; % 0%
%s = [1 -1 -1 1];
s = s8;
m = 20;
p = 10;
v = ones(size(A,1),1);
info = true;
% convergence tolerance
tol = 1e-8;
[ rval, res_init, err_init, res_rest, err_rest, nbRest, last_sz ] = CTIREK( @funcpos, @funcneg, v, s, m, p, lambda, tol, info );
%[ rval, conv_init, conv_restart, ls ] = CTIREK( @funcpos,@funcneg, s, m, k, l, p);

%% Eigenvalue plot
eigA = evalin('base','eigA');
figure;
plot(real(eigA),imag(eigA),'.');
hold on;
plot(real(rval),imag(rval),'o');
%xlim([-6 1])

%% Convergence plot

%% write to file
% downsampled spectrum
fileID = fopen('/home/daanc/phd/Copy/Writing/ExtendedKrylov/data/numexp2/eigA.dat','w');
fprintf('x\ty\n');
for i=1:size(eigA,1)
    if (abs(imag(eigA(i))) > 0) || (mod(i,40) == 1)
    fprintf(fileID','%12.6f %12.6f\n',[real(eigA(i)) imag(eigA(i))]');
    end
end
fclose(fileID);
% righ side of spectrum
I = find(real(eigA)>-7);
eigAr = eigA(I);
fileID = fopen('/home/daanc/phd/Copy/Writing/ExtendedKrylov/data/numexp2/eigAr.dat','w');
fprintf(fileID,'x\ty\n');
for i=1:size(eigAr,1)
    if (i < 30)|| (abs(imag(eigAr(i))) > 0) || (mod(i,40) == 1)
    fprintf(fileID','%12.6f %12.6f\n',[real(eigAr(i)) imag(eigAr(i))]');
    end
end
fclose(fileID);