% In this test we construct an (artificial) real matrix with a bunch of
% eigenvalues along the negative real axis and a couple of complex
% conjugated eigenvalues on the imaginary axis. The idea is that an inverse
% operation has a hard time to find the imaginary eigenvalues but an MV
% doesn't
clc; clear all; close all;
format short e
%% Set-up matrix
tail = linspace(-100,-1,100);
ccpair_real = 0;
ccpair_imag = 3000;
A = sparse(zeros(length(tail)+2,length(tail)+2));
for i=1:length(tail)
    A(i,i) = tail(i);
end
A(end-1,end-1) = ccpair_real;
A(end,end) = ccpair_real;
A(end-1,end) = ccpair_imag;
A(end,end-1) = -ccpair_imag;
eigA = eig(full(A));

%% Set-up Krylov
n = size(A,1);
global Lfac Ufac pfac;
[Lfac,Ufac,pfac] = lu(A,'vector');
v=ones(n,1); v=v/norm(v);
V = zeros(n,1);
V(:,1) = v/norm(v,2);
KLrot = zeros(2,0); KLidx = zeros(1,0);
KR = zeros(1,0); LR = zeros(1,0);
% We choose a sequence of negative operations, this should do bad in
% approximating the rightmost eigenvalues
pattern = -ones(1,6);
pattern(end) = 1;
%pattern = [-1,1,-1,1,-1,1];
[V,KLrot,KLidx,KR,LR] = CTEK(@funcpos,@funcneg,V,KLrot,KLidx,KR,LR,pattern);
[ PC ] = CTPC( KLrot, KLidx, KR, LR, [] ); 
eigPC1 = eig(PC);
% We filter out the 6 left most eigenvalues
shifts = sort(eigPC1);
shifts = shifts(1:4);
[V,KLrot,KLidx,KR,LR] = CTIR(V,KLrot,KLidx,KR,LR,shifts);
[ PC ] = CTPC( KLrot, KLidx, KR, LR, [] ); 
eigPCfilt = eig(PC);
% We expand the space to its original size, but with positive operations
pattern = ones(1,4);
[V,KLrot,KLidx,KR,LR] = CTEK(@funcpos,@funcneg,V,KLrot,KLidx,KR,LR,pattern);
[ PC ] = CTPC( KLrot, KLidx, KR, LR, [] ); 
eigPC2 = eig(PC);

%% determine residual
res_PC1_1 = min(abs( ccpair_real + 1i * ccpair_imag - eigPC1))
res_PC1_2 = min(abs( ccpair_real - 1i * ccpair_imag - eigPC1))

res_PCfilt_1 = min(abs( ccpair_real + 1i * ccpair_imag - eigPCfilt))
res_PCfilt_2 = min(abs( ccpair_real - 1i * ccpair_imag - eigPCfilt))

res_PC2_1 = min(abs( ccpair_real + 1i * ccpair_imag - eigPC2))
res_PC2_2 = min(abs( ccpair_real - 1i * ccpair_imag - eigPC2))
%% plot data
figure;
plot(real(eigA),imag(eigA),'o');
figure;
plot(real(eigPC1),imag(eigPC1),'d');
figure;
plot(real(eigPCfilt),imag(eigPCfilt),'x');
figure;
plot(real(eigPC2),imag(eigPC2),'^');

%% write data to file
fileID = fopen('/home/daanc/phd/Copy/Writing/ExtendedKrylov/data/numexp1/eigA.dat','w');
fprintf(fileID','%12.6f %12.6f\n',[real(eigA) imag(eigA)]');
fclose(fileID);

fileID = fopen('/home/daanc/phd/Copy/Writing/ExtendedKrylov/data/numexp1/eigPC1.dat','w');
fprintf(fileID','%12.6f %12.6f\n',[real(eigPC1) imag(eigPC1)]');
fclose(fileID);

fileID = fopen('/home/daanc/phd/Copy/Writing/ExtendedKrylov/data/numexp1/eigPCfilt.dat','w');
fprintf(fileID','%12.6f %12.6f\n',[real(eigPCfilt) imag(eigPCfilt)]');
fclose(fileID);

fileID = fopen('/home/daanc/phd/Copy/Writing/ExtendedKrylov/data/numexp1/eigPC2.dat','w');
fprintf(fileID','%12.6f %12.6f\n',[real(eigPC2) imag(eigPC2)]');
fclose(fileID);
