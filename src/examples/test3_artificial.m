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
ccpair_imag = 25;
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
tol = sqrt(eps);
n = size(A,1);
global Lfac Ufac pfac;
[Lfac,Ufac,pfac] = lu(A,'vector');
v=ones(n,1); v=v/norm(v);
V = zeros(n,1);
V(:,1) = v/norm(v,2);
KLrot = zeros(2,0); KLidx = zeros(1,0);
KR = zeros(1,0); LR = zeros(1,0);
%% option 1 - all -1
pattern_init = -ones(1,7);
pattern_restart = -ones(1,6);
%% option 2 - all 1
pattern_init = ones(1,7);
pattern_restart = ones(1,6);
%% option 3 - alternating
pattern_init = [-1 1 -1 1 -1 1 -1];
pattern_restart = [-1 1 -1 1 -1 1];
%% option 4
pattern_init = [-1 -1 -1 1 1 1 1];
pattern_restart = [-1 -1 -1 1 1 1];
%% option 5
pattern_init = -ones(1,7);
pattern_restart_1 = [-1 -1 -1 -1 -1 -1];
pattern_restart_2 = [1 1 1 1 1 1];
%% We choose a sequence of negative operations, this should be bad in
% approximating the rightmost eigenvalues

res_PC1_1 = [];
res_PC1_2 = [];
for i=1:length(pattern_init)
    [V,KLrot,KLidx,KR,LR] = CTEK(@funcpos,@funcneg,V,KLrot,KLidx,KR,LR,pattern_init(i));
    %[ PC ] = CTPC( KLrot, KLidx, KR, LR, [] );
    PC = V'*A*V;
    eigPC1 = eig(PC);
    % determine residual
    res_PC1_1(i) = min(abs( ccpair_real + 1i * ccpair_imag - eigPC1))/abs(ccpair_real + 1i * ccpair_imag);
    res_PC1_2(i) = min(abs( ccpair_real - 1i * ccpair_imag - eigPC1))/abs(ccpair_real + 1i * ccpair_imag);
end
res = min(res_PC1_1);
restart = 1;
res_PC2_1 = [];
res_PC2_2 = [];
eigPC2 = eigPC1;
while res > tol,
    %% We filter the eigenvalues furthest from the wanted
    [~,I1] = min(abs( ccpair_real + 1i * ccpair_imag - eigPC2));
    [~,I2] = min(abs( ccpair_real - 1i * ccpair_imag - eigPC2));
    shifts = eigPC2;
    shifts(min(I1,I2)) = [];
    shifts(max(I1,I2)-1) = [];
    [V,KLrot,KLidx,KR,LR] = CTIR(V,KLrot,KLidx,KR,LR,shifts);
    PC = V'*A*V;
    %[ PC ] = CTPC( KLrot, KLidx, KR, LR, [] ); 
    eigPCfilt = eig(PC);
    %% We expand the space to its original size, but with positive operations
%     if restart > 1,
%         pattern_restart = pattern_restart_2;
%     else
%         pattern_restart = pattern_restart_1;
%     end
    for i=1:length(pattern_restart)
        [V,KLrot,KLidx,KR,LR] = CTEK(@funcpos,@funcneg,V,KLrot,KLidx,KR,LR,pattern_restart(i));
        PC = V'*A*V;
        %[ PC ] = CTPC( KLrot, KLidx, KR, LR, [] ); 
        eigPC2 = eig(PC);
        res_PC2_1(restart,i) = min(abs( ccpair_real + 1i * ccpair_imag - eigPC2))/abs(ccpair_real + 1i * ccpair_imag);
        res_PC2_2(restart,i) = min(abs( ccpair_real - 1i * ccpair_imag - eigPC2))/abs(ccpair_real + 1i * ccpair_imag);
    end
    res = min(res_PC2_1(restart,:));
    restart = restart + 1
end
%% determine residual
res_PC1_1 = min(abs( ccpair_real + 1i * ccpair_imag - eigPC1))
res_PC1_2 = min(abs( ccpair_real - 1i * ccpair_imag - eigPC1))

res_PCfilt_1 = min(abs( ccpair_real + 1i * ccpair_imag - eigPCfilt))
res_PCfilt_2 = min(abs( ccpair_real - 1i * ccpair_imag - eigPCfilt))

res_PC2_1 = min(abs( ccpair_real + 1i * ccpair_imag - eigPC2))
res_PC2_2 = min(abs( ccpair_real - 1i * ccpair_imag - eigPC2))
%% plot eigenvalues
figure;
plot(real(eigA),imag(eigA),'o');
figure;
plot(real(eigPC1),imag(eigPC1),'d');
figure;
plot(real(eigPCfilt),imag(eigPCfilt),'x');
figure;
plot(real(eigPC2),imag(eigPC2),'^');

%% plot convergence
initx = [2:8];
restartx = [3:8];
figure;
semilogy(initx, res_PC1_1)
hold on
for i=1:restart-1
    semilogy(restartx, res_PC2_1(i,:),'o-')
end

%% save convergence to file
filename = '/home/daanc/phd/Copy/Writing/ExtendedKrylov/data/numexp1/s4.dat';
fileID = fopen(filename,'w');
for i=1:length(initx)
    if i>length(restartx)
        i
        fprintf(fileID,'%d %23.16e\n',[initx(i) res_PC1_1(i)]');
    else
        i
        fprintf(fileID,'%d %23.16e %d %23.16e ',[initx(i) res_PC1_1(i) restartx(i) res_PC2_1(:,i)']');
        fprintf(fileID,'\n');
    end
end
fclose(fileID);
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
