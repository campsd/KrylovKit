format short e;
m = 2904; B = 5.45; C = 2; Du = 0.004; Dv = 0.008; L = 1;
A = brusselator(m,C,B,Dv,Du,L);
A = sparse(A);
%v=ones(size(A,1),1); v=v/norm(v);
v=randn(size(A,1),1); v=v/norm(v);
global Lfacamd Ufacamd pfacamd Pamd;
global Lfac Ufac pfac;
tic;
[Lfac,Ufac,pfac] = lu(A,'vector');
toc
tic;
Pamd = amd(A);
[Lfacamd,Ufacamd,pfacamd] = lu(A(Pamd,Pamd),'vector');
toc;
tic;
sol1 = mldivide(Ufac,mldivide(Lfac,v(pfac)));
toc;
tic;
sol2(Pamd) = mldivide(Ufacamd,mldivide(Lfacamd,v(Pamd(pfacamd))));
toc;
norm(sol1-sol2')
%%
figure;
plot(sol1);
hold on
plot(sol2)

%% Krylov
n = size(A,1);
V = zeros(n,1);
V(:,1) = v/norm(v,2);
KLrot = zeros(2,0); KLidx = zeros(1,0);
KR = zeros(1,0); LR = zeros(1,0);
pattern = [-1,-1,+1,-1,-1,+1];
%[V,KLrot,KLidx,KR,LR] = CTEK(@funcpos,@funcneg,V,KLrot,KLidx,KR,LR,pattern);
[V,KLrot,KLidx,KR,LR] = CTEK(@funcpos,@funcnegamd,V,KLrot,KLidx,KR,LR,pattern);
        