
A = randn(10,10);
v = randn(10,1);
s = [1 1 1 1 1];
global A Lfac Ufac pfac;
tic;
[Lfac,Ufac,pfac] = lu(A,'vector');
toc;
%%
n = size(A,1);
V = zeros(n,1);
V(:,1) = v/norm(v,2);
KLrot = zeros(2,0); KLidx = zeros(1,0);
KR = zeros(1,0); LR = zeros(1,0);
[V,KLrot, KLidx,KR,LR] = CTEK(@funcpos,@funcneg,V,KLrot,KLidx,KR,LR,s);
%%
I = eye(length(s)+1);
ei = I(:,end);
%I = I(:,end);
for i=length(s):-1:1
I(i:i+1,:) = CreateRotMat(KLrot(:,i)) *  I(i:i+1,:);
end

%%
A = [1 2 3 4 5;
     6 7 8 9 10;
     11 12 13 14 15;
     16 17 18 19 20;
     21 22 23 24 25];
I = eye(5);
ei = I(:,end);
I(:,end) = [];
ei*ei'
A*ei*ei'*A