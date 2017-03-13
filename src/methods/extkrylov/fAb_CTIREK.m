function [ fAb, res ] = fAb_CTIREK( A, b, pattern, kmax, maxrestarts, tol, funmat, funselect )
%FAB_CTIREK - f(A)b computation by means of CTIREK
%   A,b              mat,vec
%   pattern          pattern of operations (+1,-1)
%   kmax             max times pattern can be executed
%   maxrestarts      max number of restarts
%   tol              tolerance
%   funpos(v)        operation A*v
%   funneg(v)        operation A^-1*v
%   funmat(P)        matrix function to execute on projected matrix
%   funselect(rv)    select relevant Ritz values

% (1) Initialize
% LU factorization
global Lfac Ufac pfac;
[Lfac,Ufac,pfac] = lu(A,'vector');
% EK variables
V = zeros(size(A,1),1);
nrmb = norm(b);
b = b/nrmb;
V(:,1) = b;
KLrot = zeros(2,0); KLidx = zeros(1,0);
KR = zeros(1,0); LR = zeros(1,0);
% f(A)b approximation
notConverged = true;
restarts = 0;
PC = [];
xm = V(:,1);
res = [];

% (2) Iterate
while notConverged,
   for i=1:kmax 
        [V,KLrot,KLidx,KR,LR] = CTEK(@funcpos,@funcneg,V,KLrot,KLidx,KR,LR,pattern);
        [ PC ] = CTPC( KLrot, KLidx, KR, LR, [] ); % PC as argument doesn't work correctly!
        xmj = V(:,1:end-1)*funmat(PC)*V(:,1:end-1)'*b;
        
        %delta
        deltamj = norm(xmj-xm)/norm(xm);
        resj = abs(deltamj/(1-deltamj));
        res(end+1) = resj;
        if resj <= tol,
            notConverged = 0;
            break
        end
        xm = xmj;
   end
   
   if notConverged && restarts < maxrestarts, %restart
       shifts = funselect(PC);
       kmax = ceil(length(shifts)/length(pattern));
       [V,KLrot,KLidx,KR,LR] = CTIR(V,KLrot,KLidx,KR,LR,shifts);
       PC = [];
       restarts = restarts + 1;
   end
   
   if restarts >= maxrestarts,
       notConverged = 0;
   end
end

fAb = xmj;
end
