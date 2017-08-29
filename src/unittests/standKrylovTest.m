% unit tests for functionality in "standkrylov"
% last edit: June 30, 2017
ntest = 100;
%% Test DNS_SK
%
n = 1000;
A = sprand(n,n,0.05);
mv = @(x) A*x;
% test 1: assert that start vector will be normalized, compute no vectors
V = randn(n,1);
m =0;
[ V, ~ ] = DNS_SK( mv, V, m );
assert((abs(norm(V,2))-1)<10*eps,'failed at: DNS_SK case 1.1');
assert(size(V,2)==1,'failed at: DNS_SK case 1.2');
% test 2: take steps, validate size, validate recurrence
V = randn(n,1);
m =10;
[ V, H ] = DNS_SK( mv, V, m );
assert(size(V,2)==m+1,'failed at: DNS_SK case 2.1');
assert(norm(mv(V(:,1:m))-V*H,'fro')/norm(V*H,'fro') < 2*eps, 'failed at: DNS_SK case 2.2');
% test 3: enlarge current basis, validate size, validate recurrence
cs = size(V,2);
m = 120;
[ V, H ] = DNS_SK( mv, V, H, m );
assert(size(V,2)==cs+m,'failed at: DNS_SK case 2.1');
assert(norm(mv(V(:,1:cs+m-1))-V*H,'fro')/norm(V*H,'fro') < 2*eps, 'failed at: DNS_SK case 2.2');


%% Test DNS_SK_BLK
% test 1: no steps
n = 1000;
A = sprand(n,n,0.05);
mv = @(x) A*x;
bs = 3; m = 0;
V = (randn(n,bs));
[ V, ~ ] = DNS_SK_BLK( mv, V, bs, m );
assert((abs(norm(V,2))-1)<10*eps,'failed at: DNS_SK_BLK case 1.1');
% test 2: some steps
m = 5;
[ V, H ] = DNS_SK_BLK( mv, V, bs, m );
assert((abs(norm(V,2))-1)<10*eps,'failed at: DNS_SK_BLK case 2.1');
assert(size(V,2)==(m+1)*bs,'failed at: DNS_SK_BLK case 2.2');
assert(norm(mv(V(:,1:m*bs))-V*H,'fro')/norm(V*H,'fro') < 2*eps, 'failed at: DNS_SK_BLK case 2.3');
% test 3: enlarge the subspace
m2 = 17;
[ V, H ] = DNS_SK_BLK( mv, V, H, bs, m2 );
assert((abs(norm(V,2))-1)<10*eps,'failed at: DNS_SK_BLK case 3.1');
assert(size(V,2)==(m+m2+1)*bs,'failed at: DNS_SK_BLK case 3.2');
assert(norm(mv(V(:,1:(m+m2)*bs))-V*H,'fro')/norm(V*H,'fro') < 2*eps, 'failed at: DNS_SK_BLK case 3.3');

%% Test CT_SK, CT_TO_DNS_SK_HESS, DNS_TO_CT_SK_HESS
n = 1000;
A = sprand(n,n,0.05);
mv = @(x) A*x;
V = randn(n,1);
% take no steps
m =0;
[ V, ~ ] = CT_SK( mv, V, m );
assert((abs(norm(V,2))-1)<10*eps,'failed at: CT_SK case 1.1');
assert(size(V,2)==1,'failed at: CT_SK case 1.2');
% take some steps and verify result is the same as DNS
m=10;
[ Vref, Href ] = DNS_SK( mv, V, m );
[ V, Hrot, HR ] = CT_SK( mv, V, m );
assert(norm(V-Vref)==0,'faild at: CT_SK case 2.1')
H = CT_TO_DNS_SK_HESS(Hrot,HR);
assert(norm(H-Href,'fro')/norm(H,'fro') < eps, 'failed at: CT_SK case 2.2')
[Hrotref,HRref] = DNS_TO_CT_SK_HESS(Href);
assert(norm(Hrot-Hrotref,'fro')/norm(Hrot,'fro') < eps, 'failed at: CT_SK case 2.3')
assert(norm(HR-HRref,'fro')/norm(HR,'fro') < eps, 'failed at: CT_SK case 2.4')
% enlarge the space
m=50;
[ Vref, Href ] = DNS_SK( mv, Vref, Href, m );
[ V, Hrot, HR ] = CT_SK( mv, V, Hrot, HR, m );
assert(norm(V-Vref)==0,'faild at: CT_SK case 3.1')
H = CT_TO_DNS_SK_HESS(Hrot,HR);
assert(norm(H-Href,'fro')/norm(H,'fro') < eps, 'failed at: CT_SK case 3.2')
[Hrotref,HRref] = DNS_TO_CT_SK_HESS(Href);
assert(norm(Hrot-Hrotref,'fro')/norm(Hrot,'fro') < eps, 'failed at: CT_SK case 3.3')
assert(norm(HR-HRref,'fro')/norm(HR,'fro') < eps, 'failed at: CT_SK case 3.4')
%% Test CT_SK_BLK, CT_TO_DNS_SK_HESS_BLK and DNS_TO_CT_SK_HESS_BLK
% test 1
n = 1000;
A = sprand(n,n,0.05);
mv = @(x) A*x;
bs = 8; m = 0;
V = (randn(n,bs));
[ V, ~ ] = CT_SK_BLK( mv, V, bs, m );
assert((abs(norm(V,2))-1)<10*eps,'failed at: CT_SK_BLK case 1.1');
% test 2: take some steps and verify result is same as DNS
m = 5;
[ Vref, Href ] = DNS_SK_BLK( mv, V, bs, m );
[ V, Hrot, Hrow, HR ] = CT_SK_BLK( mv, V, bs, m );
assert(norm(V-Vref)==0,'faild at: CT_SK_BLK case 2.1')
[ Hrotref, Hrowref, HRref ] = DNS_TO_CT_SK_HESS_BLK( Href );
H = CT_TO_DNS_SK_HESS_BLK(Hrot,Hrow,HR);
assert(norm(Hrowref(:)-Hrow(:),'fro')==0,'failed at: CT_SK_BLK case 2.2')
assert(norm(HR-HRref,'fro')/norm(HR,'fro') < eps,'failed at: CT_SK_BLK case 2.3')
assert(norm(Hrot(:)-Hrotref(:),'fro')/norm(Hrot(:),'fro') < eps,'failed at: CT_SK_BLK case 2.4')
assert(norm(H-Href,'fro')/norm(Href,'fro') < 3*eps,'failed at: CT_SK_BLK case 2.5')
% test 3: enlarge
m = 7;
[ Vref, Href ] = DNS_SK_BLK( mv, Vref, Href, bs, m );
[ V, Hrot, Hrow, HR ] = CT_SK_BLK( mv, V, Hrot, Hrow, HR, bs, m );
assert(norm(V-Vref)==0,'faild at: CT_SK_BLK case 3.1')
[ Hrotref, Hrowref, HRref ] = DNS_TO_CT_SK_HESS_BLK( Href );
H = CT_TO_DNS_SK_HESS_BLK(Hrot,Hrow,HR);
assert(norm(Hrowref(:)-Hrow(:),'fro')==0,'failed at: CT_SK_BLK case 3.2')
assert(norm(HR-HRref,'fro')/norm(HR,'fro') < eps,'failed at: CT_SK_BLK case 3.3')
assert(norm(Hrot(:)-Hrotref(:),'fro')/norm(Hrot(:),'fro') < eps,'failed at: CT_SK_BLK case 3.4')
assert(norm(H-Href,'fro')/norm(Href,'fro') < 3*eps,'failed at: CT_SK_BLK case 3.5')

%% Test CT_SK_TO_EK_LEFT
m
%% Test CT_SK_TO_EK_RIGHT

%% Test CT_SK_TO_EK_RIGHT_BLK

%% Test CT_SK_IR_SS
n = 1000;
A = sprand(n,n,0.05);
mv = @(x) A*x;
m = 30;
V = randn(n,1);
[V,Hrot,HR] = CT_SK(mv,V,m);
H = CT_TO_DNS_SK_HESS(Hrot,HR);
eigH = eig(H(1:m,:));
eigH = cplxpair(eigH);
realstart = find(imag(eigH)==0,1,'first');
if ~isempty(realstart)
    % remove a single real Ritz value from the subspace
    mu = eigH(realstart);
    tmp = eigH; tmp(realstart) = [];
    mindif = min(abs(tmp-mu));
    [V,Hrot,HR] = CT_SK_IR_SS(V,Hrot,HR,mu);
    H = CT_TO_DNS_SK_HESS(Hrot,HR);
    eigHf = eig(H(1:m-1,:));
    assert(norm(mv(V(:,1:m-1))-V*H,'fro')/norm(V*H,'fro') < 3*eps,'failed at: CT_SK_IR_SS 1.1')
    assert(abs(min(abs(eigHf-mu))-mindif)/abs(mindif) < 100*eps,'failed at: CT_SK_IR_SS 1.2')
    assert(norm(sort(tmp)-sort(eigHf),'fro')/norm(eigHf,'fro') < 100*eps,'failed at: CT_SK_IR_SS 1.3')
else
    warning('test CT_SK_IR 1 not executed')
end
if realstart > 1
    % remove a single complex Ritz value from the subspace
    idx = ceil(rand*(realstart-1));
    mu = eigHf(idx);
    tmp = eigHf; tmp(idx) = [];
    mindif = min(abs(tmp-mu));
    [V,Hrot,HR] = CT_SK_IR_SS(V,Hrot,HR,mu);
    H = CT_TO_DNS_SK_HESS(Hrot,HR);
    eigHf2 = eig(H(1:m-2,:));
    assert(norm(mv(V(:,1:m-2))-V*H,'fro')/norm(V*H,'fro') < 3*eps,'failed at: CT_SK_IR_SS 2.1')
    assert(abs(min(abs(eigHf2-mu))-mindif)/abs(mindif) < 100*eps,'failed at: CT_SK_IR_SS 2.2')
    assert(norm(sort(abs(tmp))-sort(abs(eigHf2)),'fro')/norm(eigHf2,'fro') < 100*eps,'failed at: CT_SK_IR_SS 2.3')
else
    warning('test CT_SK_IR 2 not executed')
end
% remove multiple Ritz values from the subspace
v = randn(n,1);
V = v;
[V,Hrot,HR] = CT_SK(mv,V,m);
H = CT_TO_DNS_SK_HESS(Hrot,HR);
eigH = eig(H(1:m,:));
eigH = cplxpair(eigH);
realstart = find(imag(eigH)==0,1,'first');
if ~isempty(realstart) && realstart > 1 && realstart < m
    mu = eigH(realstart-1:realstart+1);
    tmp = eigH; tmp(realstart-1:realstart+1) = [];
    mindif1 = min(abs(tmp-mu(1)));
    mindif2 = min(abs(tmp-mu(2)));
    mindif3 = min(abs(tmp-mu(3)));
    [V,Hrot,HR] = CT_SK_IR_SS(V,Hrot,HR,mu);
    H = CT_TO_DNS_SK_HESS(Hrot,HR);
    eigHf = eig(H(1:m-3,:));
    assert(norm(mv(V(:,1:m-3))-V*H,'fro')/norm(V*H,'fro') < 10*eps,'failed at: CT_SK_IR_SS 3.1')
    assert(abs(min(abs(eigHf-mu(1)))-mindif1)/abs(mindif1) < 100*eps,'failed at: CT_SK_IR_SS 3.2.1')
    assert(abs(min(abs(eigHf-mu(2)))-mindif2)/abs(mindif2) < 100*eps,'failed at: CT_SK_IR_SS 3.2.2')
    assert(abs(min(abs(eigHf-mu(3)))-mindif3)/abs(mindif3) < 100*eps,'failed at: CT_SK_IR_SS 3.2.3')
    assert(norm(sort(abs(tmp))-sort(abs(eigHf)),'fro')/norm(eigHf,'fro') < 100*eps,'failed at: CT_SK_IR_SS 3.3')
else
    warning('test CT_SK_IR_SS 3 not executed')
end
%% Test CT_SK_IR_DS

