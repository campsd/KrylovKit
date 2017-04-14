% unit tests for functionality in "standkrylov"
% last edit: April 13, 2017
ntest = 100;
%% Test DNS_SK
%
n = 1000;
A = sprand(n,n,0.05);
mv = @(x) A*x;
% test 1: assert that start vector will be normalized
v = ones(n+1,1);
V = v; H = zeros(1,0);
m = 1;
[ V, H ] = DNS_SK( mv, V, H, m );
assert(abs(norm(V,2))-1<2*eps,'failed at: DNS_SK case 1.1');
assert(size(V,2)==1,'failed at: DNS_SK case 1.1');
%% Test DNS_SK_BLK

%% Test CT_SK

%% Test CT_SK_BLK

%% Test CT_SK_HESS

%% Test CT_SK_HESS_BLK

%% Test CT_SK_IR_SS

%% Test CT_SK_TO_EK_LEFT

%% Test CT_SK_TO_EK_RIGHT

%% Test CT_SK_TO_EK_RIGHT_BLK

