% unit tests for functionality in "common"
% last edit: April 13, 2017
ntest = 10000;
%% test CT_GIV
% The function to compute a core transformation from a 2-vector is tested
% on a few examples and the results are validated.
% case 1: the second element is already zero
[cos, sin, r] = CT_GIV( 1, 0);
assert((cos==1) && (sin==0) && (r==1),'failed at: CT_GIV case 1');
% case 2: the first element is zero
[cos, sin, r] = CT_GIV( 0, 1);
assert((cos==0) && (sin==1) && (r==1),'failed at: CT_GIV case 2');
% case 3: random 2-vec - test if norm is conserved
% case 4: test if c,s are correct
% REAL
for kk=1:ntest
v = randn(2,1);
[cos, sin, r] = CT_GIV( v(1), v(2));
assert(abs(norm(v)-r)<=5*eps,'failed at: CT_GIV case 3');
assert((cos^2+sin^2-1)<=5*eps,'failed at: CT_GIV case 4');
end
% COMPLEX
for kk=1:ntest
v = randn(2,1) + rand(2,1)*1i; 
[cos, sin, r] = CT_GIV( v(1), v(2));
assert(abs(norm(v)-r)<=10*eps,'failed at: CT_GIV case 5');
assert((cos^2+sin^2-1)<=10*eps,'failed at: CT_GIV case 6');
end
%% test CT_H
% The function to compute the Hermitian conjugate of a single or multiple
% core transformations is tested.
% case 1: a single, random rotation
% REAL
v = randn(2,1); 
[cos, sin, ~] = CT_GIV( v(1), v(2));
CH = CT_H([cos;sin]);
%assert
assert(CH(1) == conj(cos),'failed at: CT_H case 1');
assert(CH(2) == -sin,'failed at: CT_H case 2');
% COMPLEX
v = randn(2,1) + randn(2,1)*1i; 
[cos, sin, ~] = CT_GIV( v(1), v(2));
CH = CT_H([cos;sin]);
%assert
assert(CH(1) == conj(cos),'failed at: CT_H case 3');
assert(CH(2) == -sin,'failed at: CT_H case 4');
% case 2: multiple rotations
G = [];
for kk=1:ntest
    v = randn(2,1);
    [cos, sin, ~] = CT_GIV(v(1), v(2));
    G = [G [cos;sin]];
end
CH = CT_H(G);
assert(sum(CH(1,:) == conj(G(1,:)))==ntest,'failed at: CT_H case 5');
assert(sum(CH(2,:) == -G(2,:))==ntest,'failed at: CT_H case 6');
assert(CH(1,1) == conj(G(1,1)),'failed at: CT_H case 7');
assert(CH(2,1) == -G(2,1),'failed at: CT_H case 8');
    
%% test CT_TO_MAT
% Tests for the function to compute the 2x2 matrix representation of a core
% transformation
% case 1: correctness
% REAL
v = randn(2,1); 
[cos, sin, ~] = CT_GIV( v(1), v(2));
M = CT_TO_MAT([cos;sin]);
assert(M(1,1)==cos,'failed at: CT_TO_MAT case 1.1');
assert(M(1,2)==sin,'failed at: CT_TO_MAT case 1.2');
assert(M(2,1)==-conj(sin),'failed at: CT_TO_MAT case 1.3');
assert(M(2,2)==conj(cos),'failed at: CT_TO_MAT case 1.4');
% COMPLEX
v = randn(2,1) + randn(2,1)*1i; 
[cos, sin, ~] = CT_GIV( v(1), v(2));
M = CT_TO_MAT([cos;sin]);
assert(M(1,1)==cos,'failed at: CT_TO_MAT case 2.1');
assert(M(1,2)==sin,'failed at: CT_TO_MAT case 2.2');
assert(M(2,1)==-conj(sin),'failed at: CT_TO_MAT case 2.3');
assert(M(2,2)==conj(cos),'failed at: CT_TO_MAT case 2.4');
% case 2 (combination with CT_GIV): does it introduce a zero?
% REAL
for kk=1:ntest
v = randn(2,1); 
[cos, sin, r] = CT_GIV( v(1), v(2));
M = CT_TO_MAT([cos;sin]);
Mv = M*v;
assert(abs(Mv(1)-r) <=10*eps,'failed at: CT_TO_MAT case 3.1');
assert(abs(Mv(2)) <=2*eps,'failed at: CT_TO_MAT case 3.2');
% COMPLEX
v = randn(2,1) + randn(2,1)*1i; 
[cos, sin, r] = CT_GIV( v(1), v(2));
M = CT_TO_MAT([cos;sin]);
Mv = M*v;
assert(abs(Mv(1)-r) <=10*eps,'failed at: CT_TO_MAT case 4.1');
assert(abs(Mv(2)) <=4*eps,'failed at: CT_TO_MAT case 4.2');
end
%% test CT_FUSE
% Test the fusion of two core transformations
% test 1: correctness
% REAL
v = randn(2,1); 
[cos, sin, ~] = CT_GIV( v(1), v(2));
G1 = [cos; sin];
v = randn(2,1); 
[cos, sin, ~] = CT_GIV( v(1), v(2));
G2 = [cos; sin];
H = CT_FUSE(G1,G2);
assert(H(1) == G1(1)*G2(1)-G1(2)*conj(G2(2)),'failed at: CT_FUSE case 1.1');
assert(H(2) == G1(1)*G2(2)+G1(2)*conj(G2(1)),'failed at: CT_FUSE case 1.2');
% COMPLEX
v = randn(2,1) +randn(2,1)*1i; 
[cos, sin, ~] = CT_GIV( v(1), v(2));
G1 = [cos; sin];
v = randn(2,1) +randn(2,1)*1i; 
[cos, sin, ~] = CT_GIV( v(1), v(2));
G2 = [cos; sin];
H = CT_FUSE(G1,G2);
assert(H(1) == G1(1)*G2(1)-G1(2)*conj(G2(2)),'failed at: CT_FUSE case 2.1');
assert(H(2) == G1(1)*G2(2)+G1(2)*conj(G2(1)),'failed at: CT_FUSE case 2.2');
% REAL & COMPLEX
v = randn(2,1);
[cos, sin, ~] = CT_GIV( v(1), v(2));
G1 = [cos; sin];
v = randn(2,1) +randn(2,1)*1i; 
[cos, sin, ~] = CT_GIV( v(1), v(2));
G2 = [cos; sin];
H = CT_FUSE(G1,G2);
assert(H(1) == G1(1)*G2(1)-G1(2)*conj(G2(2)),'failed at: CT_FUSE case 3.1');
assert(H(2) == G1(1)*G2(2)+G1(2)*conj(G2(1)),'failed at: CT_FUSE case 3.2');
% test 2: correct result
% REAL
for kk=1:ntest
v = randn(2,1); 
[cos, sin, ~] = CT_GIV( v(1), v(2));
G1 = [cos; sin];
v = randn(2,1); 
[cos, sin, ~] = CT_GIV( v(1), v(2));
G2 = [cos; sin];
H = CT_FUSE(G1,G2);
v = randn(2,1);
G1G2v = CT_TO_MAT(G1) * ( CT_TO_MAT(G2) * v);
Hv = CT_TO_MAT(H) *v;
assert(abs(Hv(1)-G1G2v(1))<=5*eps,'failed at: CT_FUSE case 4.1');
assert(abs(Hv(2)-G1G2v(2))<=5*eps,'failed at: CT_FUSE case 4.2');
end
% COMPLEX
for kk=1:ntest
v = randn(2,1) + randn(2,1)*1i; 
[cos, sin, ~] = CT_GIV( v(1), v(2));
G1 = [cos; sin];
v = randn(2,1) + randn(2,1)*1i; 
[cos, sin, ~] = CT_GIV( v(1), v(2));
G2 = [cos; sin];
H = CT_FUSE(G1,G2);
v = randn(2,1);
G1G2v = CT_TO_MAT(G1) * ( CT_TO_MAT(G2) * v);
Hv = CT_TO_MAT(H) *v;
assert(abs(Hv(1)-G1G2v(1))<=10*eps,'failed at: CT_FUSE case 5.1');
assert(abs(Hv(2)-G1G2v(2))<=10*eps,'failed at: CT_FUSE case 5.2');
end
%% test CT_TURNOVER
% test validity of result
% REAL CASE
for kk=1:ntest
M = randn(3,3);
v = randn(2,1); 
[cos, sin, ~] = CT_GIV( v(1), v(2));
G1 = [cos; sin];
v = randn(2,1); 
[cos, sin, ~] = CT_GIV( v(1), v(2));
G2 = [cos; sin];
v = randn(2,1); 
[cos, sin, ~] = CT_GIV( v(1), v(2));
G3 = [cos; sin];
[H1, H2, H3] = CT_TURNOVER(G1,G2,G3);
% from udu to dud
GM = M; 
GM(1:2,1:3) = CT_TO_MAT(G3) * GM(1:2,1:3);
GM(2:3,1:3) = CT_TO_MAT(G2) * GM(2:3,1:3);
GM(1:2,1:3) = CT_TO_MAT(G1) * GM(1:2,1:3);
HM = M;
HM(2:3,1:3) = CT_TO_MAT(H3) * HM(2:3,1:3);
HM(1:2,1:3) = CT_TO_MAT(H2) * HM(1:2,1:3);
HM(2:3,1:3) = CT_TO_MAT(H1) * HM(2:3,1:3);
assert(norm(GM-HM,'fro') < 12*eps,'failed at: CT_TURNOVER case 1.1');
% from dud to udu
HM = M; 
HM(1:2,1:3) = CT_TO_MAT(G3) * HM(1:2,1:3);
HM(2:3,1:3) = CT_TO_MAT(G2) * HM(2:3,1:3);
HM(1:2,1:3) = CT_TO_MAT(G1) * HM(1:2,1:3);
GM = M;
GM(2:3,1:3) = CT_TO_MAT(H3) * GM(2:3,1:3);
GM(1:2,1:3) = CT_TO_MAT(H2) * GM(1:2,1:3);
GM(2:3,1:3) = CT_TO_MAT(H1) * GM(2:3,1:3);
assert(norm(GM-HM,'fro') < 12*eps,'failed at: CT_TURNOVER case 1.2');
end
% COMPLEX CASE
for kk=1:ntest
M = randn(3,3) + randn(3,3)*1i;
v = randn(2,1) + randn(2,1)*1i; 
[cos, sin, ~] = CT_GIV( v(1), v(2));
G1 = [cos; sin];
v = randn(2,1) + randn(2,1)*1i; 
[cos, sin, ~] = CT_GIV( v(1), v(2));
G2 = [cos; sin];
v = randn(2,1) + randn(2,1)*1i; 
[cos, sin, ~] = CT_GIV( v(1), v(2));
G3 = [cos; sin];
[H1, H2, H3] = CT_TURNOVER(G1,G2,G3);
% from udu to dud
GM = M; 
GM(1:2,1:3) = CT_TO_MAT(G3) * GM(1:2,1:3);
GM(2:3,1:3) = CT_TO_MAT(G2) * GM(2:3,1:3);
GM(1:2,1:3) = CT_TO_MAT(G1) * GM(1:2,1:3);
HM = M;
HM(2:3,1:3) = CT_TO_MAT(H3) * HM(2:3,1:3);
HM(1:2,1:3) = CT_TO_MAT(H2) * HM(1:2,1:3);
HM(2:3,1:3) = CT_TO_MAT(H1) * HM(2:3,1:3);
assert(norm(GM-HM,'fro') < 20*eps,'failed at: CT_TURNOVER case 2.1');
% from dud to udu
HM = M; 
HM(1:2,1:3) = CT_TO_MAT(G3) * HM(1:2,1:3);
HM(2:3,1:3) = CT_TO_MAT(G2) * HM(2:3,1:3);
HM(1:2,1:3) = CT_TO_MAT(G1) * HM(1:2,1:3);
GM = M;
GM(2:3,1:3) = CT_TO_MAT(H3) * GM(2:3,1:3);
GM(1:2,1:3) = CT_TO_MAT(H2) * GM(1:2,1:3);
GM(2:3,1:3) = CT_TO_MAT(H1) * GM(2:3,1:3);
assert(norm(GM-HM,'fro') < 20*eps,'failed at: CT_TURNOVER case 2.2');
end