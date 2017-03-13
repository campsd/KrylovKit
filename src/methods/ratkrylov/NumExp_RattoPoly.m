% Numerical example to illustrate the relation between the rational Krylov
% and the polynomial Krylov
clc; clear all; close all;
format short e;
n = 50;
Ab = diag(-1:3/(n-1):2);
Bb = diag(-1:-1/(n-1):-2);
Xi = {[1 0.21 1 0],[1 0.06 0 1],[1 -0.11 1 0],[1 -0.21 0 1],[1 -0.31 1 1]};
%% CASE 1: A, B NONSINGULAR
A = Ab; B = Bb;
v = ones(n,1); v = v/norm(v);
%RKS
T = eye(6,5);
%T = zeros(6,5); T(1,:) = ones(1,5);
Vrab = zeros(n,1); 
Vrab(:,1) = v; 
K = zeros(1,0); L = zeros(1,0);
[ Vrab, K, L ] = RKS( A, B, @funcpos_generalized, @funcneg_generalized, Vrab, K, L, Xi,T);
%Arnoldi 1
vs1 = v; 
for i=1:length(Xi)
    pole = Xi{i};
    vs1 = ((pole(1)*A+pole(2)*B)\(B)) * vs1;
end
[Vp1ab,H1,flag] = arnoldiRef(B\A,vs1,length(Xi));
%Arnoldi 2
vs2 = v; 
for i=1:length(Xi)
    pole = Xi{i};
    vs2 = ((pole(1)*A+pole(2)*B)\(A)) * vs2;
end
[Vp2ab,H2,flag] = arnoldiRef(A\B,vs2,length(Xi));

%Figure
figure;

for i=1:size(Vrab,2)
   subplot(1, size(Vrab,2), i)
   %axis square
   hold on;
   plot(Vrab(:,i),'Linewidth',2);
   plot(Vp1ab(:,i),'Linewidth',2);
   plot(Vp2ab(:,i),'Linewidth',2);
end

% CHECK SPAN
norm(Vrab - Vp1ab*Vp1ab'*Vrab)
norm(Vp1ab - Vrab*Vrab'*Vp1ab)
norm(Vrab - Vp2ab*Vp2ab'*Vrab)
norm(Vp2ab - Vrab*Vrab'*Vp2ab)
norm(Vp1ab - Vp2ab*Vp2ab'*Vp1ab)
norm(Vp2ab - Vp1ab*Vp1ab'*Vp2ab)
%% CASE 2: A NON SINGULAR, B SINGULAR
A = Ab; B = Bb;
B(10,10) = 0; B(20,20) = 0; B(30,30) = 0; B(40,40) = 0; B(50,50) = 0;
%RKS
Vra = zeros(n,1); 
Vra(:,1) = v; 
K = zeros(1,0); L = zeros(1,0);
[ Vra, K, L ] = RKS( A, B, @funcpos_generalized, @funcneg_generalized, Vra, K, L, Xi,T);
%Arnoldi 2
vs2 = v; 
for i=1:length(Xi)
    pole = Xi{i};
    vs2 = ((pole(1)*A+pole(2)*B)\(A)) * vs2;
end
[Vp2a,H2,flag] = arnoldiRef(A\B,vs2,length(Xi));

%Figure
figure;

for i=1:size(Vra,2)
   subplot(1, size(Vra,2), i)
   %axis square
   hold on;
   plot(Vra(:,i),'Linewidth',2);
   plot(Vp2a(:,i),'Linewidth',2);
end
% CHECK SPAN
norm(Vra - Vp2a*Vp2a'*Vra)
norm(Vp2a - Vra*Vra'*Vp2a)

%% CASE 3: A SINGULAR, B NONSINGULAR
A = Ab; B = Bb;
A(10,10) = 0; A(20,20) = 0; A(30,30) = 0; A(40,40) = 0; A(50,50) = 0;
%RKS
Vrb = zeros(n,1); 
Vrb(:,1) = v; 
K = zeros(1,0); L = zeros(1,0);
[ Vrb, K, L ] = RKS( A, B, @funcpos_generalized, @funcneg_generalized, Vrb, K, L, Xi,T);
%Arnoldi 1
vs1 = v; 
for i=1:length(Xi)
    pole = Xi{i};
    vs1 = ((pole(1)*A+pole(2)*B)\(B)) * vs1;
end
[Vp1b,H1,flag] = arnoldiRef(B\A,vs1,length(Xi));

%Figure
figure;

for i=1:size(Vrb,2)
   subplot(1, size(Vrb,2), i)
   %axis square
   hold on;
   plot(Vrb(:,i),'Linewidth',2);
   plot(Vp1b(:,i),'Linewidth',2);
end
% CHECK SPAN
norm(Vrb - Vp1b*Vp1b'*Vrb)
norm(Vp1b - Vrb*Vrb'*Vp1b)

%% CASE 4: A SINGULAR, B SINGULAR
A = Ab; B = Bb;
A(10,10) = 0; B(20,20) = 0; A(30,30) = 0; B(40,40) = 0; A(50,50) = 0;
%RKS
Vr = zeros(n,1); 
Vr(:,1) = v; 
K = zeros(1,0); L = zeros(1,0);
[ Vr, K, L ] = RKS( A, B, @funcpos_generalized, @funcneg_generalized, Vr, K, L, Xi,T);

%Figure
figure;

for i=1:size(Vr,2)
   subplot(1, size(Vr,2), i)
   %axis square
   hold on;
   plot(Vr(:,i),'Linewidth',2);
end

%%
figure;
for i=1:size(Vrab,2)
   subplot(4, size(Vrab,2), i)
   %axis square
   hold on;
   plot(Vrab(:,i),'Linewidth',2,'Color',[0 0 0]);
   plot(Vp1ab(:,i),'Linewidth',2,'Color',[66/255 134/255 244/255]);
   plot(Vp2ab(:,i),'Linewidth',2,'Color',[244/255 66/255 83/255]);
   ylim([-1 1]);
   xlim([0 50]);
end
for i=1:size(Vra,2)
   subplot(4, size(Vra,2), i+6)
   %axis square
   hold on;
   plot(Vra(:,i),'Linewidth',2,'Color',[0 0 0]);
   plot(Vp2a(:,i),'Linewidth',2,'Color',[244/255 66/255 83/255]);
   ylim([-1 1]);
   xlim([0 50]);
end
for i=1:size(Vrb,2)
   subplot(4, size(Vrb,2), i+12)
   %axis square
   hold on;
   plot(Vrb(:,i),'Linewidth',2,'Color',[0 0 0]);
   plot(Vp1b(:,i),'Linewidth',2,'Color',[66/255 134/255 244/255]);
   ylim([-1 1]);
   xlim([0 50]);
end
for i=1:size(Vr,2)
   subplot(4, size(Vr,2), i+18)
   %axis square
   hold on;
   plot(Vr(:,i),'Linewidth',2,'Color',[0 0 0]);
   ylim([-1 1]);
   xlim([0 50]);
end
% %%
% eigAB = eig(A,B);
% eigH1 = eig(H1(1:5,:));
% eigH2i = 1./eig(H2(1:5,:));
% eigLK = eig(L(1:5,:),K(1:5,:));
% eigVr = eig(Vrab'*(B\A)*Vrab);
% eigVp1 = eig(Vp1ab'*(B\A)*Vp1ab);
% eigVp2 = eig(Vp2ab'*(B\A)*Vp2ab);
% 
% figure;
% plot(real(eigAB),imag(eigAB),'.');
% hold on;
% plot(real(eigLK),imag(eigLK),'o');
% plot(real(eigH1),imag(eigH1),'x');
% plot(real(eigH2i),imag(eigH2i),'s');

% %% OPTION 1
% % Confirm that the span of Vp1 and Vr is the same
% nrmsVp1Vr = zeros(size(Vr,2),1);
% for i=1:size(Vr,2)
%     nrmsVp1Vr(i) = norm(Vr(:,i)-Vp1*Vp1'*Vr(:,i),2);
% end
% % Span Vp2 and Vr ( this is much worse since it is numerically unstable )
% nrmsVp2Vr = zeros(size(Vr,2),1);
% for i=1:size(Vr,2)
%     nrmsVp2Vr(i) = norm(Vr(:,i)-Vp2*Vp2'*Vr(:,i),2);
% end
% % Span Vp1 and Vp2
% nrmsVp1Vp2 = zeros(size(Vr,2),1);
% for i=1:size(Vr,2)
%     nrmsVp1Vp2(i) = norm(Vp1(:,i)-Vp2*Vp2'*Vp1(:,i),2);
% end
% figure;
% semilogy(nrmsVp1Vr,'Linewidth',2);
% hold on
% semilogy(nrmsVp2Vr,'Linewidth',2);
% semilogy(nrmsVp1Vp2,'Linewidth',2);
% ylim([1e-18, 1e-0])
% legend('Vp1 vs Vr','Vp2 vs Vr', 'Vp1 vs Vp2');

