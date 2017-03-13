I = eye(size(K,2)); ekt = I(end,:);
pert = V(:,1:end-1)'*A*V(:,end)*ekt*K(end,end);

figure, colormap('summer')
imagesc(log10(abs(pert)))
colorbar, set(gca,'CLim',[-15,0]); axis square

Lpert = L(1:end-1,:) - pert ;
eigA = eig(full(A));
eigLK = eig(L(1:end-1,:),K(1:end-1,:));
eigLpK = eig(Lpert,K(1:end-1,:));
eigLpKm1 = eig(Lpert/K(1:end-1,:));

figure;
plot(eigA,'.');
hold on;
plot(eigLK,'x');
plot(eigLpK,'o');
plot(eigLpKm1,'d');
legend('A','(L,K)','(Lp,K)','Lp/K');