shifts = randn(4,1);
shifts(end+1) = randn + 1i * randn;
shifts(end+1) = conj(shifts(end));
shifts(end+1) = randn + 1i * randn;

for k=1:10
    k
    if k<5
        k = k + 1;
    end
    k
end


kk = 1; mu_DS = []; mu_SS = [];
        while kk < length(shifts)
            if (imag(shifts(kk)) > 0)
                mu_DS(:,end+1) = [shifts(kk); shifts(kk+1)];
                kk = kk + 2;
            else
                mu_SS(end+1) = shifts(kk);
                kk = kk + 1;
            end
        end

%%
v=ones(2*n,1); v=v/norm(v);
V = zeros(2*n,1);
V(:,1) = v/norm(v,2);
KLrot = zeros(2,0); KLidx = zeros(1,0);
KR = zeros(1,0); LR = zeros(1,0);
recerrEK = [];
for i=1:200
    [V,KLrot,KLidx,KR,LR] = CTEK(@funcpos,@funcneg,V,KLrot,KLidx,KR,LR,pattern(mod(i,length(pattern))+1));
    [K,L] = CONS_CTEK_PENCIL(KLrot,KLidx, KR,LR);
    % valid recurrence?
    %norm(A*V*K - B*V*L,'fro') %cavity
    recerrEK(i) = norm(A*V*K - V*L,'fro') %olmstead
    % orthogonal V?
    norm(V'*V)-1
end       



v=ones(2*n,1); v=v/norm(v);
V = zeros(2*n,1);
V(:,1) = v/norm(v,2);
H = zeros(0,0);
r = V(:,1);
recerrorA = [];
for i=1:200
    %Arnoldi v1
    [U,H,flag] = arnoldiRef(A,v,i);
    recerrorA(i) = norm(A*U(:,1:end-1) - U*H,'fro') %olmstead
    % Arnoldi v2
    %[H,V,r]=arnoldi_sorensen(A,H,V,r,i,1)
    % Arnoldi v3
%     [V,H,f] = Arnoldi(A,i,v);
%     ei = eye(i); ei = ei(end,:)
%     recerrorA(i) = norm(A*V - V*H - f*ei,'fro') %olmstead
end
figure;
semilogy(recerrEK);
hold on
semilogy(recerrorA)
%%
[K,L] = CONS_CTEK_PENCIL(KLrot,KLidx, KR,LR);
% valid recurrence?
norm(A*V*K - B*V*L,'fro') %cavity
%norm(A*V*K - V*L,'fro') %olmstead
%%
[ PC ] = CTPC( KLrot, KLidx, KR, LR, [] );
[eigPC] = eig(PC);
figure;
plot(eigAB(:,1),eigAB(:,2),'.');
hold on;
plot(real(eigPC),imag(eigPC),'o');
xlim([-2.5 0.5])
ylim([-5 5])