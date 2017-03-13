function [ rval, conv_init, conv_restart, ls ] = CTIREK( funpos, funneg, pattern, m, k, l, p )
% CTIREK - Core-Transformed Implicitly Restarted Extended Krylov
% for determining rightmost eigenvalues.
% m = maximum size Krylov subspace
% k = number of desired rightmost Ritz values
% l = buffer to be kept
% p = number of undesired Ritz values
% m = p + l + k
    global A LfacA UfacA pfacA B LfacB UfacB pfacB;
    tol = sqrt(eps);
    %tol = 1e-1;
    assert(m== p + l +  k,'incorrect restart specs');
    
    n = size(A,1);
    p_init = p;
    l_init = l;
    v=ones(n,1); v=v/norm(v);
    %v=randn(n,1); v=v/norm(v);
    V = zeros(n,1);
    V(:,1) = v/norm(v,2);
    KLrot = zeros(2,0); KLidx = zeros(1,0);
    KR = zeros(1,0); LR = zeros(1,0);
    
    % convergence and Ritz values
    conv_init = zeros(k, m-k+1);
    conv_restart = zeros(k,p,0);
    notConverged = true;
    
    % build initial Krylov subspace
    for i=1:m
        % Operate
        s = pattern(mod(i,length(pattern))+1);
        [V,KLrot,KLidx,KR,LR] = CTEK(funpos,funneg,V,KLrot,KLidx,KR,LR,s);
        % Construct pencil & determine residuals
        [K,L] = CONS_CTEK_PENCIL(KLrot, KLidx, KR, LR);
        if s == 1, %positive step
            [rvec, rval] = eig(L(1:end-1,:), K(1:end-1,:));
            rval = diag(rval);
            fact =  CreateRotMat(KLrot(:,end))*[LR(end-1,end); 0];
            fact = abs(fact(2))
            res = fact * abs(rvec(end,:));
            for kk=1:length(res)
                res(kk) = res(kk) / norm(V(:,1:end-1) * K(1:end-1,:)  * rvec(:,kk),2);
            end
        else % negative step
            %[rvec, rval] = eig(L(1:end-1,:), K(1:end-1,:)); % option 1 via
            % (L,K) -> large residual
            %[rvec, rval] = eig(V(:,1:end-1)'*A*V(:,1:end-1)); % option 2
            %via V'AV (expensive)
            %option 3 via Ltilde
            fact =  CreateRotMat(KLrot(:,end))*[KR(end-1,end); 0];
            Ltilde = L(1:end-1,:);  
            Ltilde(:,end) = Ltilde(:,end) - fact(2) * V(:,1:end-1)' * (A * V(:,end));
            [rvec, rval] = eig(Ltilde, K(1:end-1,:));
            rval = diag(rval); %rval = 1./rval;
            fact = abs(fact(2))
            %res = fact * norm(A * V(:,end)) * abs(rvec(end,:));
            res = fact * abs(rvec(end,:));
            for kk=1:length(res)
                res(kk) = res(kk) / norm(V(:,1:end-1) * K(1:end-1,:)  * rvec(:,kk),2);
                %res(kk) = res(kk) / norm(V(:,1:end-1) * rvec(:,kk),2);
            end
        end
        % Check convergence
        if i >= k
            [~,ind] = sort(real(rval));
            ind_desired = ind(i-k+1:i); % k rightmost Ritz values
            conv_init(:,i) = res(ind_desired)';
            maxres = max(res(ind_desired))

            if maxres < tol,
                notConverged = false;
                break;
            end
        end    
        
    end
    % mth positive step
    %[V,KLrot,KLidx,KR,LR] = CTEK(funpos,funneg,V,KLrot,KLidx,KR,LR,1);
    nbRestarts = 1;
    
    while notConverged
%         compute projected counterpart and Ritz pairs
%         [ PC ] = CTPC( KLrot, KLidx, KR, LR, [] );
%         [rvec, rval] = eig(PC);
%         rval = diag(rval);
%         Residuals
%         Lnorm = CreateRotMat(KLrot(:,end))*[LR(end-1,end); 0];
%         Lnorm = abs(Lnorm(2));
%         res = Lnorm * abs(rvec(end,:));

        % Split up Ritz values in desired, buffer and shifts
        [~,ind] = sort(real(rval));
        ind_shifts = ind(1:p);      % p left-most Ritz values
        ind_buffer = ind(p+1:p+l);  % l buffer Ritz values
        ind_desired = ind(p+l+1:m); % k right-most Ritz values

        % Determine maximum residual of p rightmost Ritzvalues
        maxres = max(res(ind_desired))

        if maxres < tol,
            notConverged = false;
            break;
        end

        % Excecute restart
        % Check if shifts can be paired
        shifts = rval(ind_shifts);
        try
            shifts = cplxpair(shifts);
        catch ME
            if (strcmp(ME.identifier,'MATLAB:cplxpair:ComplexValuesPaired'))
                disp('changing the p,l');
                p = p_init - 1;
                l = l_init + 1;
                ind_buffer = [ind_shifts(end); ind_buffer];
                ind_shifts(end) = [];
                shifts = rval(ind_shifts);
                shifts = cplxpair(shifts);
            
            end
            %rethrow(ME)
        end 
        
        % Get complex conjugate pairs
        kk = 1; mu_DS = []; mu_SS = [];
        while kk <= length(shifts)
            if (abs(imag(shifts(kk))) > 0)
                mu_DS(:,end+1) = [shifts(kk); shifts(kk+1)];
                kk = kk + 2;
            else
                mu_SS(end+1) = shifts(kk);
                kk = kk + 1;
            end
        end
        %min(min(abs(KLrot)))
        % Apply single shifts
        if ~isempty(mu_SS)
            [V,KLrot,KLidx,KR,LR] = CTIR(V,KLrot,KLidx,KR,LR,mu_SS);
        end
        %min(min(abs(KLrot)))
        % Apply double shifts
        if ~isempty(mu_DS)
            [V,KLrot,KLidx,KR,LR] = DSCTIR(V,KLrot,KLidx,KR,LR,mu_DS);
        end
        %min(min(abs(KLrot)))
        
        % Expand subspace
        for i=1:p
            s = pattern(mod(i,length(pattern))+1);
            [V,KLrot,KLidx,KR,LR] = CTEK(funpos,funneg,V,KLrot,KLidx,KR,LR,s);
            [K,L] = CONS_CTEK_PENCIL(KLrot, KLidx, KR, LR);
            if s == 1, %positive step
                [rvec, rval] = eig(L(1:end-1,:), K(1:end-1,:));
                rval = diag(rval);
                fact =  CreateRotMat(KLrot(:,end))*[LR(end-1,end); 0];
                fact = abs(fact(2));
                res = fact * abs(rvec(end,:));
                for kk=1:length(res)
                    res(kk) = res(kk) / norm(V(:,1:end-1) * K(1:end-1,:)  * rvec(:,kk),2);
                end
            else % negative step
                %[rvec, rval] = eig(L(1:end-1,:), K(1:end-1,:)); %option 1
                %-> large residual
                %[rvec, rval] = eig(V(:,1:end-1)'*A*V(:,1:end-1)); %option
                %2 -> expensive
                %option 3 via Ltilde
                fact =  CreateRotMat(KLrot(:,end))*[KR(end-1,end); 0];
                Ltilde = L(1:end-1,:);  
                Ltilde(:,end) = Ltilde(:,end) - fact(2) * V(:,1:end-1)' * (A * V(:,end));
                [rvec, rval] = eig(Ltilde, K(1:end-1,:));
                rval = diag(rval); %rval = 1./rval;
                %fact =  CreateRotMat(KLrot(:,end))*[KR(end-1,end); 0];
                fact = abs(fact(2));
                res = fact * abs(rvec(end,:));
                %res = fact * norm(A * V(:,end)) * abs(rvec(end,:));
                for kk=1:length(res)
                    res(kk) = res(kk) / norm(V(:,1:end-1) * K(1:end-1,:)  * rvec(:,kk),2);
                    %res(kk) = res(kk) / norm(V(:,1:end-1) * rvec(:,kk),2);
                end
            end
            % Check convergence
            [~,ind] = sort(real(rval));
            ind_desired = ind(end-k+1:end); % k rightmost Ritz values
            conv_restart(:,i,nbRestarts) = res(ind_desired)';
            maxres = max(res(ind_desired))

            if maxres < tol,
                notConverged = false;
                break;
            end
            
        end
        nbRestarts = nbRestarts + 1;
         % pth positive step
        %[V,KLrot,KLidx,KR,LR] = CTEK(funpos,funneg,V,KLrot,KLidx,KR,LR,1);

    end
    ls = i;
    rval = rval(ind_desired);
end

