function [ rval, conv ] = CTIREK( funpos, funneg, pattern, m, k, l, p )
% CTIREK - Core-Transformed Implicitly Restarted Extended Krylov
% for determining rightmost eigenvalues.
% m = maximum size Krylov subspace
% k = number of desired rightmost Ritz values
% l = buffer to be kept
% p = number of undesired Ritz values
% m = p + l + k
    global A LfacA UfacA pfacA B LfacB UfacB pfacB;
    tol = sqrt(eps);
    tol = 1e-1;
    assert(m== p + l +  k,'incorrect restart specs');
    
    n = size(A,1);
    p_init = p;
    l_init = l;
    %v=ones(n,1); v=v/norm(v);
    v=randn(n,1); v=v/norm(v);
    V = zeros(n,1);
    V(:,1) = v/norm(v,2);
    KLrot = zeros(2,0); KLidx = zeros(1,0);
    KR = zeros(1,0); LR = zeros(1,0);
    
    % build initial Krylov subspace
    for i=1:m-1
        s = pattern(mod(i,length(pattern))+1);
        [V,KLrot,KLidx,KR,LR] = CTEK(funpos,funneg,V,KLrot,KLidx,KR,LR,s);
        
    end
    % mth positive step
    [V,KLrot,KLidx,KR,LR] = CTEK(funpos,funneg,V,KLrot,KLidx,KR,LR,1);
    notConverged = true;
    
    while notConverged
        % compute projected counterpart and Ritz pairs
        [ PC ] = CTPC( KLrot, KLidx, KR, LR, [] );
        [rvec, rval] = eig(PC);
        rval = diag(rval);
        % Residuals
        Lnorm = CreateRotMat(KLrot(:,end))*[LR(end-1,end); 0];
        Lnorm = abs(Lnorm(2));
        res = Lnorm * abs(rvec(end,:));
        % Determine p rightmost Ritzvalues
        [~,ind] = sort(real(rval));
        ind_shifts = ind(1:p);      % p left-most Ritz values
        ind_buffer = ind(p+1:p+l);  % l buffer Ritz values
        ind_desired = ind(p+l+1:m); % k right-most Ritz values

        % Determine maximum residual of p - 2 rightmost Ritzvalues
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
        min(min(abs(KLrot)))
        % Apply single shifts
        if ~isempty(mu_SS)
            [V,KLrot,KLidx,KR,LR] = CTIR(V,KLrot,KLidx,KR,LR,mu_SS);
        end
        min(min(abs(KLrot)))
        % Apply double shifts
        if ~isempty(mu_DS)
            [V,KLrot,KLidx,KR,LR] = DSCTIR(V,KLrot,KLidx,KR,LR,mu_DS);
        end
        min(min(abs(KLrot)))
        % Expand subspace
        for i=1:p-1
            [V,KLrot,KLidx,KR,LR] = CTEK(funpos,funneg,V,KLrot,KLidx,KR,LR,pattern(mod(i,length(pattern))+1));
        end
         % pth positive step
        [V,KLrot,KLidx,KR,LR] = CTEK(funpos,funneg,V,KLrot,KLidx,KR,LR,1);

    end
    rval = rval(ind_desired);
end

