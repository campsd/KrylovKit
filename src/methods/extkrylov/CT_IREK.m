function [ rval, res_init, err_init, res_rest, err_rest, nbRest, last_sz ] = CT_IREK( funpos, funneg, v, s, m, p, l, lambda, tol, info, GEP )
% CTIREK - Core-Transformed Implicitly Restarted Extended Krylov algorithm
% for determining rightmost eigenvalues.
%
% IN
% funpos    --  evaluates the MV
% funneg    --  SYSSLV
% v         --  start vector
% s         --  selection vector
% m         --  maximum dimension subspace
% p         --  number of undesired Ritz values that can be filtered out
% l         --  number of rightmost Ritz values that need to have converged
%               to stop the iteration 
% lambda    --  (complex conjugate) pair of target eigenvalues
% info      --  boolean to print information
% GEP       --  boolean that indicates GEP or not
%
% OUT
% rval      --  converged Ritzvalues
% res_init  --  residual during the initial run
% err_init  --  error during the initial run
% res_rest  --  residual during restarts
% err_rest  --  error during restarts
% nbRest    --  number of restarts
% last_sz   --  dimension of final subspace upon convergence
%
% daan.camps@cs.kuleuven.be, December 14, 2016.
    
    if GEP % Generalized eigenvalue problem
        global A LfacA UfacA pfacA B LfacB UfacB pfacB;
        normA = norm(A,inf);
        normB = norm(B,inf);
    else
        %global A Lfac Ufac pfac;
        normA = norm(A,inf);
        normB = 0;
    end
    
    % Deter
    % check if problem makes sense
    n = size(A,1);
    assert(n == length(v), 'incompatible matrix-vector');
    p_init = p;
    
    % initialize variables for building EK subspace
    V = zeros(n,1);
    V(:,1) = v/norm(v,2);
    KLrot = zeros(2,0); KLidx = zeros(1,0);
    KR = zeros(1,0); LR = zeros(1,0);
    
    % initialize cell arrays for tracking convergence
    err_init = cell(2,1); % we track the convergence of two eigenvalues
    res_init = cell(2,1);
    res_rest = cell(2,0); % number of columns will be dependent on nbRest
    err_rest = cell(2,0);
    
    % parameters for flow
    notConverged = true;
    nbRest = 0;
    
    if info
        fprintf('********************************************\n');
        fprintf('Starting CTIREK for problem of size %d\n',n);
        fprintf('********************************************\n');
        fprintf('Parameters:\n');
        fprintf('--------------------------------------------\n');
        fprintf('Tolerance:\t%e\t for %d Ritz values\n',tol,l);
        fprintf('Maximum size m EK subspace:\t%d\n',m);
        fprintf('Restart size p:\t%d\n',p);
        fprintf('--------------------------------------------\n');
        fprintf('\n');
    end
    
    % build initial Krylov subspace
    for i=1:m
        
        % Expand size by one
        sc = s(mod(i,length(s))+1);
        if (i==m), sc = 1; end % We enforce a positive step at the end
        [V,KLrot,KLidx,KR,LR] = CTEK(funpos,funneg,V,KLrot,KLidx,KR,LR,sc);
        
        % If we performed an MV, check convergence
        if sc == 1 && i > l-1 % i larger than 1 because we check for a complex conjugate pair
            % Construct pencil
            [K,L] = CT_EK_PENCIL(KLrot, KLidx, KR, LR);
            % Determine Ritz values and residuals
            [rval, rvec, res] = checkRitzRes(V,K,L,GEP,normA,normB,A,B);
            % Check convergence
            %detShifts = false;
            [res_desired, err_desired, notConverged] = checkConvergence(rval, res, l, lambda, tol);
            %[res_desired, err_desired, notConverged, ~, ~, ~] = checkConvergence(rval, res, lambda, tol, detShifts, p_init);
            res_init{1} = [res_init{1}; [i res_desired(1)]];
            err_init{1} = [err_init{1}; [i err_desired(1)]];
            res_init{2} = [res_init{2}; [i res_desired(2)]];
            err_init{2} = [err_init{2}; [i err_desired(2)]];
            % If converged, we quit
            if ~notConverged
                last_sz = i; % set the final size
                break;
            end
        end
    end    
    
    if info
        fprintf('After initial run:\n')
        fprintf('--------------------------------------------\n');
        fprintf('Min res lambda1\t lambda2\n');
        l1 = res_init{1}; l2 = res_init{2};
        fprintf('%e\t%e\n',min(l1(:,2)),min(l2(:,2)));
        fprintf('Min err lambda1\t lambda2\n');
        l1 = err_init{1}; l2 = err_init{2};
        fprintf('%e\t%e\n',min(l1(:,2)),min(l2(:,2)));
        fprintf('--------------------------------------------\n');
        fprintf('*** MAX RESIDUAL: %e\n', max(res_desired));
        fprintf('--------------------------------------------\n');
        fprintf('\n');
    end
    
    % Restart the iteration and keep doing so until convergence
    while notConverged
        nbRest = nbRest + 1;
        % Increase size of variables
        res_rest{1,nbRest} = [];
        err_rest{1,nbRest} = [];
        % Determine shifts
        %detShifts = true;
        %[~, ~, ~, shifts_s, shifts_cc, p] = checkConvergence(rval, res, lambda, tol, detShifts, p_init);
        [shifts_s, shifts_cc, p] = detShifts(rval, p_init);
        % Apply single shifts
        if ~isempty(shifts_s)
            if info
                fprintf('--------------------------------------------\n');
                fprintf('Removing %d real Ritz values\n',length(shifts_s));
                fprintf('--------------------------------------------\n');
                fprintf('\n');
            end
            [V,KLrot,KLidx,KR,LR] = CTIR(V,KLrot,KLidx,KR,LR,shifts_s);
        end
        % Apply double shifts
        if ~isempty(shifts_cc)
            if info
                fprintf('--------------------------------------------\n');
                fprintf('Removing %d pairs of complex conjugate Ritz values\n',size(shifts_cc,2));
                fprintf('--------------------------------------------\n');
                fprintf('\n');
            end
            [V,KLrot,KLidx,KR,LR] = DSCTIR(V,KLrot,KLidx,KR,LR,shifts_cc);
        end
        
        % Expand subspace
        for i=1:p
             % Expand size by one
            sc = s(mod(i,length(s))+1);
            if (i==p), sc = 1; end % We enforce a positive step at the end
            [V,KLrot,KLidx,KR,LR] = CTEK(funpos,funneg,V,KLrot,KLidx,KR,LR,sc);
            
            % If we performed an MV, check convergence
            if sc == 1
                % Construct pencil
                [K,L] = CT_EK_PENCIL(KLrot, KLidx, KR, LR);
                % Determine Ritz values and residuals
                [rval, rvec, res] = checkRitzRes(V,K,L,GEP,normA,normB,A,B);
                % Check convergence
                %detShifts = false;
                %[res_desired, err_desired, notConverged, ~, ~, ~] = checkConvergence(rval, res, lambda, tol, detShifts, p_init);
                [res_desired, err_desired, notConverged] = checkConvergence(rval, res, l, lambda, tol);
                res_rest{1,nbRest} = [res_rest{1,nbRest}; [i+m-p res_desired(1)]];
                err_rest{1,nbRest} = [err_rest{1,nbRest}; [i+m-p err_desired(1)]];
                res_rest{2,nbRest} = [res_rest{2,nbRest}; [i+m-p res_desired(2)]];
                err_rest{2,nbRest} = [err_rest{2,nbRest}; [i+m-p err_desired(2)]];
                % If converged, we quit
                if ~notConverged
                    last_sz = i; % set the final size
                    break;
                end
            end   
        end
        
        if info
            fprintf('After restart %d:\n',nbRest)
            fprintf('--------------------------------------------\n');
            fprintf('Min res lambda1\t lambda2\n');
            l1 = res_rest{1,nbRest}; l2 = res_rest{2,nbRest};
            fprintf('%e\t%e\n',min(l1(:,2)),min(l2(:,2)));
            fprintf('Min err lambda1\t lambda2\n');
            l1 = err_rest{1,nbRest}; l2 = err_rest{2,nbRest};
            fprintf('%e\t%e\n',min(l1(:,2)),min(l2(:,2)));
            fprintf('--------------------------------------------\n');
            fprintf('*** MAX RESIDUAL: %e\n', max(res_desired));
            fprintf('--------------------------------------------\n');
            fprintf('\n');
        end
    
    end
end

function [rval, rvec, res] = checkRitzRes(V,K,L,GEP, normA, normB,A,B)
    % Checks the Ritz values and relative residuals
    % ! IT IS ASSUMED SC = 1
    [rvec, rval] = eig(L(1:end-1,:), K(1:end-1,:));
    rval = diag(rval);
    fact =  abs(L(end,end));
    %res = fact * abs(rvec(end,:));
    if GEP % Generalized Eigenvalue Problem
        res = [];
        for kk=1:length(rval)
            x = V(:,1:end-1) * K(1:end-1,:)  * rvec(:,kk);
            res(kk) = norm(A*x-rval(kk)*B*x,inf) /(normA + normB + abs(rval(kk)));
        end
    else
        for kk=1:length(res)
            res(kk) = res(kk) / (normA + abs(rval(kk)));
        end 
    end
%     OBSOLETE
%     for kk=1:length(res)
%         res(kk) = res(kk) / norm(V(:,1:end-1) * K(1:end-1,:)  * rvec(:,kk),2);
%     end
end

function  [res_desired, err_desired, notConverged] = checkConvergence(rval, res, l, lambda, tol)
    % Analyzes the Ritz values, determines residuals and error for desired
    % eigenvalues. Checks for convergence.
    [~,ind] = sort(real(rval),'descend');
    res_desired = res(ind(1:l));
    err_desired = [(min(abs(lambda(1)-rval(ind(1:l)))) / abs(lambda(1))); ...
                   (min(abs(lambda(2)-rval(ind(1:l)))) / abs(lambda(2)))];
    
    notConverged = true;
    if max(res_desired) < tol
        notConverged = false;
    end
end

function [shifts_s, shifts_cc, p] =detShifts(rval, p)
    % Determines p exact shifts from the Ritz values split up in complex
    % conjugate pairs and real shifts.
    [~,ind] = sort(real(rval),'descend');
    shifts = rval(ind(end-p+1:end));
    try
        shifts = cplxpair(shifts);
    catch ME
        if (strcmp(ME.identifier,'MATLAB:cplxpair:ComplexValuesPaired'))
            % if this error is thrown, the shifts can't be paired in
            % complex conjugate pairs. This most likely means that only
            % one ritz value of a complex pair is included. We remove
            % one and try again.
            fprintf('Changing p from %d to %d\n',p,p-1);
            p = p - 1;
            shifts = rval(ind(end-p+1:end));
            shifts = cplxpair(shifts);
        end
    end 

    % Sort the shifts in complex conjugate pairs and real shifts
    kk = 1; shifts_cc = []; shifts_s = [];
    while kk <= length(shifts)
        if (abs(imag(shifts(kk))) > 0)
            shifts_cc(:,end+1) = [shifts(kk); shifts(kk+1)];
            kk = kk + 2;
        else
            shifts_s(end+1) = shifts(kk);
            kk = kk + 1;
        end
    end
end
