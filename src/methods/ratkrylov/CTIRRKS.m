function [ rval, rres, theta ] = CTIRRKS( A, B, v, Xi, m, p, l, funselect, lambda, tol, info )
%CTIRRKS Core-Transformed Implicitly Restarted Rational Krylov Algorithm
%
% IN
% A             --  first matrix of pencil
% B             --  second matrix of pencil
% v             --  start vector (normalized)
% Xi            --  poleset (cell array) !! TODO: function for pole selection
% m             --  maximum dimension subspace
% p             --  number of undesired Ritz values that can be filtered out
% l             --  number of converged Ritz values to keep
% funselect     --  function that orders the Ritz values to criterion
% (DESIRED to UNDESIRED)
% lambda        --  eigenvalues of the pencil
% tol           --  convergence tolerance
% info          --  boolean to print information
%
% OUT
% 
% daan.camps@cs.kuleuven.be

global Lc Uc pc;

assert(sum(size(A) == size(B))==2, 'incompatible pencil');
assert(size(A,2) == size(v,1), 'incompatible matrix-vector');

normA = norm(A,inf);
normB = norm(B,inf);

n = size(A,1);
p_init = p;

% Initialize variables for rational Krylov subspace
% Basis vectors
V = zeros(n,1);
V(:,1) = v/norm(v,2);
% Recurrence matrices
Lrot = zeros(2,0);
KR = zeros(1,0); LR = zeros(1,0);
% Matrix factorizations
ID = []; % ID of the factorization is theta(1)
LU = struct('L',[],'U',[],'p',[]);
% Algorithm paramaters
notConverged = true ;
nbRest = 0;
% Convergence tracking
rres = {}; theta = {};

if info
    fprintf('********************************************\n');
    fprintf('Starting CTIRRKS for problem of size %d\n',n);
    fprintf('********************************************\n');
    fprintf('Parameters:\n');
    fprintf('--------------------------------------------\n');
    fprintf('Tolerance:\t%e\t for %d Ritz values\n',tol,l);
    fprintf('Maximum size m EK subspace:\t%d\n',m);
    fprintf('Restart size p:\t%d\n',p);
    fprintf('--------------------------------------------\n');
    fprintf('\n');
end

% Build initial Krylov subspace
for i=1:m
    
    % Expand the subspace by one
    Xic = Xi{mod(i,length(Xi))+1};
    [Xiconv, thetac, facc] = convPole(Xic);
    % Check if we already computed the factorization for the pole
    IDc = find(ID==thetac(1));
    if isempty(IDc) % the requested factorization doesn't exist
        [Lc,Uc,pc] = lu(Xiconv(1)*A + Xiconv(2)*B,'vector');
        LU.L = cat(3,LU.L,Lc);
        LU.U = cat(3,LU.U,Uc);
        LU.p = cat(1,LU.p,pc);
        ID = cat(1,ID,thetac(1));
    else % it does exist, we get it from memory
        Lc = LU.L(:,:,IDc);
        Uc = LU.U(:,:,IDc);
        pc = LU.p(IDc,:);
    end
    T = zeros(i+1,1); T(i) = 1;
    [V,KR,Lrot,LR] = CTRKS1(A,B,@funcPole,facc,Xic,V,KR,Lrot,LR,T); %TODO implement and set T, XiPos
    % Construct RKS pencil
    L = LR;
    for kk=size(Lrot,2):-1:1
        L(kk:kk+1,:) = CreateRotMat(Lrot(:,kk)) * L(kk:kk+1,:);
    end
    norm(A*V*KR-B*V*L,'fro')
    % Compute Ritzvalues and vectors
    [rvec,rval] = eig(L(1:end-1,:),KR(1:end-1,:));
    rval = diag(rval);
    theta{i} = rval;
    rres{i} = RitzResiduals(lambda,rval);
    % Order the Ritz values according to selection criterion
    [ordrval,ordidx] = funselect(rval); 
    % Compute residuals based on Ritz vectors
    res = zeros(size(rval));
    for kk=1:length(rval)
        x = V(:,1:end-1) * KR(1:end-1,:)  * rvec(:,kk);
        res(kk) = norm(A*x-rval(kk)*B*x,inf) /(normA + normB + abs(rval(kk)));
    end
    if i>=l
        ordres = res(ordidx);
        % Check convergence
        notConverged = checkConvergence(ordres,l,tol);
        % If converged, we quit
        if ~notConverged
            last_sz = i; % set the final size
            break;
        end
    end
end

% Restart the iteration and keep doing so until convergence
while notConverged
    nbRest = nbRest + 1;
    % Determine shifts
    [shifts_s, shifts_cc, p] = detShifts(ordrval, p_init);
    % Apply single shifts
    if ~isempty(shifts_s)
        if info
            fprintf('--------------------------------------------\n');
            fprintf('Removing %d real Ritz values\n',length(shifts_s));
            fprintf('--------------------------------------------\n');
            fprintf('\n');
        end
        [V,KR,Lrot,LR] = CTIR(V,KR,Lrot,LR,shifts_s);
    end
    % Apply double shifts
    if ~isempty(shifts_cc)
        if info
            fprintf('--------------------------------------------\n');
            fprintf('Removing %d pairs of complex conjugate Ritz values\n',size(shifts_cc,2));
            fprintf('--------------------------------------------\n');
            fprintf('\n');
        end
        [V,KR,Lrot,LR] = DSCTIR(V,KR,Lrot,LR,shifts_cc); %TODO implement
    end
    
    % Expand subspace again
    for i=1:p
        % Expand the subspace by one
        Xic = Xi{mod(i,length(Xi))+1};
        [Xiconv, thetac, facc] = convPole(Xic);
        % Check if we already computed the factorization for the pole
        IDc = find(ID==thetac(1));
        if isempty(IDc) % the requested factorization doesn't exist
            [Lc,Uc,pc] = lu(Xiconv(1)*A - Xiconv(2)*B,'vector');
            LU.L = cat(3,LU.L,Lc);
            LU.U = cat(3,LU.U,Uc);
            LU.p = cat(1,LU.p,pc);
            ID = cat(1,ID,thetac(1));
        else % it does exist, we get it from memory
            Lc = LU.L(:,:,IDc);
            Uc = LU.U(:,:,IDc);
            pc = LU.p(IDc,:);
        end
        T = zeros(m-p+i+1,1); T(m-p+i) = 1;
        [V,KR,Lrot,LR] = CTRKS1(A,B,@funcPole,facc,Xic,V,KR,Lrot,LR,T); %TODO implement
        % Construct RKS pencil
        L = LR;
        for kk=size(Lrot,2):-1:1
            L(kk:kk+1,:) = CreateRotMat(Lrot(:,kk)) * L(kk:kk+1,:);
        end
        norm(A*V*KR-B*V*L,'fro')
        % Compute Ritzvalues and vectors
        [rvec,rval] = eig(L(1:end-1,:),KR(1:end-1,:));
        rval = diag(rval);
        theta{end+1} = rval;
        rres{end+1} = RitzResiduals(lambda,rval);
        % Order the Ritz values according to selection criterion
        [ordrval,ordidx] = funselect(rval); 
        % Compute residuals based on Ritz vectors
        res = zeros(size(rval));
        for kk=1:length(rval)
            x = V(:,1:end-1) * KR(1:end-1,:)  * rvec(:,kk);
            res(kk) = norm(A*x-rval(kk)*B*x,inf) /(normA + normB + abs(rval(kk)));
        end
        ordres = res(ordidx);
        % Check convergence
        notConverged = checkConvergence(ordres,l,tol);
        % If converged, we quit
        if ~notConverged
            last_sz = i; % set the final size
            break;
        end
    end
end


end

function [Xiconv,theta,fac] = convPole(Xi)
    % convert the pole Xi to a "pole" on the unit circle
    % Xiconv -- the resulting pole
    % theta -- vector with two angles on the unit circle (Pole ID)
    % fac -- vector with scaling factors
    norm1 = norm(Xi(1:2),2);
    norm2 = norm(Xi(3:4),2);
    Xiconv = [Xi(1:2)./norm1 Xi(3:4)./norm2];
    theta = [angle(Xiconv(1) + 1i*Xiconv(2)) angle(Xiconv(3) + 1i*Xiconv(4))];
    fac = [norm1 norm2];
end

function [w]  = funcPole(poleFac, w)
    % Solves the system with the factorization in the global variables
    global Uc Lc pc;
    w = poleFac * mldivide(Uc,mldivide(Lc,w(pc)));
end

function  [notConverged] = checkConvergence(ordres, l, tol)
    res_desired = ordres(1:l);    
    notConverged = true;
    if max(res_desired) < tol
        notConverged = false;
    end
end

function [shifts_s, shifts_cc, p] =detShifts(ordrval, p)
    % Determines p exact shifts from the Ritz values split up in complex
    % conjugate pairs and real shifts.
    shifts = ordrval(end-p+1:end);
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
            shifts = ordrval(end-p+1:end);
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
