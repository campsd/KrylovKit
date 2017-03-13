% CTEK - Core-Transformations factorised Extended Krylov
%
% INPUT
% funpos	function that excecutes the matvec product
% funneg	function that solves the system
% V      	Krylov Basis to start from
% KLrot  	n-1 c,s Givens rotations K and L Hessenbergs
% KLidx  	side of Givens rotations 0 = K, 1 = L
% KR     	Upper triangular matrix K Hessenberg
% LR     	Upper triangular matrix L Hessenberg
% selection	defines the number and type of operations to perform
%
%
% OUTPUT
% V     Adjusted Krylov basis of dim-len(mu)
% KLrot  c,s Givens rotations adjusted K Hessenberg
% KLidx  order Givens rotations adjusted K Hessenberg
% KR     Upper triangular matrix adjusted K Hessenberg
% LR     Upper triangular matrix adjusted L Hessenberg
%
% The recursion relation that holds before and after the restart is:
%       A*V*K = V*L
%
%
% daan.camps@cs.kuleuven.be
% September 2, 2016

function [V,KLrot,KLidx,KR,LR] = CTEK(funpos,funneg,V,KLrot,KLidx,KR,LR,selection)
    % selection(i) > 0 : we use funpos
    % selection(i) < 0 : we use funneg
    	
    k = length(selection);
    start_idx = size(V,2);
    
    % zero padding arrays to appropriate size
    V = padarray( V , [0 k], 0, 'post');
    KR = padarray( KR , [k k], 0, 'post');
    LR = padarray( LR , [k k], 0, 'post');
    
    % main loop
    for i=1:k
        curr_idx = start_idx + i;
        % (1) operate
        if selection(i) > 0, %A*v
            vec = find(KLidx==1,1,'last') + 1;
            if isempty(vec), vec = 1; end
            w = funpos(V(:,vec));
            %w = A*V(:,lp);
        else
            vec = find(KLidx==0,1,'last') + 1;
            if isempty(vec), vec = 1;end
            w = funneg(V(:,vec));
            %v = V(:,ln);
            %y = mldivide(Lf,v(p));
            %w = mldivide(Uf,y);%A\V(:,ln);
        end

        % (2) orthogonalise
        h = zeros(curr_idx,1);
        for kk=1:2
            for j=1:curr_idx-1
                    hc = V(:,j)'*w;
                    h(j) = h(j) + hc;
                    w = w - hc*V(:,j);
            end
        end
        h(curr_idx) = norm(w,2);

        % (3) compute rotations and store results
        V(:,curr_idx) = w/h(curr_idx);
        if selection(i) > 0, %A*v
            KR(vec,curr_idx-1) = -1;
            % apply previous rotation
            for kk = 1:size(KLrot,2)
                if ( KLidx(kk) == 0) % K
                    KR(kk:kk+1,curr_idx-1) = CreateRotMat(RotH(KLrot(:,kk))) * KR(kk:kk+1,curr_idx-1) ;
                else % L
                    h(kk:kk+1) = CreateRotMat(RotH(KLrot(:,kk))) * h(kk:kk+1);
                end
            end

            % compute new rotation
            [c,s,r] = RotGIV(h(curr_idx-1),h(curr_idx));
            h(curr_idx-1:curr_idx) = [r,0];
            LR(1:curr_idx,curr_idx-1) = -h;
            KLrot(:,curr_idx-1) = [conj(c); -s]; 
            KLidx(curr_idx-1) = 1 ; 
        else %A\v
            LR(vec,curr_idx-1) = 1;
            % apply previous rotations
            for kk = 1:size(KLrot,2)
                if ( KLidx(kk) == 1 ) % L
                    LR(kk:kk+1,curr_idx-1) = CreateRotMat(RotH(KLrot(:,kk))) * LR(kk:kk+1,curr_idx-1) ;
                else % K
                    h(kk:kk+1) = CreateRotMat(RotH(KLrot(:,kk))) * h(kk:kk+1);
                end
            end

            % compute new rotation
            [c,s,r] = RotGIV(h(curr_idx-1),h(curr_idx));
            h(curr_idx-1:curr_idx) = [r,0];
            KR(1:curr_idx,curr_idx-1) = h;
            KLrot(:,curr_idx-1) = [conj(c); -s]; 
            KLidx(curr_idx-1) = 0 ; 

        end
    end
end

function [v] = ACCT(rot)
    % Accumulates the Hermitian conjegates of the core transformation in 
    % rot
    v = zeros(size(rot,2)+1,1);
    sinprod = 1;    % product of sines
    v(1) = rot(1,1);
    for i = 1:size(rot,2)-1
        sinprod = sinprod * (-conj(rot(2,i)));
        v(i+1) = rot(1,i+1) * sinprod ;
    end
    v(end) = sinprod;
end