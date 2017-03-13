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
    plus_ops = find(KLidx==1);
    min_ops = find(KLidx==0);
    
    % zero padding arrays to appropriate size
    V = padarray( V , [0 k], 0, 'post');
    KR = padarray( KR , [k k], 0, 'post');
    LR = padarray( LR , [k k], 0, 'post');
    KLrot = padarray( KLrot, [0 k], -1, 'post');
    KLidx = padarray( KLidx, [0 k], -1, 'post');
    
    % main loop
    for i=1:k
        curr_ip1 = start_idx + i;
        curr_i = curr_ip1 - 1;
        
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
        h = zeros(curr_ip1,1);
        for kk=1:2
            for j=1:curr_i
                    hc = V(:,j)'*w;
                    h(j) = h(j) + hc;
                    w = w - hc*V(:,j);
            end
        end
        h(curr_ip1) = norm(w,2);

        % (3) compute rotations and store results
        V(:,curr_ip1) = w/h(curr_ip1);
        if selection(i) > 0, %A*v
%             % old
%             KR(vec,curr_i) = -1;
%             % apply previous rotation
%             for kk = 1:size(KLrot,2)
%                 if ( KLidx(kk) == 0) % K
%                     KR(kk:kk+1,curr_i) = CreateRotMat(RotH(KLrot(:,kk))) * KR(kk:kk+1,curr_i) ;
%                 elseif ( KLidx(kk) == 1) % L
%                     h(kk:kk+1) = CreateRotMat(RotH(KLrot(:,kk))) * h(kk:kk+1);
%                 end
%             end
            
            
            % new -> works for + and - but not yet for consecutive runs
            %L
            for kk = plus_ops 
                h(kk:kk+1) = CreateRotMat(RotH(KLrot(:,kk))) * h(kk:kk+1);
            end
            %K
            if (curr_i>1)
                if isempty(plus_ops) % change occurred
                    ln = 1; % always in this case
                    KR(ln:curr_i,curr_i) = -ACCT(KLrot(:,ln:curr_i-1));
                else
                    if (plus_ops(end) == curr_i-1)  % no change
                        KR(curr_i,curr_i) = -1 ;
                    else % change from A\v to A*v, accumulate last rotations
                        % determine length of series of consecutive A\v ops
                        ln = 1; % default
                        for kk = curr_i-1:-1:1
                            %if selection(kk) > 0
                            if any(plus_ops==kk)
                                ln = kk + 1;
                                break;
                            end
                        end
                        KR(ln:curr_i,curr_i) = -ACCT(KLrot(:,ln:curr_i-1));
                    end
                end
            else
                KR(1,1) = -1 ;
            end

            % compute new rotation
            [c,s,r] = RotGIV(h(curr_i),h(curr_ip1));
            h(curr_i:curr_ip1) = [r,0];
            LR(1:curr_ip1,curr_i) = -h;
            KLrot(:,curr_i) = [conj(c); -s]; 
            KLidx(curr_i) = 1 ; 
            plus_ops(end+1) = curr_i;
            
        else %A\v
%             % old
%             LR(vec,curr_i) = 1;
%             % apply previous rotations
%             for kk = 1:size(KLrot,2)
%                 if ( KLidx(kk) == 1 ) % L
%                     LR(kk:kk+1,curr_i) = CreateRotMat(RotH(KLrot(:,kk))) * LR(kk:kk+1,curr_i) ;
%                 elseif ( KLidx(kk) == 0 ) % K
%                     h(kk:kk+1) = CreateRotMat(RotH(KLrot(:,kk))) * h(kk:kk+1);
%                 end
%             end
            
            % new
            %K
            for kk = min_ops 
                h(kk:kk+1) = CreateRotMat(RotH(KLrot(:,kk))) * h(kk:kk+1);
            end
            %L
            if (curr_i>1)
                if isempty(min_ops), % change occurred
                     lp = 1; % always in this case
                     LR(lp:curr_i,curr_i) = ACCT(KLrot(:,lp:curr_i-1));
                else
                    if (min_ops(end) == curr_ip1-2) % no change
                        LR(curr_i,curr_i) = 1 ;
                    else % change from Av, accumulate last rotations
                        % determine length of series of consecutive A\v ops
                        lp = 1; %default
                        for kk = curr_i-1:-1:1
                            %if selection(kk) < 0 % change to minops !!
                            if any(min_ops==kk)
                                lp = kk + 1;
                                break;
                            end
                        end
                        %ln = find(diff(min_ops)>1,1,'last')
                        LR(lp:curr_i,curr_i) = ACCT(KLrot(:,lp:curr_i-1));
                    end
                end
            else
                LR(1,1) = 1 ;
            end
            
            % compute new rotation
            [c,s,r] = RotGIV(h(curr_i),h(curr_ip1));
            h(curr_i:curr_ip1) = [r,0];
            KR(1:curr_ip1,curr_i) = h;
            KLrot(:,curr_i) = [conj(c); -s]; 
            KLidx(curr_i) = 0 ; 
            min_ops(end+1) = curr_i;
        end
    end
end

function [v] = ACCT(rot)
    % Accumulates the Hermitian conjegates of the core transformation in 
    % rot
    v = zeros(size(rot,2)+1,1);
    sinprod = 1;    % product of sines
    v(1) = conj(rot(1,1));
    for i = 1:size(rot,2)-1
        sinprod = sinprod * conj(rot(2,i));
        v(i+1) = conj(rot(1,i+1)) * sinprod ;
    end
    v(end) = conj(rot(2,end)) * sinprod;
end