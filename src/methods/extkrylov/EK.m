% EK -  Extended Krylov
%
% INPUT
% funpos	function that excecutes the matvec product
% funneg	function that solves the system
% V      	Krylov Basis to start from
% K     	Upper Hessenberg
% LR    	Upper Hessenberg
% selection	defines the number and type of operations to perform
%
%
% OUTPUT
% V     Krylov basis
% K     Upper K Hessenberg
% L     Upper triangular matrix adjusted L Hessenberg
%
% The recursion relation that holds:
%       A*V*K = V*L
%
%
% daan.camps@cs.kuleuven.be
% September 12, 2016

function [ V,K,L ] = EK( funpos,funneg,V,K,L,selection )

    k = length(selection);
    start_idx = size(V,2);
    
    % zero padding arrays to appropriate size
    V = padarray( V , [0 k], 0, 'post');
    K = padarray( K , [k k], 0, 'post');
    L = padarray( L , [k k], 0, 'post');
    
    % main loop
    for i=1:k
        curr_idx = start_idx + i;
        % (1) operate
        if selection(i) > 0, %A*v
            vec = find(diag(L,-1)<0,1,'last') + 1;
            if isempty(vec), vec = 1; end
            w = funpos(V(:,vec));
            %w = A*V(:,lp);
        else
            vec = find(diag(K,-1)>0,1,'last') + 1;
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
            K(vec,curr_idx-1) = -1;
            L(1:curr_idx,curr_idx-1) = -h;

        else %A\v
            L(vec,curr_idx-1) = 1;
            K(1:curr_idx,curr_idx-1) = h; 

        end
    end

end

