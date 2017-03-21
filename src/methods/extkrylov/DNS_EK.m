function [ V,K,L ] = DNS_EK( funpos,funneg,V,K,L,selection )
% DNS_EK -  (Dense) Extended Krylov without factorization
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
        if selection(i) > 0 %A*v
            vec = find(diag(L,-1)<0,1,'last') + 1;
            if isempty(vec), vec = 1; end
            w = funpos(V(:,vec));
        else
            vec = find(diag(K,-1)>0,1,'last') + 1;
            if isempty(vec), vec = 1;end
            w = funneg(V(:,vec));
        end

        % (2) orthogonalise
        h = zeros(curr_idx,1);
        for kk=1:2 % two passes to ensure orthogonal vectors
            hc = V(:,1:curr_idx-1)'*w;
            h(1:curr_i) = h(1:curr_idx-1) + hc;
            w = w - V(:,1:curr_idx-1)*hc;
        end
        h(curr_idx) = norm(w,2);

        % (3) Store results
        V(:,curr_idx) = w/h(curr_idx);
        if selection(i) > 0 %A*v
            K(vec,curr_idx-1) = -1;
            L(1:curr_idx,curr_idx-1) = -h;

        else %A\v
            L(vec,curr_idx-1) = 1;
            K(1:curr_idx,curr_idx-1) = h; 
        end
    end
end

