function [ V, K, L ] = RKS( A, B, funpos, funneg, V, K, L, Xi, T )
%RKS - rational Krylov sequence alogorithm
%   Detailed explanation goes here

% TODO there is a problem in this implementation that leads to spurious
% eigenvalues, compare with rat_krylov
% non optimal continuation vector? no
% normalize start vector
% ANSWER: the problem appears to be if the matrix B appears in the
% numerator! The gamma and delta parameters appear to influence the result
    k = length(Xi);
    start_idx = size(V,2);
    
    % zero padding arrays to appropriate size
    V = padarray( V , [0 k], 0, 'post');
    K = padarray( K , [k k], 0, 'post');
    L = padarray( L , [k k], 0, 'post');
    
    % main loop
    for i=1:k
        curr_ip1 = start_idx + i;
        curr_i = curr_ip1 - 1;
        
        % (1) Execute the operation
        if Xi{i} == Inf %AB^-1v
            w = funpos(V(:,1:curr_i)*T(1:curr_i,i));
        elseif Xi{i} == 0 %BA^-1v
            w = funneg(V(:,1:curr_i)*T(1:curr_i,i));
        else %[alpha, beta, gamma, delta]
            pole = Xi{i};
            w = (pole(1)*A + pole(2)*B) \ ((pole(3)*A + pole(4)*B)*(V(:,1:curr_i)*T(1:curr_i,i)));
        end
        
        % (2) Orthogonalise the new vector
        h = zeros(curr_ip1,1);
        for kk=1:2
            for j=1:curr_i
                    hc = V(:,j)'*w;
                    h(j) = h(j) + hc;
                    w = w - hc*V(:,j);
            end
        end
        h(curr_ip1) = norm(w,2); %TODO check for breakdown
        
        % (3) Update the recursion
        V(:,curr_ip1) = w/h(curr_ip1);
        if Xi{i} == Inf %AB-1v
            % In this case there only appears a rotation at the L side 
            % Vectors
            ki = -T(1:curr_ip1,i);
            li = -h;
        elseif Xi{i} == 0 %BA-1v
            % In this case there only appears a rotation at the K side
            % Vectors
            ki = h;
            li = T(1:curr_ip1,i);
        else %[alpha, beta, gamma, delta]
            % In this case there appear rotations at both sides
            % Vectors
            ki = pole(1)*h - pole(3)*T(1:curr_ip1,i);
            li = -pole(2)*h + pole(4)*T(1:curr_ip1,i);
        end
            K(1:curr_ip1,curr_i) = ki;
            L(1:curr_ip1,curr_i) = li;
    end

end

