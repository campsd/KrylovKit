function [res] = RitzResiduals(lambda, theta)
    % function that computes the true norm of the residuals of the Ritz
    % values 'theta' compared to the true eigenvalues 'lambda'.
    % daan.camps@cs.kuleuven.be
    
    res = [];
    for i=1:length(theta)
       diff = abs(lambda-theta(i));
       res(i) = min(diff);
    end
end