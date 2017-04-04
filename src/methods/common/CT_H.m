function [Grot] = CT_H(Grot)
    % Converts a core transformation to its conjugate transpose
    % s.t. CT_TO_MAT(Grot)*CT_TO_MAT(CT_H(Grot)) = eye(2)
    Grot = [conj(Grot(1,:)); -Grot(2,:)];
end