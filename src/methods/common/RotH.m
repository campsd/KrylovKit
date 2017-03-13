function [Grot] = RotH(Grot)
    % Converts a Givens rotation in Grot format to its conjugate transpose
    % s.t. CreateRotMat(Grot)*CreateRotMat(RotH(Grot)) = eye(2)
    Grot = [conj(Grot(1)); -Grot(2)];
end