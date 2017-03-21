function [Rot] = CT_TO_MAT(Grot)
    % Creates 2x2 rotational matrix based on the input core transformation
    % which stores the first row of the transformation
    Rot = [Grot(1) Grot(2); -conj(Grot(2)) conj(Grot(1))];
end