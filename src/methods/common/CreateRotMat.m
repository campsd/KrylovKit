function [Rot] = CreateRotMat(Grot)
    % Creates 2x2 rotational matrix based on Grot which stores the first
    % row
    Rot = [Grot(1) Grot(2); -conj(Grot(2)) conj(Grot(1))];
end