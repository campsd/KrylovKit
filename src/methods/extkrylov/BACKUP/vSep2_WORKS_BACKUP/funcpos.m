function [ Av ] = funcpos( v )
% Matvec of input vector v with matrix A from 'base' workspace
    A = evalin('base','A');
    Av = A*v;
end

