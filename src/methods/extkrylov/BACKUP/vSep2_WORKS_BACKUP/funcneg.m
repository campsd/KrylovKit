function [ Amv ] = funcneg( v )
% Syssolve with LU factorization of matrix A extracted from base workspace
% (Lfac,Ufac,pfac) with input vector v
    L = evalin('base','Lfac');
    U = evalin('base','Ufac');
    p = evalin('base','pfac');
    y = mldivide(L,v(p));
    Amv = mldivide(U,y);
end

