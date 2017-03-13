function [ Amv ] = funcneg( v )
% Syssolve with LU factorization of matrix A extracted from base workspace
% (Lfac,Ufac,pfac) with input vector v
    %L = evalin('base','Lfac');
    %U = evalin('base','Ufac');
    %p = evalin('base','pfac');
    global Lfac Ufac pfac
    Amv = mldivide(Ufac,mldivide(Lfac,v(pfac)));
end

