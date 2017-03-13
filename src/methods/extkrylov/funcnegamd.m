function [ Amv ] = funcnegamd( v )
% Syssolve with amd ordered LU factorization of matrix A extracted from base workspace
% (Lfac,Ufac,pfac, Pamd) with input vector v
    %L = evalin('base','Lfac');
    %U = evalin('base','Ufac');
    %p = evalin('base','pfac');
    global Lfacamd Ufacamd pfacamd Pamd
    Amv(Pamd) = mldivide(Ufacamd,mldivide(Lfacamd,v(Pamd(pfacamd))));
    Amv = transpose(Amv);
end

