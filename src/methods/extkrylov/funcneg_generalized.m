function [ Amv ] = funcneg_generalized( v )
% Returns A^-1*B*v for the generalized eigenvalue problem
    global LfacA UfacA pfacA
    B = evalin('base','B');
    Amv = B*v; 
    Amv = mldivide(UfacA,mldivide(LfacA,Amv(pfacA)));
end
