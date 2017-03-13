function [ Amv ] = funcpos_generalized( v )
% Returns B^-1*A*v for the generalized eigenvalue problem
    global LfacB UfacB pfacB
    A = evalin('base','A');
    Amv = A*v; 
    Amv = mldivide(UfacB,mldivide(LfacB,Amv(pfacB)));
end
