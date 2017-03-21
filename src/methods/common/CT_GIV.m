function [c,s,r] = CT_GIV(x,y)
% CT_GIV -- compute a core transformation based on a 2-vector
%
% GivensR Compute a Givens transformation matrix
%
%    [c,s,r] = GivensR(x,y) returns the cosine and sine of the complex
%    Givens rotation matrix and the resulting r which is positive and
%    real in this case.
%    
%        | c       s       |                  | x |     | r |
%    G = |                 |   such that  G * |   |  =  |   |
%        |-conj(s) conj(c) |                  | y |     | 0 |
% 
%    where c, s are complex, and c^2 + |s|^2 = 1. 
%
%
%    Software written by Raf Vandebril
%    raf.vandebril-at-cs.kuleuven.be
%    Revision Date: August 30, 2010
%
%
%  Modified by Thomas Mach
%  
%  Reducing the necessary flops
%  partially explicit computations have been used
%  r positive real
if (y == 0)
  c = 1; s = 0; r = x;
else
  if (abs(x) >= abs(y))
    theta=conj(sign(x));
    t = y/x; r = sqrt(1 + abs(t)^2);
    c = theta/r;
    s = conj(t)*c;
    r = theta*x*r;
  else
    theta=conj(sign(y));
    t = x/y; r = sqrt(1 + abs(t)^2 );
    s = theta/r;
    c = conj(t)*s;
    r = theta*y*r;
  end
end
end