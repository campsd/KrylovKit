function G=CT_FUSE(G1,G2)
% CT_FUSE   Combines two core transformations to obtain a single core 
%          transformation
%
%    G=CT_FUSE(G1,G2)
% 
%    G1,G2 represent two rotations.
%    G represents the product G1*G2
%
%    All three Givens transformations are of the form [c,s]^T, 
%    where
%        
%        | c       s       |   
%    G = |                 |   
%        |-conj(s) conj(c) |   
%
%    Software written by Raf Vandebril
%    raf.vandebril-at-cs.kuleuven.be
%    Revision Date: September 1, 2010
%
% COPY OF Fusion.m
% Based on two Givens transformations G1 and G2, we fuse these two
% matrices together    
  G(1,1)=G1(1)*G2(1)-G1(2)*conj(G2(2));
  G(2,1)=G1(1)*G2(2)+G1(2)*conj(G2(1));
end
