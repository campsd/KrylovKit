function [G1,G2,G3]=CT_TURNOVER(H1,H2,H3,pos,tol)
% CT_TURNOVER  Performs the shift-through (turnover) operation, 
% with some additional parameters
%
%    Transitions in both directions are performed. This means:
%
%      /\  -> \  /    OR \  / -> /\
%     /  \     \/         \/    /  \
%
%
%
%    [G1,G2,G3]=ST(H1,H2,H3,pos,tol)
% 
%    All involved rotations are stored in the following format [c,s]^T, 
%    where
%        
%        | c       s       |   
%    G = |                 |   
%        |-conj(s) conj(c) |   
%
%
%    - pos equals 'l' or 'r', meaning that energy is pushed to the left, or
%    to the right. (If H1 and H3 are close to the identity, there is
%    namely flexibility in choosing the position of H2 left or right.)
%    - tol is the tolerance to decide whether a rotation is zero or not
%
%
%    Software written by Raf Vandebril
%    raf.vandebril-at-cs.kuleuven.be
%    Revision Date: September 12, 2011
%
%  Modified by Thomas Mach
%  
%  Reducing the necessary flops
%  partially explicit computations have been used
%
%       SLIMMED DOWN VERSION OF STlr4b
%
% This function performs the shift through-lemma
% It changes slightly the order in which the rotations are computed,
% to retain maximum accuracy.
  G1 = zeros(2,1); G2 = zeros(2,1); G3 = zeros(2,1);
  
  cH12H32=conj(H1(2))*H3(2);
  H31H11=H3(1)*H1(1);
  % explicit computed first rotation
  G1(1)=conj(H1(1)*H3(2))*H2(1)+H3(1)*conj(H1(2)); %p1
  G1(2)=H2(2)*H3(2);
  % explicit computed second rotation
  G2(2)=norm(G1);
  G1=G1./G2(2);
  % explicit compute cosine second rotation
  G2(1)=H31H11-H2(1)*conj(cH12H32);
 
  % Compute the lower right 2x2 block
  B(1,1)=-cH12H32+conj(H31H11)*H2(1);
  B(1,2)=H2(2)*conj(H1(1));
  B(2,1)=-conj(H2(2))*conj(H3(1));
  B(2,2)=conj(H2(1));
  
  % Compute only the relevant parts of G1*B:
  [c,s,th3]=CT_GIV(B(1,2)*conj(G1(2))+G1(1)*B(2,2),B(1,1)*conj(G1(2))+G1(1)*B(2,1));
  % Store G3
  G3(1)=c;
  G3(2)=-s;
end