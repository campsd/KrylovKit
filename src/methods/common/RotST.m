% ready for 5
%
% STlr   Performs the shift-through operation, with some additional parameters
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
%
%
%  Modified by Thomas Mach
%  
%  Reducing the necessary flops
%  partially explicit computations have been used
%
%  maple code:
%  with(LinearAlgebra);
%  G1:=Matrix(3,3,[[1,0,0],[0,c,s],[0,-conjugate(s),conjugate(c)]]);
%  H1:=Matrix(3,3,[[d,t,0],[-conjugate(t),conjugate(d),0],[0,0,1]]);
%  J1:=Matrix(3,3,[[1,0,0],[0,e,u],[0,-conjugate(u),conjugate(e)]]); 
%  F:=G1.H1.J1;
%  V2:=SubMatrix(F,[1,2,3],[3]);
%  J3:=GivensRotationMatrix(V2,2,1);
%  J3[2,1]:=conjugate(J3[2,1]); 
%  J3[2,2]:=conjugate(J3[2,2]);
%  F2:=J3.F;
%  V3:=SubMatrix(F2,[1,2,3],[3]);
%  J4:=GivensRotationMatrix(V3,3,2);
%  J4[3,2]:=conjugate(J4[3,2]);
%  J4[3,3]:=conjugate(J4[3,3]); 
%  F3:=J4.F2;
%  V4:=SubMatrix(F3,[1,2,3],[2]);
%  J5:=GivensRotationMatrix(V4,2,1);
%  J5[2,1]:=conjugate(J5[2,1]);
%  J5[2,2]:=conjugate(J5[2,2]);
%  F4:=J5.F3;
%%
%    G1      ,              G2                      ,               G3
%
% conj(p1)                                            H1(2)H2(1)conj(H3(1)) + H1(1)H3(2)
% --------   ,   H3(1)H1(1) - conj(H3(2))H2(1)H1(2) , ----------------------------------
%    q2                                                             q2
%
% H2(2)H3(2)                                                 H1(2)H2(2)
% ---------- ,              q2                      ,        ----------
%    q2                                                          q2
%
%
%  p1 = H1(1)conj(H2(1))H3(2) + H1(2)conj(H2(2))
%  q2 = sqrt( |p1|^2 + |H2(2)|^2|H3(2)|^2 )
%
%
%
%%%
%% COPY OF STlr4b
function [G1,G2,G3]=RotST(H1,H2,H3,pos,tol)
% This function performs the shift through-lemma
% It changes slightly the order in which the rotations are computed,
% to retain maximum accuracy.
    
% It appears that almost identity rotations can create problems, hence
% we set them explicitly equal to the identity
% we take sqrt(3) as 3 is the size of the matrix
% if ((abs(H1(1))-1)<sqrt(3)*eps) H1(2)=0; end;
% if ((abs(H2(1))-1)<sqrt(3)*eps) H2(2)=0; end;
% if ((abs(H3(1))-1)<sqrt(3)*eps) H3(2)=0; end;


tol=1e-400;
pos='l';
% testing is too expensive
% if (~exist('tol','var'))
%   tol=eps;
% end
% if (~exist('pos','var'))
%   pos='l';
% end

if ((size(H1,1)>2)&&(size(H2,1)>2)&&(size(H3,1)>2))
  % [   [       [          [       [   [
  % [ [ [  => [ [ [ resp [ [ [  => [ [ [
  %   [       [   [      [   [       [
  G1=zeros(5,1);
  G2=G1;
  G3=G1;
  G1(3:4)=H2(3:4);
  G2(3:4)=H1(3:4);
  G3(3:4)=H2(3:4);
  G1(5)=H1(5);
  G2(5)=H2(5);
  G3(5)=H3(5);
end

% unchanged
% if ((abs(H1(2))<tol) && (abs(H3(2))<tol))
%   if (pos=='l')
%     % Move the energy to the left side that is G1
%     % implying that the identities are pushed in G3
%     % disp('left')
%     
%     % Push the energy of H2 to G1
%     G1(1)=H2(1);
%     G1(2)=H2(2)*conj(H1(1));
%     
%     % Update G2 and G3
%     G2(1)=H1(1)*H3(1);
%     G2(2)=0;
%     G3(1)=1;
%     G3(2)=0;
%        
%   elseif (pos=='r')
%     % Move the energy to the right side that is G3
%     % implying that the identities are pushed in G1
%     
%     %disp('right')
%     % Push the energy of H2 to G3
%     G3(1)=H2(1);
%     G3(2)=H2(2)*H3(1);
%       
%     % Update G2 and G1
%     G2(1)=H1(1)*H3(1);
%     G2(2)=0;
%     G1(1)=1;
%     G1(2)=0;
%   else
%     error ('invalid choice for pos in ST.m')
%   end
% else % Turnover
  
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
  % B=[conj(G1(1)), -G1(2); conj(G1(2)), G1(1)]*B;
  % [c,s,th3]=GivensR(B(2,2),B(2,1));
  [c,s,th3]=GivensR2(B(1,2)*conj(G1(2))+G1(1)*B(2,2),B(1,1)*conj(G1(2))+G1(1)*B(2,1));
  % Store G3
  G3(1)=c;
  G3(2)=-s;
%end

end
%%% Local Variables: 
%%% flyspell-mode:nil
%%% mode:flyspell-prog
%%% ispell-local-dictionary: "american"
%%% End: 
