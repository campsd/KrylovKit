% ready for 5
%
%  ConsMatRight  Reconstructs a matrix using rotational patterns
%
%    A=ConsMatRight(Grot,Gind,R) Constructs the matrix A=RQ, where the unitary
%    matrix Q is represented by the factorization Grot, Gind.
%
%    A=ConsMatRight(Grot,R,colset,rowset) constructs A=RQ, where only 
%    rotations on the cols in colset are considered. 
%    Only the rows in rowset are computed.
%    A \in \R^{|rowset|,|colset|}
%    R \in \R^{n,m}
%    rowset=1:n resp. colset=1:m if rowset resp. colset unset
%    scale 's' apply scaling of crossing rotations 
%          'n' do not apply scaling, DEFAULT
%
%    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%    !! order of rowset and colset interchanged !!
%    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
%    R is an upper triangular matrix, on which the rotations will be applied.
% 4Row data format: 
%    Grot stores a sequence of rotations.
%    It is a 4 x m vector, each column representing an individual rotation.
%    The first 2 components store the cosine and the sine. 
%    For instance Grot(1:2,i)=[c;s] corresponds to
%    G_i=|       c            s    |
%        |    -conj(s)     conj(c) |
%    The third component j identifies rows j,j+1 affected by the rotation.
%    The fourth component k specifies the position in the sequence,
%    going from left to right.
%  
% 
%
% 5Row data format: (allows rotations altering nonsuccessive rows)
%    Grot stores a sequence of rotations.
%    It is a 5 x m vector, each column representing an individual rotation.
%    The first 2 components store the cosine and the sine. 
%    For instance Grot(1:2,i)=[c;s] corresponds to
%    G_i=|       c            s    |
%        |    -conj(s)     conj(c) |
%    The third and fourth component j1,j2 identify the rows j1,j2 affected by the rotation.
%    The fifth component k specifies the position in the sequence,
%    going from left to right.
%
%
%    Software written by Raf Vandebril
%    raf.vandebril-at-cs.kuleuven.be
%    Revision Date: July 5, 2011
%

function R=ConsMatRight(Grot,R,rowset,colset,scale)

[m,n]=size(R);
if (~exist('colset'))
  colset=1:n;
end
if (~exist('rowset'))
  rowset=1:m;
end
if (~exist('scale'))
  scale = 'n';
end



  
  
% Check if the input is not empty
if ~isempty(Grot)
  if (size(Grot,1)==4)
    % For plotting we use the 5 row dataformat
    [Grot]=DataConvert(Grot);
  end
  
  % Sorting the rotations acting on colset
  Gselect=[];
  for k=1:size(Grot,2)
    if (sum(Grot(3,k)==colset)>=1 && sum(Grot(4,k)==colset)>=1)
      Gselect=[Gselect, Grot(:,k)];
    end
  end

  if (scale=='s') 
    for k=1:size(Grot,2)
      if (sum(Grot(3,k)==colset)>=1 && sum(Grot(4,k)==colset)<1)
	R(:,Grot(3,k)) = Grot(1,k)*R(:,Grot(3,k));
      end
      if (sum(Grot(4,k)==colset)>=1 && sum(Grot(3,k)==colset)<1)
	R(:,Grot(4,k)) = conj(Grot(1,k))*R(:,Grot(4,k));
      end
    end
  end  
  
  if (~isempty(Gselect))
    % By default we assume that the matrix Gind is not ordered,
    % Most likely this matrix will be ordered, but for consistency 
    % we first order the rotations in the order that they can be performed
    [SG,I]=sort(Gselect(5,:));
    Gselect=Gselect(:,I); % Rorder the rotations themselves
    
    for i=1:size(Gselect,2)
      % Assign the rotation and execute it
      Giv=[Gselect(1,i) Gselect(2,i) ; -Gselect(2,i)' Gselect(1,i)'];
      R(rowset,[Gselect(3,i),Gselect(4,i)])=R(rowset,[Gselect(3,i),Gselect(4,i)])*Giv;
    end
  end
end    
%R
%colset
%rowset
R=R(rowset,colset);