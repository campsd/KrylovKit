% ready for 5
%
%  ConsMatLeft  Reconstructs a matrix using rotational patterns
%
%    A=ConsMatLeft(Grot,R) Constructs the matrix A=QR, where the unitary
%    matrix Q is represented by the factorization Grot, Gind.
%
%    A=ConsMatLeft(Grot,R,rowset,colset) constructs A=QR, where only 
%    rotations on the rows in rowset are considered. 
%    Only the cols in colset are computed.
%    A \in \R^{|rowset|,|colset|}
%    R \in \R^{n,m}
%    rowset=1:n resp. colset=1:m if rowset resp. colset unset
%    scale 's' apply scaling of crossing rotations 
%          'n' do not apply scaling, DEFAULT
%
%    R is an upper triangular matrix, on which the rotations will be applied.
%    Grot can appear in two formats.
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


function R=ConsMatLeft(Grot,R,rowset,colset,scale)

[m,n]=size(R);
if (~exist('rowset'))
  rowset=1:m;
end
if (~exist('colset'))
  colset=1:n;
end
if (~exist('scale'))
  scale = 'n';
end



  
% Check if the input is not empty
if (~isempty(Grot))
  if (size(Grot,1)==4)
    [Grot]=DataConvert(Grot);
  end
  
  % Sorting the rotations acting on rowset

  Gselect=[];
  if (scale=='s')
    for k=1:size(Grot,2)
      if (sum(Grot(3,k)==rowset)>=1 || sum(Grot(4,k)==rowset)>=1)
	Gselect=[Gselect, Grot(:,k)];
      end
    end
  else    
    for k=1:size(Grot,2)
      if (sum(Grot(3,k)==rowset)>=1 && sum(Grot(4,k)==rowset)>=1)
	Gselect=[Gselect, Grot(:,k)];
      end
    end
  end
    
  if (~isempty(Gselect))
    % The rotations are re-ordered sot that they can be performed
    % immediately on the matrix R
    [~,I]=sort(Gselect(5,:));
    Gselect=Gselect(:,I); % Reorder the rotations themselves
    
    % Perform all operations on the matrix.
% $$$     if (scale=='s')
% $$$       for i=size(Gselect,2):-1:1
% $$$ 	if (sum(Gselect(3,i)==rowset)>=1 && sum(Gselect(4,i)==rowset)>=1)
% $$$ 	  % Assign the rotation to Giv and execute it.
% $$$ 	  Giv=[Gselect(1,i) Gselect(2,i) ; -conj(Gselect(2,i)) conj(Gselect(1,i))];
% $$$ 	  R([Gselect(3,i),Gselect(4,i)],colset)=Giv*R([Gselect(3,i),Gselect(4,i)],colset);
% $$$ 	else
% $$$ 	  if (sum(Gselect(3,i)==rowset)>=1)
% $$$ 	    R(Gselect(3,i),colset)=Gselect(1,i)*R(Gselect(3,i),colset);
% $$$ 	  else 
% $$$ 	    R(Gselect(4,i),colset)=Gselect(1,i)*R(Gselect(4,i),colset);
% $$$ 	  end
% $$$ 	end
% $$$       end        
% $$$     else
      for i=size(Gselect,2):-1:1
	% Assign the rotation to Giv and execute it.
	Giv=[Gselect(1,i) Gselect(2,i) ; -conj(Gselect(2,i)) conj(Gselect(1,i))];
	R([Gselect(3,i),Gselect(4,i)],colset)=Giv*R([Gselect(3,i),Gselect(4,i)],colset);
      end    
% $$$     end
  end
end
R=R(rowset,colset);


% wrong implementation
% Check if the input is not empty
% $$$ if (~isempty(Grot))
% $$$   if (size(Grot,1)==4)
% $$$     [Grot]=DataConvert(Grot);
% $$$   end
% $$$   
% $$$   % Sorting the rotations acting on rowset
% $$$   Gselect=[];
% $$$   for k=1:size(Grot,2)
% $$$     if (sum(Grot(3,k)==rowset)>=1 && sum(Grot(4,k)==rowset)>=1)
% $$$       Gselect=[Gselect, Grot(:,k)];
% $$$     end
% $$$   end
% $$$ 
% $$$   if (scale=='s') 
% $$$     for k=1:size(Grot,2)
% $$$       if (sum(Grot(3,k)==rowset)>=1 && sum(Grot(4,k)==rowset)<1)
% $$$ 	R(:,Grot(3,k)) = Grot(1,k)*R(:,Grot(3,k));
% $$$       end
% $$$       if (sum(Grot(4,k)==rowset)>=1 && sum(Grot(3,k)==rowset)<1)
% $$$ 	R(:,Grot(4,k)) = conj(Grot(1,k))*R(:,Grot(4,k));
% $$$       end
% $$$     end
% $$$   end  
% $$$ 
% $$$   
% $$$   if (~isempty(Gselect))
% $$$     % The rotations are re-ordered sot that they can be performed
% $$$     % immediately on the matrix R
% $$$     [~,I]=sort(Gselect(5,:));
% $$$     Gselect=Gselect(:,I); % Reorder the rotations themselves
% $$$     
% $$$     % Perform all operations on the matrix.
% $$$     for i=size(Gselect,2):-1:1
% $$$       % Assign the rotation to Giv and execute it.
% $$$       Giv=[Gselect(1,i) Gselect(2,i) ; -conj(Gselect(2,i)) conj(Gselect(1,i))];
% $$$       R([Gselect(3,i),Gselect(4,i)],colset)=Giv*R([Gselect(3,i),Gselect(4,i)],colset);
% $$$     end    
% $$$   end

