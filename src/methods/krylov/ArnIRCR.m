% ArnIRCR Computes an implicit restart of a compressed represented Arnoldi
% decomposition. 
%
% INPUT
% V     Krylov Basis
% Grot  c,s Givens rotations Arnoldi fact.
% Gind  order Givens rotations Arnoldi fact.
% R     Upper triangular matrix
% mu    shifts
%
% OUTPUT
% V     Adjusted Krylov basis of dim-len(mu)
% Grot  c,s Givens rotations adjusted Arnoldi fact.
% Gind  order Givens rotations adjusted Arnoldi fact.
% R     Upper triangular matrix
%
% We make use of the fact that the compressed format of the matrix in an
% Arnoldi factorisation is of Hessenberg-type. The rotations are chased
% down. 
%
% daan.camps@cs.kuleuven.be
% July 5, 2016

% On Rotations
% -------------------------------------------------------------------------
% The convention we take with respect to Givens rotations is the following:
%  (1)  whenever a rotation is determined by [c,s,~] = Givens(x,y) this
%       means that
%                   |   c           s       |
%           Grot =  |                       |
%                   |   -conj(s)    conj(c) |
%  (2) The similarity transformations we excecute are assumed to be of the 
%       form Htilde = G^H * H * G, so the conjugate transpose rotation is
%       multiplied from the left
%                       | conj(c)       -s  |
%           Grot^H =    |                   |
%                       | conj(s)       c   |

% Main function
% -------------------------------------------------------------------------
function [V,Grot,Gind,R] = ArnoldiCTIR(V,Grot,Gind,R,mu)
    k = size(R,1);
    for kk=1:length(mu)    

        % Determine initial rotation of shifted H
        ChaRot = InitRot(Grot(:,1), R(1:2,1), mu(kk));
        V = ApplyRotToV(RotH(ChaRot), V, 1);
        [WorkRot, R] = ShiftRotThroughR(ChaRot, R, 1, k-1);
	
        % Initial fusion
        Grot(:,1) = Fusion(ChaRot,Grot(:,1));
        ChaRot = RotH(WorkRot);

        % Main loop
        for i=1:k-3
            [ChaRot, Grot(:,i),Grot(:,i+1)] = RotST(Grot(:,i),Grot(:,i+1), ChaRot);
            V = ApplyRotToV(ChaRot, V, i+1);
            [ChaRot, R] = ShiftRotThroughR(RotH(ChaRot), R, i+1, k-1);
            ChaRot = RotH(ChaRot);
        end

        % Last rotation
        i = k-2;
        [ChaRot, Grot(:,i),Grot(:,i+1)] = RotST(Grot(:,i),Grot(:,i+1), ChaRot);
        V = ApplyRotToV(ChaRot, V, i+1);
        
        % Discard last basis vector, rotation and part of R
        V(:,end) = []; 
        Grot(:,end) = []; Gind(end,:) = [];
        R(:,end) = []; R(end,:) = [];
    end
end

% Additional functionality to manipulate rotations
% -------------------------------------------------------------------------
function [Grot] = InitRot(Grot, R, mu)
    % Determines initial rotation R = [R11, R21], Grot = Grot(:,1)
    R = CreateRotMat(Grot)*R - [mu;0];
    [c,s,~]=RotGIV(R(1),R(2));
    Grot=[c, s];  %G(1,:)
end

function [Grot, R] = ShiftRotThroughR(Grot, R, i, k)
    % Shifts a rotation from the right of the upper triangular matrix to
    % the left
    % INPUT
    % Grot  c,s of right Givens rotation 
    % R     Upper Triangular
    % i     rotation acts on cols/rows i&i+1
    % k     problem size
    %
    % OUTPUT
    % Grot  c,s of left Givens rotation
    % R     Upper Triangular
    R(1:i+1,i:i+1) = R(1:i+1,i:i+1)*CreateRotMat(RotH(Grot));
    [c,s,~]=RotGIV(R(i,i),R(i+1,i));
    Grot = [c,s];
    R(i:i+1,i:k) = CreateRotMat(Grot)*R(i:i+1,i:k);
end

function [V] = ApplyRotToV(Grot, V, i)
    % Applies a rotation to the Krylov basis
    V(:,i:i+1) = V(:,i:i+1)*CreateRotMat(Grot);
end


