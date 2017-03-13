function [ rot, ord, R, D ] = ChaseSSpD( rot, ord, R, D, mu)
% -------------------------------------------------------------------------
%CHASE_SSPD -- Chasing algorithm for a semi-separable plus diagonal matrix.
%
% Executes a single shift step on a matrix of condensed semi-separable plus
% diagonal form. This structure arises in a (special) case of the rational
% Krylov method. The chasing procedure can be used for the implicit restart
% of the rational Krylov method.
%
% INPUT
% rot   n-1 rotations stored by c,s
% ord   mutual ordering of rotations: 0 = ,1 = 
% R     upper triangular matrix of size nxn
% D     n diagonal entries (poles)
% mu    shift to apply
%
% OUTPUT
% rot
% ord
% R
% D
%
% daan.camps@cs.kuleuven.be
% January 2, 2017
% -------------------------------------------------------------------------

% Assumptions
%
% * For descending parts in the rotational pattern, we assume the
% corresponding diagonal entries are equal to the next ascending rotation
% (next finite pole)
%
% * ....

descval = 1; % value of ord(..) in descending case
ascval = 0; % value of ord(..) in ascending case

% -------------------------------------------------------------------------
% Determine initial perturbing core transformation
% -------------------------------------------------------------------------
if (ord(1) == descval) % First rotations are descending
    x = CreateRotMat(rot(:,1)) * R(1:2,1) + [D(1); 0] - [mu; 0];
else % First rotations are ascending
    x = mu*CreateRotMat(RotH(rot(:,1)))*[1/R(1,1); 0];
    x = [1;0] - x;
end
[c,s,~] = RotGIV(R(1),R(2));
ChaRot = [c;s];

% -------------------------------------------------------------------------
% Step 0
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Main loop
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Final step
% -------------------------------------------------------------------------
end

function [OutputArgs] = UpdateDiagonal( InputArgs )
    % Updates the diagonal and upper triangular such that it is ready for a
    % similarity transformation
    
end

function [Grot, R] = ShiftRotLeftToRight(Grot, R, i, k)
	% Shifts a rotation from left to right through the upper triangular
	R(i:i+1,i:k) = CreateRotMat(Grot)*R(i:i+1,i:k);
	[c,s,~]=RotGIV(R(i+1,i+1),R(i+1,i));
	Grot = [c,s];
	R(1:i+1,i:i+1) = R(1:i+1,i:i+1)*CreateRotMat(Grot);
end

function [Grot, R] = ShiftRotRightToLeft(Grot, R, i, k)
	% Shifts a rotation from right to left through the upper triangular
	R(1:i+1,i:i+1) = R(1:i+1,i:i+1)*CreateRotMat(Grot);
	[c,s,~]=RotGIV(R(i,i),R(i+1,i));
	Grot = [c,s];
	R(i:i+1,i:k) = CreateRotMat(Grot)*R(i:i+1,i:k);
end
