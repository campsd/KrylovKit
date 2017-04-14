function [V,Hrot,HR] = CT_SK_IR_SS(V,Hrot,HR,mu)
% [V,Hrot,HR] = CT_SK_IR_SS(V,Hrot,HR,mu)
% -- Single-shift implicit restart for CT_SK 
%
% INPUT
% V     standard Krylov Basis (V x m+1)
% Hrot  rotations upper Hessenberg (2 x m)
% HR    upper triangular matrix Hessenberg (m+1 x m)
% mu    array with shifts (k)
%
% OUTPUT
% V     Adjusted Krylov basis with len(mu) fewer vectors (V x m-k+1)
% Hrot  core transformations adjusted Hessenberg (2 x m-k)
% HR    Upper triangular matrix adjusted Hessenberg (m-k+1 x m-k)
%
% The rotations are chased down one by one, this is the single shift
% variant
%
% daan.camps@cs.kuleuven.be
% last edit: April 10, 2017
%
% See also: CT_SK
    k = size(HR,1);
    for kk=1:length(mu)    

        % Determine initial rotation of shifted H
        ChaRot = InitRot(Hrot(:,1), HR(1:2,1), mu(kk));
        V = ApplyRotToV(CT_H(ChaRot), V, 1);
        ChaRot = CT_H(ChaRot);
        [WorkRot, HR] = ShiftRotThroughR(ChaRot, HR, 1, k-1);
	
        % Initial fusion
        Hrot(:,1) = CT_FUSE(ChaRot,Hrot(:,1));
        ChaRot = CT_H(WorkRot);

        % Main loop
        for i=1:k-3
            [ChaRot, Hrot(:,i),Hrot(:,i+1)] = CT_TURNOVER(Hrot(:,i),Hrot(:,i+1), ChaRot);
            V = ApplyRotToV(ChaRot, V, i+1);
            [ChaRot, HR] = ShiftRotThroughR(ChaRot, HR, i+1, k-1);
            ChaRot = CT_H(ChaRot);
        end

        % Last rotation
        i = k-2;
        [ChaRot, Hrot(:,i),Hrot(:,i+1)] = CT_TURNOVER(Hrot(:,i),Hrot(:,i+1), ChaRot);
        V = ApplyRotToV(ChaRot, V, i+1);
        
        % Discard last basis vector, rotation and part of R
        V(:,end) = []; 
        Hrot(:,end) = [];
        HR(:,end) = []; HR(end,:) = [];
    end
    
    % Nested functionS (shared workspace)
    function ShiftRotThroughR(i)
        % Shifts a rotation from the right of the upper triangular matrix to
        % the left.
        % INPUT
        % i     rotation acts on cols/rows i&i+1
        HR(1:i+1,i:i+1) = HR(1:i+1,i:i+1)*CT_TO_MAT(ChaRot);
        [c,s,~]=RotGIV(HR(i,i),HR(i+1,i));
        ChaRot = [c;s];
        HR(i:i+1,i:k-1) = CT_TO_MAT(ChaRot)*HR(i:i+1,i:k-1);
    end

    function ApplyRotToV(Grot, i)
        % Applies a rotation to the Krylov basis
        V(:,i:i+1) = V(:,i:i+1)*CT_TO_MAT(Grot);
    end
end

% Additional functionality to manipulate rotations
% -------------------------------------------------------------------------
function [Grot] = InitRot(Grot, R, mu)
    % Determines initial rotation R = [R11, R21], Grot = Grot(:,1)
    R = CT_TO_MAT(Grot)*R - [mu;0];
    [c,s,~]=RotGIV(R(1),R(2));
    Grot=[c; s];  %G(1,:)
end
