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
% See also: CT_SK, CT_SK_IR_DS

% TODO: detect when purging is required
% mu = 0 shift fails
    ns = length(mu);
    for kk=1:ns
%         ChaRots = [];
        k = size(HR,1);
        % Determine initial rotation of shifted H
        ChaRot = InitRot(Hrot(:,1), HR(1:2,1), mu(kk));
%         ChaRots = [ChaRots ChaRot];
        ApplyRotToV(CT_H(ChaRot), 1);
	
        % Initial fusion
        Hrot(:,1) = CT_FUSE(ChaRot,Hrot(:,1));
        
        ChaRot = CT_H(ChaRot);
        ShiftRotThroughR(1);
        ChaRot = CT_H(ChaRot);

        % Main loop
        for i=1:k-3
            [ChaRot, Hrot(:,i),Hrot(:,i+1)] = CT_TURNOVER(Hrot(:,i),Hrot(:,i+1), ChaRot);
            ApplyRotToV(ChaRot, i+1);
            ShiftRotThroughR(i+1);
            ChaRot = CT_H(ChaRot);
%             ChaRots = [ChaRots ChaRot];
            
%             figure(10)
%             clf
%             semilogy(abs(Hrot(1,:)))
%             hold on
%             semilogy(abs(Hrot(2,:)))
%             legend('cos','sin')
%             
%             figure(11)
%             clf
%             semilogy(abs(ChaRots(1,:)))
%             hold on
%             semilogy(abs(ChaRots(2,:)))
%             legend('cos','sin')
        end

        % Last rotation
        i = k-2;
        [ChaRot, Hrot(:,i),Hrot(:,i+1)] = CT_TURNOVER(Hrot(:,i),Hrot(:,i+1), ChaRot);
        ApplyRotToV(ChaRot, i+1);

%         ChaRots = [ChaRots ChaRot];
%         figure(10)
%         clf
%         semilogy(abs(Hrot(1,:)))
%         hold on
%         semilogy(abs(Hrot(2,:)))
%         legend('cos','sin')
% 
%         figure(11)
%         clf
%         semilogy(abs(ChaRots(1,:)))
%         hold on
%         semilogy(abs(ChaRots(2,:)))
%         legend('cos','sin')
            
%         % Discard last basis vector, rotation and part of R
        V(:,end) = []; 
        Hrot(:,end) = [];
        HR(:,end) = []; HR(end,:) = [];
    end
    
%     % Discard last basis vector, rotation and part of R
%     V(:,end-ns+1:end) = []; 
%     Hrot(:,end-ns+1:end) = [];
%     HR(:,end-ns+1:end) = []; HR(end-ns+1:end,:) = [];
    
    % Nested functions (shared workspace)
    function ShiftRotThroughR(i)
        % Shifts a rotation from the right of the upper triangular matrix to
        % the left.
        % INPUT
        % i     rotation acts on cols/rows i&i+1
        HR(1:i+1,i:i+1) = HR(1:i+1,i:i+1)*CT_TO_MAT(ChaRot);
        [c,s,HR(i,i),]=CT_GIV(HR(i,i),HR(i+1,i));
        ChaRot = [c;s];
        HR(i+1,i) = 0;
        HR(i:i+1,i+1:end) = CT_TO_MAT(ChaRot)*HR(i:i+1,i+1:end);
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
    [c,s,~]=CT_GIV(R(1),R(2));
    Grot=[c; s];  %G(1,:)
end
