% Test for ShiftBranchThroughPencil
% This test is for the case we have the AscBranch as a descending branch at
% K side, bring it through KR to the right and remove it from the right so
% that it shows up at the right of LR where we bring it through the upper
% triangular to the left

% Initial 
AscPattInit = [AscBranch; [1,2]; [1,2]];
KRinit=ConsMatLeft(AscPattInit,KR);
LRinit=LR;

% Intermediate
DescPattInt = [TBranch; [1,2]; [1,2]];
KRint=ConsMatRight(DescPattInt,KR);
norm(KRint-KRinit,'fro') %Correct

%KRint=ConsMatLeft(AscPattInit(:,1),KR);
%KRint=ConsMatRight(DescPattInt(:,2),KRint);

%Final
AscPattFin = [AscBranch; [1,2]; [2,1]];
Q = ConsMatLeft(AscPattFin,eye(size(LR,1)));
LRfinleft = ConsMatLeft(AscPattFin,LR);
LRfinleft = Q*LR;
Patt = [[RotH(TBranch(:,1)), RotH(TBranch(:,2))];[1,2];[2,1]];
Q = ConsMatLeft(Patt,eye(size(LR,2)));
LRfinright = ConsMatRight(Patt,LRinit);
LRfinright = LRinit*Q;
LRfinleft-LRfinright;
norm(LRfinleft-LRfinright,'fro') % error is quite big...

%% Test branch switch
%before
Patt = [[DescBranch, AscBranch];[1,2,3,1,2];[2,3,4,2,3];[1,2,3,5,4]];
Gplot(Patt(1:2,:),Patt(3:5,:)',4)
Q = ConsMatLeft(Patt,eye(4));
%after
Patt = [[AscBranch, DescBranch];[2,3,1,2,3];[3,4,2,3,4];[2,1,3,4,5]];
Gplot(Patt(1:2,:),Patt(3:5,:)',4)
Q2 = ConsMatLeft(Patt,eye(4));

%% Test recurrence
upto = 3;
[K,L] = CONS_CTEK_PENCIL(KLrot(:,1:upto-1),KLidx(1:upto-1),KR(1:upto,1:upto-1),LR(1:upto,1:upto-1));
norm(A*V(:,1:upto)*K-V(:,1:upto)*L,'fro')


%% Test augmented recurrence
% At line 48
[K,L] =  CONS_CTEK_PENCIL(KLrot,KLidx,KR,LR);
Patt = [[AscBranch];[1,2];[2,3];[1,2]];
Q = ConsMatLeft(Patt,eye(size(K,1)));
norm(A*V*Q*K-V*Q*L,'fro')
norm(A*V*K-V*L,'fro')

%%
tmpBranch(:,1) = RotH(AscBranch(:,1));
tmpBranch(:,2) = RotH(AscBranch(:,2));
Patt = [[AscBranch, tmpBranch];[1,2,1,2];[2,3,2,3];[1,2,4,3]];
Q = ConsMatLeft(Patt,eye(size(KR,1)));