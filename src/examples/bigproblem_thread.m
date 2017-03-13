load('thread.mat');
A = Problem.A;
display('Matrix loaded, starting eigenvalues computation');
tic;
eigA = eig(full(A));
toc;
% write to file
display('writing to file')
fileID = fopen('eigA_thread.dat','w');
fprintf(fileID','%23.16e %23.16e\n',[real(eigA) imag(eigA)]');
fclose(fileID);