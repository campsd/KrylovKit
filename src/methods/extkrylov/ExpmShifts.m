function [ shifts ] = ExpmShifts( PC )
%EXPMSHIFTS Shifts for the matrix exponential
RV = eig(PC);
frac = 0.2;
[~,I] = sort(real(RV));
shifts = RV(I(1:ceil(frac*length(RV))));
end

