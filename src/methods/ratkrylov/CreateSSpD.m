function [ R ] = CreateSSpD( rot, ord, R, D )
%CREATESSPD Creates a dense matrix from SSpD format
%   January 10, 2017

descval = 1; % value of ord(..) in descending case
ascval = 0; % value of ord(..) in ascending case

i =1;
while i<=length(ord)
    if ord(i) == ascval
        R(i:i+1,:) = CreateRotMat(rot(:,i))*R(i:i+1,:);
        i = i + 1;
    else
        I = find(ord(i:end)==ascval, 1, 'first'); % find the end of the descending sequence
        if isempty(I) %descending till the end
            I = length(ord);
        else
            I = I - 1;
        end
        for j=I:-1:i
            R(j:j+1,:) = CreateRotMat(rot(:,j))*R(j:j+1,:);
        end
        i = I + 1;
    end
end
R = R + D;

end