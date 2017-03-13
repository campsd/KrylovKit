function [ ordrval, ordidx ] = selectLeftmost( rval )
%SELECTLEFTMOST Summary of this function goes here
    [ordrval,ordidx] = sort(real(rval),'ascend');
end

