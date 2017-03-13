function [ ordrval, ordidx ] = selectRightmost( rval )
%SELECTRIGHTMOST - orders the rvals such that rightmost come first,
%leftmost last
    [ordrval,ordidx] = sort(real(rval),'descend');

end

