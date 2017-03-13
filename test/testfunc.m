function [ o ] = testfunc( a )
%TESTFUNC Summary of this function goes here
%   Detailed explanation goes here
    testval = evalin('base','testval');
    o = testval + a ;

end

