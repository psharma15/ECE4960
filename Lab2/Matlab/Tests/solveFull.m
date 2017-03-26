% Finds x in Ax = b, for a full matrix
% Pragya Sharma, March 2017

function [x] = solveFull(value,rowPtr,colInd,b)
    A = retrieveElement(value,rowPtr,colInd);
    x = A\b';
end