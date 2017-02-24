%% Pragya Sharma, ps847, 20th Feb 2017
% This code takes row-compressed sparse matrix and gives full matrix o/p
function [A] = retrieveElement(value,rowPtr,colInd)
n = max(colInd);
A = zeros(n,n);
count = 1;
for i = 1:length(rowPtr)-1
    for j = 1:rowPtr(i+1)-rowPtr(i)
        A(i,colInd(count)) = value(count);
        count = count + 1;
    end
end
end
