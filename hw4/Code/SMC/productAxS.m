%% Pragya Sharma, ps847, 20th Feb 2017
% This function takes row compressed matrix input, and carries out product
% with a vector
function [b] = productAxS(value,rowPtr,colInd,x)
n = max(colInd);
count = 1;
b = zeros(n,1);
for i = 1:length(rowPtr)-1
    for j = 1:rowPtr(i+1)-rowPtr(i)
        b(i) = b(i) + value(count)*x(colInd(count));
        count = count + 1;
    end
end
clearvars -except b
end
