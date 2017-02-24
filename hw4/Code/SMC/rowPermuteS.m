%% Pragya Sharma, ps847, 20th Feb 2017
% This code interchanges rows i and j, where i<j for Sparse Row-Compressed
% Matrix
function [valueNew,rowPtrNew,colIndNew] = rowPermuteS(value,rowPtr,colInd,i,j)
valueNew = value;
colIndNew = colInd;
m = length(rowPtr);
diffRow = zeros(1,m);
for iter1 = 1:m-1
    diffRow(iter1+1) = rowPtr(iter1+1)-rowPtr(iter1);    
end
rowPtrNew = [1 zeros(1,m-1)];
% Keeping everything initially same, bringing row j to i-th row
counter1 = rowPtr(i);
counter2 = rowPtr(i+1);
for iter2 = 1:diffRow(j+1)
    colIndNew(counter1) = colInd(rowPtr(j)+iter2-1);
    valueNew(counter1) = value(rowPtr(j)+iter2-1);
    counter1 = counter1 + 1;
end
% Before jt-th row, keeping everything else same, adjusting indices
while counter2 < rowPtr(j)
    colIndNew(counter1) = colInd(counter2);
    valueNew(counter1) = value(counter2);
    counter1 = counter1 + 1;
    counter2 = counter2 + 1;
end
% At j-th row, bringing  contents of i-th row, rest everything will be same
for iter3 = 1:diffRow(i+1)
    colIndNew(counter1) = colInd(rowPtr(i)+iter3-1);
    valueNew(counter1) = value(rowPtr(i)+iter3-1);
    counter1 = counter1 + 1;
end

temp = diffRow(i+1);
diffRow(i+1) = diffRow(j+1);
diffRow(j+1) = temp;
for iter4 = 1:m-1
    rowPtrNew(iter4+1) = rowPtrNew(iter4) + diffRow(iter4+1);    
end
clearvars -except valueNew rowPtrNew colIndNew
end

