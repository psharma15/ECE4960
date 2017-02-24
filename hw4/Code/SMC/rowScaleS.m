%% Pragya Sharma, ps847, 20th Feb 2017.
% This code scales i-th row by a and adds it to j-th row of square full
% rank matrices in Row-Compressed Sparse format.
function [valueNew,rowPtrNew,colIndNew] = rowScaleS(value,rowPtr,colInd,i,j,a)
m = length(rowPtr);
diffrow = zeros(1,m);
for iter0 = 1:m-1
    diffrow(iter0+1) = rowPtr(iter0+1)-rowPtr(iter0);    
end
% Getting indicies of non-zero columns in i and j rows and converting from
% sparse to full.
indI = rowPtr(i):rowPtr(i+1)-1;
indJ = rowPtr(j):rowPtr(j+1)-1;
colindI = zeros(1,m-1);
colindJ = zeros(1,m-1);
colindI(colInd(indI)) = 1;
colindJ(colInd(indJ)) = 1;
sum = 0;
% Getting number of final non-zero columns in j-th row after calculation
for iter1 = 1:m-1
    if colindI(iter1)|| colindJ(iter1)
        sum = sum + 1;
    end
end
coljNew = zeros(1,sum);
valuejNew = zeros(1,sum);
counter = 1;
counteri = 1;
counterj = 1;
% Getting new j-th row by looping over the column elements. Do calculation
% by checking possible cases of both columns non-zero or either one
% non-zero, or both zero.
for iter2 = 1:m-1
    if (colindI(iter2) == 1) && (colindJ(iter2) == 1)
        coljNew(counter) = iter2;
        valuejNew(counter) = a*value(indI(counteri))+value(indJ(counterj));
        counteri = counteri + 1;
        counterj = counterj + 1;
        counter = counter + 1;
    else if (colindI(iter2) == 1) && (colindJ(iter2) == 0)
            coljNew(counter) = iter2;
            valuejNew(counter) = a*value(indI(counteri));
            counteri = counteri + 1;
            counter = counter + 1;
        else if (colindI(iter2) == 0) && (colindJ(iter2) == 1)
                coljNew(counter) = iter2;
                valuejNew(counter) = value(indJ(counterj));
                counterj = counterj + 1;
                counter = counter + 1;
            end
        end
        
    end
end
% Concatenating old results with the new j-th row and returning in
% row-compressed format.
valueNew = [value(1:rowPtr(j)-1), valuejNew, value(rowPtr(j+1):rowPtr(end)-1)];
colIndNew = [colInd(1:rowPtr(j)-1), coljNew, colInd(rowPtr(j+1):rowPtr(end)-1)];
diffrow(j+1) = sum;
rowPtrNew = [1 zeros(1,m-1)];
for iter3 = 1:m-1
    rowPtrNew(iter3+1) = rowPtrNew(iter3) + diffrow(iter3+1);    
end
clearvars -except valueNew rowPtrNew colIndNew
end
