%% ECE 4960, HW-4
% Submitted by: Pragya Sharma, ps847

function [normTestA,normTestb] = ps847_hw4_p1()
path(path,genpath(pwd));

%% Full and Sparse Matrix operations
% Test Data
value = 1:12;
rowPtr = [0,3,6,9,10,12]+1;
colInd = [0,1,4,0,1,2,1,2,4,3,0,4]+1;
x = sort(1:5,'Descend')';
Afull = [1 2 0 0 3; 
         4 5 6 0 0;
         0 7 8 0 9;
         0 0 0 10 0;
         11 0 0 0 12];
% Operation 1: Row Permutation
[valueNew,rowPtrNew,colIndNew] = rowPermuteS(value,rowPtr,colInd,1,3);
[valueNew,rowPtrNew,colIndNew] = rowPermuteS(valueNew,rowPtrNew,colIndNew,1,5);
[AfullNew] = rowPermuteF(Afull,1,3);
[AfullNew] = rowPermuteF(AfullNew,1,5);
% Operation 2: Row Scaling
[valueNew,rowPtrNew,colIndNew] = rowScaleS(valueNew,rowPtrNew,colIndNew,1,4,3);
[valueNew,rowPtrNew,colIndNew] = rowScaleS(valueNew,rowPtrNew,colIndNew,5,2,-4.4);
[AfullNew] = rowScaleF(AfullNew,1,4,3);
[AfullNew] = rowScaleF(AfullNew,5,2,-4.4);
% Operation 3: Product
[bSparse] = productAxS(valueNew,rowPtrNew,colIndNew,x);
[bFull] = productAxF(AfullNew,x);

%% Wilkinson Test
ARetrieve = retrieveElement(valueNew,rowPtrNew,colIndNew);
normTestA = norm(AfullNew - ARetrieve);
normTestb = norm(bFull - bSparse);
end