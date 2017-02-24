%% ECE 4960, HW-4
% Submitted by: Pragya Sharma, ps847

function [CheckSum] = ps847_hw4_p2()
path(path,genpath(pwd));

%% Test data: Convert to row compressed storage format 
memplus = importfile('memplus.mtx');
[rowsort,I] = sort(memplus(:,1));
colInd = memplus(I,2);
value = memplus(I,3);
L = length(memplus);
N = rowsort(end); %Dimension
rowPtr = zeros(N+1,1);
counter = 1;
sumInd = 0;
for i = 1:N
    flag = 1;
    while flag
        if rowsort(counter) == i
            sumInd = sumInd+1;
            counter = counter + 1;
            if counter == L+1
                flag = 0;
            end
        else
            flag = 0;
        end
    end
    rowPtr(i+1) = sumInd;
end
rowPtr = rowPtr+1;
x = ones(N,1);
% All need to be horizontal 
rowPtr = rowPtr';
colInd = colInd';
value = value';

%% Run operations
% Operation 1: Row Permutation
t = cputime;
userview1 = memory;
[value,rowPtr,colInd] = rowPermuteS(value,rowPtr,colInd,1,3);
[value,rowPtr,colInd] = rowPermuteS(value,rowPtr,colInd,1,5);
[value,rowPtr,colInd] = rowPermuteS(value,rowPtr,colInd,10,3000);
[value,rowPtr,colInd] = rowPermuteS(value,rowPtr,colInd,5000,10000);
[value,rowPtr,colInd] = rowPermuteS(value,rowPtr,colInd,6,15000);
% Operation 2: Row Scaling
[value,rowPtr,colInd] = rowScaleS(value,rowPtr,colInd,2,4,3);
% Operation 3: Row Permutation
[value,rowPtr,colInd] = rowPermuteS(value,rowPtr,colInd,2,5);
% Operation 4: Row Scaling
[value,rowPtr,colInd] = rowScaleS(value,rowPtr,colInd,5,4,-3);
% Operation 5: Product
[bSparse] = productAxS(value,rowPtr,colInd,x);
userview2 = memory;
e = cputime - t;
assignin('base','bSparse',bSparse);
fprintf('CPU time taken for these operations is %f sec \n',e);
% 146432 bytes is the additional required storage due to sparsity reduction and b array storage
MemUse =  (userview2.MemUsedMATLAB-userview1.MemUsedMATLAB-146432)/1e6; 
fprintf('Memory usage increase for these operations is %f Mb \n',MemUse);

%% Test
sumA = sum(value);
sumb = sum(bSparse);
CheckSum = abs(sumA - sumb);
tol = 1e-7;
if CheckSum < tol
    fprintf('Test Successful!\n');
else
    fprintf('CheckSum failed\n');
end    
end