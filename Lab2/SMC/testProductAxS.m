% Pragya Sharma, March 2017
function testProductAxS(value,rowPtr,colInd)
%% This function tests if productAxS.m is working correctly
%  The function finds b = A*x, where x_i = 1
%  Testing is done by finding || sum(a_ij) - b ||2 < tolerance

%% Testing here
x = ones(1,length(rowPtr)-1);
b = productAxS(value,rowPtr,colInd,x);
sumA = sum(value);
sumb = sum(b);
CheckSum = abs(sumA - sumb);
tol = 1e-7;
if CheckSum < tol
    fprintf('Sparse Ax product Test Successful!\n');
else
    fprintf('Sparse Ax product Test  failed\n');
end   
end