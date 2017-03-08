%% Pragya Sharma, ps847, 04th March 2017
% This function solves for x in Ax = b using Gauss Seidel iterative method with A
% in Row compressed sparse matrix format.

function [x,residualVector,normDelx,ConvInfNorm,gaussIter] = gaussSeidelRcD(value,rowPtr,colInd,b,maxIter,epsilonDD,testOn)
% A matrix is in row compressed sparse format: [value, rowPtr, colInd]
% b is same as the matrix in Ax = b 
% maxIter is maximum number of iterations
% All MUST BE row vectors.
% ------------- Introducing Diagonal Matrix Conditioning -----------------

%% Extraction of D^-1 and L+U matrices in sparse format
m = length(rowPtr);
ARank = m - 1; 
diffRow = zeros(1,m);
for iter1 = 1:m-1
    diffRow(iter1+1) = rowPtr(iter1+1)-rowPtr(iter1);    
end
nonZeroLength = length(colInd);
valueDiagInv = zeros(1,ARank);  % D inverse
nonZeroDiagonal = 0; 
indexDiagonal = zeros(1,ARank); % Indices of nonzero diagonal in colInd and value
diffRowLU = zeros(1,m);
rowPtrLU = [1 zeros(1,ARank)];
for iter2 = 1:ARank
    nonZeroLU = diffRow(iter2+1);
    for iter3 = rowPtr(iter2):rowPtr(iter2+1)-1
        if (colInd(iter3) == iter2) 
            valueDiagInv(iter2) = 1/(value(iter3) + epsilonDD);
            nonZeroDiagonal = nonZeroDiagonal + 1;      
            indexDiagonal(nonZeroDiagonal) = iter3;               
            nonZeroLU = nonZeroLU - 1;
        else
            valueDiagInv(iter2) = 1/(epsilonDD);    
        end
    end
    diffRowLU(iter2+1) = nonZeroLU;
    rowPtrLU(iter2+1) = 1 + sum(diffRowLU); 
end
indexDiagonal = indexDiagonal(1:nonZeroDiagonal);
lengthLU = nonZeroLength - nonZeroDiagonal;
valueLU = zeros(1,lengthLU); 
colIndLU = zeros(1,lengthLU);
counter1 = 1;
counter2 = 1;
iter4 = 1;
while (iter4 <= lengthLU)
    if (indexDiagonal(counter1) ~= counter2)
        valueLU(iter4) = -value(counter2);
        colIndLU(iter4) = colInd(counter2);
        iter4 = iter4 + 1;
    else
        if (indexDiagonal(counter1) == counter2)
            counter1 = counter1 + 1;
        end
    end
    counter2 = counter2 + 1;
end

%% Gauss Seidel Iteration starts here
C = valueDiagInv.*b;
x0 = C;
TrowPtr = rowPtrLU;
TcolInd = colIndLU;
Tvalue = valueLU;
rowSum = zeros(1,ARank);
for iter5 = 1:ARank                             
    ind = TrowPtr(iter5):TrowPtr(iter5+1) - 1; 
    Tvalue(ind) = valueDiagInv(iter5) .* Tvalue(ind);   % D^-1 * (L + U)
    absTvalue = abs(Tvalue(ind));
    rowSum(iter5) = sum(absTvalue,2);
end
ConvInfNorm = max(rowSum);      % Infinity norm of T matrix
tol = 1e-7; 
gaussIter = 1;
xold = x0;
flag = 1;
residualVector = zeros(1,maxIter);
normDelx = zeros(1,maxIter);
while (flag && (gaussIter < maxIter))
    xoldNew = xold;
    for iter6 = 1:ARank
        ind = TrowPtr(iter6):TrowPtr(iter6+1) - 1;
        tempVector = zeros(1,ARank);
        tempVector(TcolInd(ind)) = Tvalue(ind);
        xiNew =  tempVector*(transpose(xoldNew)) + C(iter6);
        xoldNew(iter6) = xiNew;
    end
    xnew = xoldNew;
    residualVector(gaussIter) = norm (b - productAxS(Tvalue,TrowPtr,TcolInd,xnew));
    normDelx(gaussIter) = norm(xnew - xold);
    xold = xnew;
    if normDelx(gaussIter) < tol
        flag = 0;
    end
    gaussIter = gaussIter + 1;
end
x = xold;
normDelx = normDelx(1:gaussIter-1);
residualVector = residualVector(1:gaussIter-1);
if (testOn)
    testProductAxS(Tvalue,TrowPtr,TcolInd);
end

%% CLear variables
clearvars -except x residualVector normDelx ConvInfNorm gaussIter
end