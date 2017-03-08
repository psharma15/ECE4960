%% ECE 4960, Lab 2
% Pragya Sharma - ps847, 04 March 2017

function [x,b,residualVector,normDelx,ConvInfNorm,iter,normRes] = solverRC()
path(path,genpath(pwd));

%% Select options here
matSize = 'small1';
maxIter = 250;
testOn = 1;
if strcmp(matSize,'large')
    bOpt = 1;
end
solverMethod = 'jacobi'; 

%% Get test data
switch matSize
    case 'small1'
        value = [10 -1 2 -1 11 -1 3 2 -1 10 -1 3 -1 8];
        rowPtr = [1 4 8 12 15];
        colInd = [1 2 3 1 2 3 4 1 2 3 4 2 3 4];
        matRank = length(rowPtr)-1;
        b = sort(1:matRank,'Descend'); 
    case 'small2'
        value = 1:12;
        rowPtr = [0,3,6,9,10,12]+1;
        colInd = [0,1,4,0,1,2,1,2,4,3,0,4]+1;
        matRank = length(rowPtr)-1;
        b = sort(1:matRank,'Descend'); 
        % Precondition with Epsilon for diagonal dominance
        epsilon = 0.001; 
    case 'large'
        % Get Memplus in row compressed format 
        [rowPtr,colInd,value,rank] = testDataMemplusRC();
        % Precondition with Epsilon for diagonal dominance
        epsilon = 10;    
        switch bOpt
            case 1
                b = [1 zeros(1,rank-1)]; 
            case 2
                b = [0 1 zeros(1,rank-2)]; 
            case 3
                b = [0 0 0 0 1 zeros(1,rank-5)]; 
            case 4
                b = [0 0 0 0 0 1 zeros(1,rank-6)]; 
            case 5
                b = [zeros(1,7) 1 zeros(1,rank-8)]; 
            case 6
                b = ones(1,rank); 
            otherwise
                warning('This b option doesn''t exist');
        end
    otherwise
        warning('Wrong test data matrix size selected');
end

%% Solve here 
switch solverMethod
    case 'jacobi'
        [x,residualVector,normDelx,ConvInfNorm,iter] = jacobiRc(value,rowPtr,colInd,b,maxIter,testOn);
    case 'gaussSeidel'
        [x,residualVector,normDelx,ConvInfNorm,iter] = gaussSeidelRc(value,rowPtr,colInd,b,maxIter,testOn);
    case 'jacobiD'
        [x,residualVector,normDelx,ConvInfNorm,iter] = jacobiRcD(value,rowPtr,colInd,b,maxIter,epsilon,testOn);
    case 'gaussSeidelD'
        [x,residualVector,normDelx,ConvInfNorm,iter] = gaussSeidelRcD(value,rowPtr,colInd,b,maxIter,epsilon,testOn);
    otherwise
        display('This solver doesn''t exist');
end
normRes = residualVector(end)/norm(b); % Normalized residual vector
if iter == maxIter
    fprintf('No Convergence. Solution not found. \n');
    return;
end
iter = iter - 1;
if testOn
    % Test 1 - Residual convergence
    for i = 1:iter-1
        temp = residualVector(i)/residualVector(i+1);
        if temp > 1.1
            fprintf('Residual Vector not converging. Test failed \n');
            return;
        end
    end
    fprintf('Residual Vector is converging. Test successful \n');
    % Test 2 - Agreement with Matlab sparse solution
    [xMat,resMat,normResMat] = solveMat(matSize,b');
    xMat = transpose(xMat);
    cmpX = norm(xMat - x); % Compare solution x
    tol = 1e-2;
    if cmpX < tol
        fprintf('The iterative solution agrees with Matlab built-in solution. \nTest Passed. \n');
    else
        fprintf('The iterative solution doesn''t agree with Matlab built-in solution. \nTest Failed. \n');
    end
    fprintf('The norm ||x-xMat||_2 is %d \n',cmpX);
    assignin('base','xMat',xMat);
    assignin('base','resMat',resMat);
    assignin('base','normResMat',normResMat);
end

end 