%% ECE 4960, Lab 2
% Pragya Sharma - ps847, March 2017
% Solves Ax = b in Row Compressed format using Jacobi or Gauss Seidel
% methods.

function [x,b,residualVector,normDelx,ConvInfNorm,iter,normRes] = solverRC()
path(path,genpath(pwd));

%% Select options here
matSize = 'large';
maxIter = 20;
testOn = 1;
if strcmp(matSize,'large')
    bOpt = 3;
end
solverMethod = 'gaussSeidelD'; % Epsilon varies with this.
tolSolver = 1e-7;              % Tolerance of Solver
tolSolution = 1e-2;            % Tolerance with 'true' Solution

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
        % Here, 0.2 for gaussSeidel, 1 for Jacobi.        
        epsilon = 1;    
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
timeInitial = cputime;
memoryInitial = memory;
% Perform pivoting
[valueNew,rowPtrNew,colIndNew,bNew,chngOrder] = doPivoting(value,rowPtr,colInd,b);
% Solve using possible four methods
switch solverMethod
    case 'jacobi'
        [x,residualVector,normDelx,ConvInfNorm,iter] = jacobiRc(valueNew,rowPtrNew,colIndNew,bNew,maxIter,tolSolver,testOn);
        fprintf('Running Jacobi Solver.\n');
    case 'gaussSeidel'
        [x,residualVector,normDelx,ConvInfNorm,iter] = gaussSeidelRc(valueNew,rowPtrNew,colIndNew,bNew,maxIter,tolSolver,testOn);
        fprintf('Running Gauss Seidel.\n');
    case 'jacobiD'
        [x,residualVector,normDelx,ConvInfNorm,iter] = jacobiRcD(valueNew,rowPtrNew,colIndNew,bNew,maxIter,epsilon,tolSolver,testOn);
        fprintf('Running Jacobi with Diagonal Pre-conditioning using Epsilon = %f.\n',epsilon);
    case 'gaussSeidelD'
        [x,residualVector,normDelx,ConvInfNorm,iter] = gaussSeidelRcD(valueNew,rowPtrNew,colIndNew,bNew,maxIter,epsilon,tolSolver,testOn);
        fprintf('Running Gauss Seidel with Diagonal Pre-conditioning using Epsilon = %f.\n',epsilon);
    otherwise
        display('This solver doesn''t exist');
end

x = x(chngOrder); % Change order back to previous.

memoryFinal = memory;
timeFinal = cputime - timeInitial;
fprintf('CPU time taken for these operations is %f sec \n',timeFinal);
MemUse =  (memoryFinal.MemUsedMATLAB-memoryInitial.MemUsedMATLAB)/1e6; 
fprintf('Memory usage increase is %d Mb \n',MemUse);

normRes = residualVector(end)/norm(b); % Normalized residual vector
if iter == maxIter
    fprintf('No Convergence. Solution not found. \n');
    return;
end
iter = iter - 1;
if testOn
    % Test 2 - Residual convergence
    for i = 1:iter-1
        temp = residualVector(i)/residualVector(i+1);
        if temp > 1.01
            fprintf('Residual Vector not converging. Test failed \n');
            return;
        end
    end
    fprintf('Residual Vector is converging. Test successful \n');
    % Test 4 - Agreement with Matlab sparse solution
    [xMat,resMat,normResMat] = solveMat(matSize,b');
    xMat = transpose(xMat);
    cmpX1 = norm(xMat - x); % Compare solution x
    if cmpX1 < tolSolution
        fprintf('The iterative solution agrees with Matlab built-in sparse solution. \nTest Passed. \n');
    else
        fprintf('The iterative solution doesn''t agree with Matlab built-in sparse solution. \nTest Failed. \n');
    end
    fprintf('The norm ||x-xMat||_2 is %d \n',cmpX1);
    assignin('base','xMat',xMat);
    assignin('base','resMat',resMat);
    assignin('base','normResMat',normResMat);
    if ~(strcmp(matSize,'large'))
        % Test 5 - For Small Matrices, comparing with Full matrix solution
        xFullMat = solveFull(value,rowPtr,colInd,b);
        cmpX2 = norm(xFullMat' - x);
        if cmpX2 < tolSolution
            fprintf('The iterative solution agrees with full matrix solution. \nTest Passed. \n');
        else
            fprintf('The iterative solution doesn''t agree with full matrix solution. \nTest Failed. \n');
        end
        % Test 6 - For Small Matrices, comparing pivoting result with full
        % matrix solution
        if strcmp(matSize,'small2')
            Apivot = [12 0 0 0 11; 0 10 0 0 0;9 0 8 7 0;0 0 6 5 4; 3 0 0 2 1];
            ApivotSparse = retrieveElement(value,rowPtr,colInd);
            if norm(Apivot - ApivotSparse) < tolSolution
                fprintf('Sparse pivoting solution matches the true solution\n');
            end
        end
                  
    end
   
    
end

end 