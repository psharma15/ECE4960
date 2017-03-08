%% Pragya Sharma, ps847, 06th March 2017
% This function solves for x in Ax = b of A, stored as sparse in matlab

function [x,res,normRes] = solveMat(matSize,b)
switch matSize
    case 'small1'
        A = [10,-1,2,0;-1,11,-1,3;2,-1,10,-1;0,3,-1,8];
        A = sparse(A);
    case 'small2'
        A = [1,2,0,0,3;4,5,6,0,0;0,7,8,0,9;0,0,0,10,0;11,0,0,0,12];
        A = sparse(A);
    case 'large'
        % Get Memplus in row compressed format 
        [A,~] = testDataMemplusMat(); 
    otherwise
        warning('Wrong test data matrix size selected');
end
x = A\b;
res = norm(b-A*x);
normRes = res/norm(b);
end
