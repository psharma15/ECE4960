%% Test data: Convert to Matlab sparse storage format
function [A,rankMemplus] = testDataMemplusMat()
    memplus = importfile('memplus.mtx');
    [colInd,I] = sort(memplus(:,2));
    rowPtr = memplus(I,1);
    value = memplus(I,3);
    rankMemplus = 17758; % Dimension
    A = sparse(rowPtr,colInd,value,rankMemplus,rankMemplus);
end