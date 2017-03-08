%% Test data: Convert to row compressed storage format 
function [rowPtr,colInd,value,rankMemplus] = testDataMemplusRC()
    memplus = importfile('memplus.mtx');
    [rowsort,I] = sort(memplus(:,1));
    colInd = memplus(I,2);
    value = memplus(I,3);
    lengthMemplus = length(memplus);
    rankMemplus = rowsort(end); %Dimension
    rowPtr = zeros(rankMemplus+1,1);
    counter = 1;
    sumInd = 0;
    for i = 1:rankMemplus
        flag = 1;
        while flag
            if rowsort(counter) == i
                sumInd = sumInd+1;
                counter = counter + 1;
                if counter == lengthMemplus+1
                    flag = 0;
                end
            else
                flag = 0;
            end
        end
        rowPtr(i+1) = sumInd;
    end
    rowPtr = rowPtr+1;
    % All need to be horizontal 
    rowPtr = rowPtr';
    colInd = colInd';
    value = value';
end