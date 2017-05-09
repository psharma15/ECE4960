%% This code takes sparse matrix in Row Compressed form, do pivoting.
% Pragya Sharma, March 09, 2017
% Bug detected: May 09 2017 - If first entry (row 1, column 1) is largest,
% it is stuck in loop without any changes.
% Bug corrected: May 09 2017 - Pragya Sharma

function [valueNew,rowPtrNew,colIndNew,bNew,chngOrder] = doPivoting(value,rowPtr,colInd,b)
% With row, change b
% With column, record indices to permute x vector values
% path(path,genpath(pwd));
valueTemp = value;
rowPtrTemp = rowPtr;
colIndTemp = colInd;
aRank = length(rowPtr) - 1;
chngOrder = 1:aRank;        % Order of x vector
aRankTemp = aRank;
valueNew = value;
rowPtrNew = rowPtr;
colIndNew = colInd;
bNew = b;
for iter2 = 1:aRank-1
    % Finding maximum value that should be a pivot
    [~,ind] = max(valueTemp);
    colNum = colIndTemp(ind);
    for iter3 = 1:aRankTemp
        if (rowPtrTemp(iter3) <= ind) && (ind < rowPtrTemp(iter3 + 1))
            rowNum = iter3;
        end
    end
    if rowNum ~= 1        
        [valueTemp,rowPtrTemp,colIndTemp] = rowPermuteS(valueTemp,rowPtrTemp,colIndTemp,1,rowNum);
        [valueNew,rowPtrNew,colIndNew] = rowPermuteS(valueNew,rowPtrNew,colIndNew,iter2,rowNum+iter2-1);
        bNew(iter2) = b(rowNum + iter2 -1);
        bNew(rowNum+ iter2 -1) = b(iter2);
    end
    if colNum ~= 1
        [valueTemp,rowPtrTemp,colIndTemp] = colPermuteS(valueTemp,rowPtrTemp,colIndTemp,1,colNum);
        [valueNew,rowPtrNew,colIndNew] = colPermuteS(valueNew,rowPtrNew,colIndNew,iter2,colNum+iter2-1);
        chngOrder(iter2) = colNum+iter2-1;
        chngOrder(colNum+iter2-1) = iter2;
    end
    % --------------------------------BUG----------------------------------
%     if (rowNum==1) && (colNum==1)
%         % This is the case when element in the first row and column is the
%         % pivot element. Do Nothing.
%         break;
%     end
    % ---------------------------------------------------------------------
    aRankTemp = aRankTemp - 1;
    diffRow = zeros(1,aRankTemp+1);
    for iter4 = 2:aRankTemp+1
        diffRow(iter4) = rowPtrTemp(iter4+1)-rowPtrTemp(iter4); 
        subCol = 0;
        for iter5 = rowPtrTemp(iter4):rowPtrTemp(iter4+1)-1
            if colIndTemp(iter5) == 1
                subCol = 1;
                continue;
            end
        end
        diffRow(iter4) = diffRow(iter4) - subCol;       
    end
    valueTemp = valueTemp(rowPtrTemp(2):end);
    colIndTemp = colIndTemp(rowPtrTemp(2):end); % This will just take out row1
    [~,indColOne] = find(colIndTemp == 1);

    rowPtrTemp = [1 zeros(1,aRankTemp)];
    for iter6 = 1:aRankTemp
        rowPtrTemp(iter6+1) = rowPtrTemp(iter6) + diffRow(iter6+1);    
    end
    

    colIndTemp(indColOne)=[];
    colIndTemp = colIndTemp - 1;

    valueTemp(indColOne) = [];
%     retrieveElement(valueNew,rowPtrNew,colIndNew)
end
end