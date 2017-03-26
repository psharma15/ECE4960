%% Pragya Sharma, ps847, March 09 2017
% This code interchanges columns i and j, where i<j for Sparse Row-Compressed
% Matrix
function [valueNew,rowPtrNew,colIndNew] = colPermuteS(value,rowPtr,colInd,i,j)
valueNew = value;
colIndNew = colInd;
aRank = length(rowPtr) - 1;
rowPtrNew = rowPtr;
for iter0 = 1:aRank
    temp1 = 0;
    temp2 = 0;
    startPtr = rowPtr(iter0);
    endPtr = rowPtr(iter0+1)-1;
    for iter1 = startPtr:endPtr
        if colInd(iter1) == i
            colIndNew(iter1) = j;
            temp1 = (iter1);
        else if colInd(iter1) == j
                colIndNew(iter1) = i;
                temp2 = (iter1);
            end
        end
    end
    colIndNew(startPtr:endPtr) = sort(colIndNew(startPtr:endPtr),'ascend');
    if (temp1 == 0) && (temp2 == 0)
        continue;
    else if (temp1 == 0) && temp2
            flag2 = 1;
            for iter2 = startPtr:temp2
                if (colInd(iter2) > i) && (flag2 == 1)
                    valueNew(iter2) = value(temp2);
                    flag2 = 0;
                    continue; % Should skip that loop
                end
                if (colInd(iter2) > i) && (flag2 == 0)
                    valueNew(iter2) = value(iter2-1);
                end
            end
        else if temp1 && (temp2 == 0)
            for iter2 = temp1:endPtr
                if iter2 ~= endPtr
                    if colInd(iter2+1)<j
                        valueNew(iter2) = value(iter2+1);
                    else
                        valueNew(iter2) = value(temp1);
                        break; % Should exit the loop
                    end
                else
                    valueNew(iter2) = value(temp1);
                end
                
            end
            else
                valueNew(temp1) = value(temp2);
                valueNew(temp2) = value(temp1);
            end
        end
    end
end
end
        
