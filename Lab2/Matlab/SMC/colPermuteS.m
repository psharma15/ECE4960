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
        % The algoritm is based on the idea that if an element exists at
        % that column index, change it's column index to the column it is
        % to be swapped with. If there's no element at one location, then
        % do nothing.
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
        % This is the case when nothing is changed because no element
        % exists in either column, do nothing.
        continue;
    else if (temp1 == 0) && temp2
            % This is the case when no element in i-th column
            flag2 = 1;
            for iter2 = startPtr:temp2
                % Do not change 'value' when column idx < i. As column idx
                % = i doesn't exist, other case is column idx > i. When
                % that happens for the first time, 'value' at column idx j
                % is shifted to current iteration, and remaining values are
                % all shifted forward by one index, until j-th column index
                % is reached, that's why iteration ends at 'temp2'
                if (colInd(iter2) > i) && (flag2 == 1)
                    valueNew(iter2) = value(temp2);
                    flag2 = 0;
                    continue; % Should skip that loop
                end
                if (colInd(iter2) > i) && (flag2 == 0)
                    % The first condition check should not be needed.
                    valueNew(iter2) = value(iter2-1);
                end
            end
        else if temp1 && (temp2 == 0)
                % This is the case when no element in j-th column
            for iter2 = temp1:endPtr
                % Shift 'value' backward by 1 after i-th column idx (temp1)
                % until j-th column idx is reached, as we are removing i-th
                % column entry. When next column idx is more than j-th idx,
                % put the 'value' from i-th column idx. This doesn't
                % include the case when no column entries after j-th
                % column. In that case, put 'value' from i-th column idx at
                % the last entry.
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
                % This is the simple case when both i-th and j-th column
                % indices are filled. Just exchange. 
                valueNew(temp1) = value(temp2);
                valueNew(temp2) = value(temp1);
            end
        end
    end
end
end
        
