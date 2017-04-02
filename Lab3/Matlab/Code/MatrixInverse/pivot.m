% ------------------------------------------------------------------------
% Pragya Sharma, March 30, 2017
% ps847@cornell.edu

% ------------------------------------------------------------------------
% This function takes augmented matrix as input and does partial pivoting
% by doing only row permutation of rowNum

% ------------------------------------------------------------------------
function Aaug = pivot(AaugTemp,rowNum)
    [~,F] = max( abs( AaugTemp(rowNum:end,rowNum) ) );
    if F > 1
        AaugTemp([F + rowNum - 1, rowNum],:) = AaugTemp([rowNum,F + rowNum - 1],:);
        Aaug = AaugTemp;
    else
        Aaug = AaugTemp;
    end
end