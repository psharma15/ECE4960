% ------------------------------------------------------------------------
% Inverse of a square non-singular matrix using Gauss Elimination
% With partial pivoting using pivot function

% ------------------------------------------------------------------------
% Pragya Sharma, March 30, 2017
% ps847@cornell.edu
% ------------------------------------------------------------------------

function invA = matInverse(A)
    % --------------------------------------------------------------------
    % Function call: invA = matInverse(A);
    % Example:
    %
    % A = [1,2,0,0,3;4,5,6,0,0;0,7,8,0,9;0,0,0,10,0;11,0,0,0,12]
    % invA = matInverse(A)
    % 
    % invA = 
    %        0.0125     0.1003    -0.0752         0     0.0533
    %        0.5110     0.0878    -0.0658         0    -0.0784
    %       -0.4342     0.0266     0.1050         0     0.0298
    %             0          0          0    0.1000          0
    %       -0.0115    -0.0920     0.0690         0     0.0345
    % --------------------------------------------------------------------
    
    %% Creating Augmented matrix
    [dim,~] = size(A);
    I = eye(dim);
    Aaug = [A,I];
    
    %% Gauss Elimination
    for i = 1:dim-1
        Aaug = pivot(Aaug,i);
        Aaug(i+1:dim,:) = Aaug(i+1:dim,:) - ( Aaug(i+1:dim,i)./Aaug(i,i) ) * Aaug(i,:);
    end
    
    %% Inverse matrix
    invA = zeros(dim);
    for i = dim : -1 : 1
        invA(i,:) = ( Aaug(i,dim+1:end) - Aaug(i,1:dim)*invA ) / Aaug(i,i);
    end
    
end


