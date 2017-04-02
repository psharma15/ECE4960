% ------------------------------------------------------------------------
% Parameter Extraction from least-square fitting
% ------------------------------------------------------------------------

% ---------------------------Output generated-----------------------------
% res.iter: nIter - Total number of iterations 
% res.v: Least Square cost function vector [1 x nIter]
% res.delv: Gradient column vector [1 x nIter]
% res.aold: Parameters [nPara x (nIter + 1)]
% res.dela: Delta paramters [nPara x nIter]
% res.normDelCheck: || dela/a ||_2 [1 x nIter]
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% Pragya Sharma, March 30 2017
% ps847@cornell.edu
% ------------------------------------------------------------------------

function [res] = paraExtractLS(par)
    
    % --------------------Set up needs input from user--------------------
    % par.sMeas: Measurement data 
    % ---------- last column: Output variable
    % ---------- First to (last - 1) columns: Variables
    % par.sMod: Function handle to Model function
    % par.a0: Initial parameter input column vector
    % par.h: Numerical differentiation step size column vector
    % par.mode: 
    % ---------'QN': Quasi-Newton 
    % ---------'SC': Secant 
    % par.maxIter: Maximum number of iterations
    % par.ls: Line Search Yes or No
    % par.tol: Tolerance on converging condition
    % --------------------------------------------------------------------
    
    % --------------------------------------------------------------------
    % Setting up problem
    % ------------------
    YES = 1;
    NO = 0;
    % ------------------
    nPara = length(par.a0);
    switch(par.mode)
        case 'QN'
            delFuncv = @delFuncvQN;
            hessFuncv = @hessFuncvQN;
        case 'SC'
            delFuncv = @delFuncvSC;
            hessFuncv = @hessFuncvSC;
        otherwise
            warning('Wrong solver mode. Choose Options: QN, SC')
    end
    if (par.ls == YES)
        t = [10, 1, 0.1];
    else 
        t = 1;
    end
    flagIter = YES;
    
    res.iter = 1;
    res.v = zeros(nPara,par.maxIter);
    res.delv = zeros(nPara,par.maxIter);
    res.dela = zeros(nPara,par.maxIter);
    res.aold = zeros(nPara,par.maxIter+1);
    res.aold(:,1) = par.a0;
    res.sol = 0;
    res.normDelCheck = zeros(1,par.maxIter);
    % --------------------------------------------------------------------
    % Suppressing warning for ill-conditioned matrix
    id = 'MATLAB:illConditionedMatrix';
    warning('off',id)
    % --------------------------------------------------------------------
    while (res.iter <= par.maxIter) && (flagIter == YES)
        res.v(res.iter) = funcV(res.aold(:,res.iter),par.sMeas,par.sMod);
        res.delv(:,res.iter) = delFuncv(res.aold(:,res.iter),par.sMeas,par.sMod,par.h);
        
        % delta(a) generated 
        res.dela(:,res.iter) = -(hessFuncv(res.aold(:,res.iter),par.sMeas,par.sMod,par.h))\res.delv(:,res.iter);
        
        % Line search to choose best step
        vLineSearch = zeros(length(t),1);
        anewLineSearch = zeros(nPara,length(t));
        for iter2 = 1:length(t)
            anewLineSearch(:,iter2) = res.aold(:,res.iter) + t(iter2)*res.dela(:,res.iter);
            vLineSearch(iter2) = funcV(anewLineSearch(:,iter2),par.sMeas,par.sMod);
        end
        [~,I] = min(abs(vLineSearch));
        anew = anewLineSearch(:,I);
        res.aold(:,res.iter+1) = anew;
        res.normDelCheck(res.iter) = sum((res.dela(:,res.iter).^2)./(res.aold(:,res.iter+1).^2));
                
        % Secant method
        if strcmp(par.mode,'SC')
            par.h = res.dela(:,res.iter)./2;
        end

        % Checking convergence by NORM(relative increment vector) < tol
        if res.normDelCheck(res.iter) < par.tol
            flagIter = NO;
            res.sol = anew;
            tempSol = sprintf('%10.2e',res.sol);
            fprintf('Solution is %s\n',tempSol);            
        end
        res.iter = res.iter + 1;
        
    end
    % --------------------------------------------------------------------
    res.iter = res.iter - 1;
    res.v = res.v(1:res.iter);
    res.delv = res.delv(:,1:res.iter);
    res.dela = res.dela(:,1:res.iter);
    res.aold = res.aold(:,1:res.iter+1);
    res.normDelCheck = res.normDelCheck(:,1:res.iter);
    if res.sol == 0
        fprintf('Convergence not achieved in %d iterations. Iterations stopped.\n', par.maxIter);
    end
    
    % --------------------------------------------------------------------
    warning('on',id)
end
