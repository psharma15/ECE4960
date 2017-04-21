% ------------------------------------------------------------------------
% Validation is performed for y = 10*x^-0.5 (sMeas = log(y)), power law.
% ------------------------------------------------------------------------
% ECE 4960, Lab 3, Task 2 - Validation.
% Pragya Sharma, ps847@cornell.edu, March 2017
% ------------------------------------------------------------------------

function [res,par] = paraExtractValidation()
    % Set up input
    % --------------------------
    YES = 1;
    NO = 0;
    currentPath = cd;
    addpath([currentPath,'/ParameterExtraction']);
    % --------------------------
    nSample = 10;
    x = logspace(0,3,nSample)';
    c0 = 10;
    m = -0.5;
    truePara = [c0; m];
    sMeas = [x,powerLawTrue(x,c0,m)];
    par.sMeas = sMeas;
    par.sMod = @powerLawModel;
    par.a0 = [0; 0];
    par.mode = 'QN';
    par.ls = YES;
    par.maxIter = 100;
    par.tol = 1e-9;
    switch (par.mode)
        case 'QN'
            par.h = [1e-6; 1e-6];
        case 'SC'
            par.h = [1e-6; 1e-6];
        otherwise
            warning('Choose correct par.mode');
    end

    % Calling function for parameter extraction 
    [res] = paraExtractLS(par);

    % Taking inverse log to revert back to exponential power law parameters
    if res.sol ~= 0
        res.sol(1) = 10^res.sol(1);
    end

    % Generating validation result
    if norm(res.sol - truePara) < par.tol
        fprintf('Parameter Extraction using Least Squares fitting');
        fprintf(' is working correctly.\n');
    end
end