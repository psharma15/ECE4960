% ------------------------------------------------------------------------
% Performing parameter extraction for EKV model of MOSFET drain current
% ------------------------------------------------------------------------
% ECE 4960, Lab 3, Task 4
% Pragya Sharma, March 2017
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% EKV model parameter extraction
% ------------------------------------------------------------------------

function [res,par] = paraExtractEKV()

    % Solution part 
    %% With ekv model (Task 4)
    % --------------------------------------------------------------------
    YES = 1;
    NO = 0;
    currentPath = cd;
    addpath([currentPath,'/ParameterExtraction']);
    
    % --------------------------------------------------------------------
    load('nmosData.mat')
    par.sMeas = sMeas;
    par.sMod = @ekvModel;
    par.a0 = [8.9e-7; 0.89; 0.75]; % Converging for: [8.9e-7; 0.89; 0.75]
    par.h = [1e-11; 1e-05; 1e-05];
    par.mode = 'QN';
    par.maxIter = 200;
    par.ls = YES;
    par.tol = 1e-9;
    
    tempa0 = sprintf('%10.2e',(par.a0)');
    fprintf('Running Parameter Extraction for EKV Model with initial point: %s\n',tempa0);
    [res] = paraExtractLS(par);
        
    % --------------------------------------------------------------------
%     %% With ekv model Normalized wrt true current in sMeas (Task 5)
%     parNorm.sMeas = sMeas;
%     parNorm.sMeas(:,end) = 1;
%     parNorm.sMod = @ekvModel;
%     parNorm.a0 = [8.9e-7; 0.89; 0.75]; % Converging for: [8.9e-7; 0.89; 0.75]
%     parNorm.h = [1e-11; 1e-05; 1e-05];
%     parNorm.mode = 'QN';
%     parNorm.maxIter = 200;
%     parNorm.ls = YES;
%     parNorm.tol = 1e-9;
% 
%     tempa0 = sprintf('%10f',(parNorm.a0)');
%     fprintf('Running Parameter Extraction for Normalized EKV Model with initial point: %f\n',tempa0);
%     [resNorm] = paraExtractLS(parNorm);
%     
%     % --------------------------------------------------------------------
%     %% Searching through initial points (Task 6)
%     Is0 = [1e-8,3e-8,1e-7,3e-7,1e-6,3e-6,1e-5,3e-5];
%     K0 = 0.2:0.1:0.8;
%     Vth0 = 1.1:0.1:2.0;
%     vMin = zeros(length(Is0)*length(K0)*length(Vth0),4);
%     count = 1;
%     parSearch = par;
%     for i = 1:length(Is0)
%         for j = 1:length(K0)
%             for k = 1:length(Vth0)
%                 parSearch.a0 = [Is0(i);K0(j);Vth0(k)];
%                 [resSearch] = paraExtractLS(parSearch);
%                 if resSearch.sol~= 0
%                     vMin(count) = [resSearch.v(end),(parSearch.a0)'];
%                 end
%                 count = count + 1;
%             end
%         end
%     end

    %% Post solution checks (for NOT normalized) (Task 4, 7)
    if res.sol ~= 0
        % --------------------------------------------------------------------
        % Observing if quadratic convergence in v and ||del||_2
        fprintf('Observing if quadratic convergence of v and ||del||_2.\n');
        fprintf('Ratio of v and ||del||_2 in consecutive iterations: \n');

        c1 = zeros(res.iter-1,1);
        c2 = zeros(res.iter-1,1);
        for i = 1:res.iter-1
            c1(i) = norm(res.v(i+1))/norm(res.v(i));
            c2(i) = res.normDelCheck(i+1)/res.normDelCheck(i);
            fprintf('norm(v(%d)/v(%d)) = %f \n',i+1,i,c1(i));
            fprintf('||del||_2(%d)/||del||_2(%d) = %f \n',i+1,i,c2(i));
        end

        % --------------------------------------------------------------------
        % Parameter Sensitivity for Vgs = 1, Vds = 3
        Vgs = 1;
        Vds = 3;
        Smeas = [Vgs,Vds];
        delIs = par.h(1);
        delK = par.h(2);
        delVth = par.h(3);
        a = res.sol + [delIs;0;0];
        delSIs = (ekvModel(a,Smeas)/ekvModel(res.sol,Smeas))/((res.sol(1)+delIs)/res.sol(1));
        a = res.sol + [0;delK;0];
        delSK = (ekvModel(a,Smeas)/ekvModel(res.sol,Smeas))/((res.sol(2)+delK)/res.sol(2));
        a = res.sol + [0;0;delVth];
        delSVth = (ekvModel(a,Smeas)/ekvModel(res.sol,Smeas))/((res.sol(3)+delVth)/res.sol(3));
        res.delS = [delSIs;delSK;delSVth];
        fprintf('Sensitivity for I_s: %f\n',delSIs);
        fprintf('Sensitivity for K: %f\n',delSK);
        fprintf('Sensitivity for Vth: %f\n',delSVth);
        
        % --------------------------------------------------------------------
        % Visualization and tests(Task 7)
        plotIVcheck(res.sol,par.sMeas)
         
        
    end  
end

