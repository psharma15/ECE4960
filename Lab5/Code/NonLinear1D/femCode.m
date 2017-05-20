function femCode()
%FEMCODE is Finite Element Code to solve Boundary Value Problem 
% ------------------------------------------------------------------------
% -----                 Equation: u'' + uu' = exp(2x)                -----
% -----                      u(0) = 1, u(1) = 0                      -----
% ------------------------------------------------------------------------
% Pragya Sharma, 05-17-2017
% Project 5, ECE 4960
% ------------------------------------------------------------------------

%% Setting up problem
n = 5;                  % Number of elements
nn = n + 1;             % Number of nodes
xmax = 1;               % Maximum limit of x
xmin = 0;               % Minimum limit of x
h = (xmax - xmin)/n;    % Step size (Length of each element)
x = xmin:h:xmax;        % Data points (independent variable)
tol = 5e-5;             % Tolerance
f = zeros(nn,1);        % Initializing FEM result vector
f(1) = exp(0);          % Implementing Boundary Conditions
f(nn) = exp(1);

%% Starting iterative method for non-linear problem
% ------------------------------------------------------------------------
iterFlag = 1;           % Flag to stop/ continue iterations
count = 0;              % Number of iterations
tic                     % Starting time measurement
testOn = 1;             % Testing individual modules in ASSEMBLY
memoryInitial = memory; % Memory usage till this point in code

while(iterFlag)
    f1 = assembly(f,n,h,testOn,tol);
    testOn = 0;
    iterFlag = 0;
    for i = 1:nn
        if (abs(f(i) - f1(i)) > tol)
            iterFlag = iterFlag + 1;
            break;
        end
    end
    f = f1;
    count = count + 1;
end

memoryFinal = memory;                 % Memory usage after the computation
err = abs(f - exp(x)');               % Calculating error
timeCalc = toc;                       % Time taken for entire operation
MemUse =  (memoryFinal.MemUsedMATLAB - ...
    memoryInitial.MemUsedMATLAB)/1e6; % Total memory usage increase

%% Verification: FEM with exact solution
% ------------------------------------------------------------------------
fprintf('Number of elements = %d\n',n);
fprintf('%12s %12s %12s %12s\n','x','FEM','Exact','Error');
fprintf('%12f %12f %12f %12f\n',x',f,exp(x)',err)
fprintf('Number of iterations = %d\n',count);

plot(x,f,'--rs','LineWidth',2)
xlabel('x')
ylabel('u(x)')
title('Solution plot using FEM');

fprintf('Time taken is = %f sec\n',timeCalc);
fprintf('Memory usage increase is %d Mb \n',MemUse);
end

