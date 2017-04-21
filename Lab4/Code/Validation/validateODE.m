% -----------------------------------------------------------------------
% This code is for validation of ODE solvers. 
% It solves dx/dt = 4exp(0.8t) - 0.5x
% True solution is x(t) = (4/1.3)(exp(0.8t)-exp(-0.5t)) + 2exp(-0.5t)

% -----------------------------------------------------------------------
% Pragya Sharma, April 2017
% ps847@cornell.edu
% -----------------------------------------------------------------------

function [h,h2] = validateODE
% validate ODE using forward and backward Euler, Trapezoidal Euler, and
% RK34 with and without time adaptivity.

%% -----------------------------------------------------------------------
% Initial Settings
tmin = 0;
tmax = 4;
delT = 1;
lengthT = uint32((tmax-tmin)/delT + 1);
t = linspace(tmin,tmax,lengthT);
xtrue = fValidateSolution(t);
epsilonR = 1e-2;
x0 = xtrue(1);
ode = @fValidate;
% dOde = @dfValidate;
xNew = @xImplicit;

%% -----------------------------------------------------------------------
% Add Path
solverPath = fullfile(cd,'..\Solvers');
addpath(solverPath);

%% -----------------------------------------------------------------------
% Solve ODE using different solvers
xFwdEuler = forwardEuler(tmin,tmax,delT,x0,ode);
% xBwdEuler = backwardEulerNewton1(tmin,tmax,delT,x0,ode,dOde);
xBwdEuler = backwardEulerSym(tmin,tmax,delT,x0,ode,xNew);
% xTrapezEuler = trapezEulerNewton1(tmin,tmax,delT,x0,ode,dOde);
xTrapezEuler = trapezEulerSym(tmin,tmax,delT,x0,ode,xNew);
[xRK34,~,~] = rk34(tmin,tmax,delT,x0,epsilonR,ode);
[xRK34Adapt,tNew,h,~] = rk34Adaptive(tmin,tmax,delT,x0,epsilonR,ode);
[xRK34Adapt2,tNew2,h2,~] = rk34Adaptive2(tmin,tmax,delT,x0,epsilonR,ode);

%% -----------------------------------------------------------------------
% Relative normalized error
xtrueNorm = sqrt(sum(xtrue.^2,1));
normErrFwdEuler = sqrt(sum((xFwdEuler - xtrue).^2,1))./xtrueNorm;
normErrBwdEuler = sqrt(sum((xBwdEuler - xtrue).^2,1))./xtrueNorm;
normErrTrapezEuler = sqrt(sum((xTrapezEuler - xtrue).^2,1))./xtrueNorm;
normErrRK34 = sqrt(sum((xRK34 - xtrue).^2,1))./xtrueNorm;
xtrue1 = fValidateSolution(tNew);
normErrRK34Adapt = sqrt(sum((xRK34Adapt - xtrue1).^2,1))./sqrt(sum(xtrue1.^2,1));
xtrue2 = fValidateSolution(tNew2);
normErrRK34Adapt2 = sqrt(sum((xRK34Adapt2 - xtrue2).^2,1))./sqrt(sum(xtrue2.^2,1));

%% -----------------------------------------------------------------------
% Plotting Error
figure
plot(t,normErrFwdEuler,t,normErrBwdEuler,t,normErrTrapezEuler,t,normErrRK34,...
    tNew,normErrRK34Adapt,'--*',tNew2,normErrRK34Adapt2,'--o');
legend('Forward Euler','Backward Euler','Trapezoidal Euler','RK34','Adaptive RK34','Adaptive RK34 - 2','FontSize',12)
title('Normalized Error vs Time','FontSize',12)
xlabel('Time','FontSize',12)
ylabel('Normalized Error','FontSize',12)

%% -----------------------------------------------------------------------
rmpath(solverPath);

end