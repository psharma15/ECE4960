% -----------------------------------------------------------------------
% This code is for circuit simulation. It takes care of ODE solution part.
% The circuit is Trasnsitor CS amplifier using EKV model for current.
%
% -----------------------------------------------------------------------
% Pragya Sharma, April 2017
% ps847@cornell.edu
% -----------------------------------------------------------------------

function [h,h2] = circuit2ODE
% CIRCUIT2ODE solves ODE using forward, backward and Trapezoidal Euler, and
% RK34 with and without time adaptivity.

%% -----------------------------------------------------------------------
% Initial Settings - Settings 1
tmin = 0;
tmax = 100e-9;
delT = 0.2e-9;
t = tmin:delT:tmax;
epsilonR = 1e-3;
x0 = [2.5;2.5];
ode = @fCircuit2;
dOde = @dfCircuit2;

%% -----------------------------------------------------------------------
% Add Path
solverPath = fullfile(cd,'..\Solvers');
addpath(solverPath);

%% -----------------------------------------------------------------------
% Solve ODE using different solvers
xFwdEuler = forwardEuler(tmin,tmax,delT,x0,ode);
xBwdEuler = backwardEulerNewton1(tmin,tmax,delT,x0,ode,dOde);
xTrapezEuler = trapezEulerNewton1(tmin,tmax,delT,x0,ode,dOde);
[xRK34,~,~] = rk34(tmin,tmax,delT,x0,epsilonR,ode);
[xRK34Adapt,tNew,h,~] = rk34Adaptive(tmin,tmax,delT,x0,epsilonR,ode);
[xRK34Adapt2,tNew2,h2,~] = rk34Adaptive2(tmin,tmax,delT,x0,epsilonR,ode);

%% -----------------------------------------------------------------------
% Plotting Error
t = t*1e9;

figure
plot(t,xFwdEuler(1,:),t,xFwdEuler(2,:))
legend('V1','V2')
xlabel('Time (ns)')
ylabel('Voltage')
title('Forward Euler')

figure
plot(t,xBwdEuler(1,:),t,xBwdEuler(2,:))
legend('V1','V2')
xlabel('Time (ns)')
ylabel('Voltage')
title('Backward Euler')

figure
plot(t,xTrapezEuler(1,:),t,xTrapezEuler(2,:))
legend('V1','V2')
xlabel('Time (ns)')
ylabel('Voltage')
title('Trapezoidal Euler')

figure
plot(t,xRK34(1,:),t,xRK34(2,:))
legend('V1','V2')
xlabel('Time')
ylabel('Voltage')
title('RK34')

tNew = tNew * 1e9;
figure
plot(tNew,xRK34Adapt(1,:),tNew,xRK34Adapt(2,:))
legend('V1','V2')
xlabel('Time (ns)')
ylabel('Voltage')
title('RK34 Adaptive')

tNew2 = tNew2 * 1e9;
figure
plot(tNew2,xRK34Adapt2(1,:),tNew2,xRK34Adapt2(2,:))
legend('V1','V2')
xlabel('Time (ns)')
ylabel('Voltage')
title('RK34 Adaptive - 2')

%% -----------------------------------------------------------------------
rmpath(solverPath);

end