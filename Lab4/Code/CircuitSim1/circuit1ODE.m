% -----------------------------------------------------------------------
% This code is for circuit simulation solved using ODEs

% -----------------------------------------------------------------------
% Pragya Sharma, April 2017
% ps847@cornell.edu
% -----------------------------------------------------------------------

function [h,h2,tNew2] = circuit1ODE
% Solve circuit ODE using forward, backward and Trapezoidal Euler, and
% RK34 with and without time adaptivity.

%% -----------------------------------------------------------------------
% Initial Settings - Settings 1
tmin = 0;
tmax = 100e-9;
delT = 1e-9;
t = tmin:delT:tmax;
epsilonR = 1e-2;
x0 = [0;0];
ode = @fCircuit1;
xNew = @xImplicit;
% dOde = @dfCircuit1;

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