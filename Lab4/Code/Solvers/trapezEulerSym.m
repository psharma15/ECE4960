% -----------------------------------------------------------------------
% This code solves ODE by Trapezoidal Euler Method
% This is an implicit method, using symbolic equation to get x_i+1
% -----------------------------------------------------------------------
% Pragya Sharma, April 2017
% ps847@cornell.edu
% -----------------------------------------------------------------------

function xt = trapezEulerSym(tmin,tmax,delT,x0,ode,xNew)
% -----------------------------------------------------------------------
% Input:
% tmin: Starting Time
% tmax: Final time
% delT: Time step (fixed for all ODEs
% x0: Vector of initial x values
% ~: ode(not used): Function handle to the ODE to be solved
% xNew: Function handle to x_i+1 calculated symbolically
%
% -----------------------------------------------------------------------
% Output:
% xt: Final x, columns are result x vectors at a time instant
% -----------------------------------------------------------------------

lengthT = uint32((tmax-tmin)/delT + 1);
t = linspace(tmin,tmax,lengthT);
lengthX = length(x0);
xt = zeros(lengthX,lengthT);
xt(:,1) = x0;
for i = 2:lengthT
    % Backward Euler to get xt_i+1
    xt(:,i) = xNew(xt(:,i-1),t(i),delT);
    xt(:,i) = xt(:,i-1) + ((ode(xt(:,i-1),t(i-1))+ode(xt(:,i),t(i)))/2)*delT;
end
end
