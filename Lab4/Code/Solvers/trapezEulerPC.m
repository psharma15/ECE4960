% -----------------------------------------------------------------------
% This code solves ODE by Trapezoidal Euler Method
% This is an implicit method, I'm using forward Euler to generate x at new
% time step and then using that to generate new x
% Thus, it is like a predictor corrector model
% -----------------------------------------------------------------------
% Pragya Sharma, April 2017
% ps847@cornell.edu
% -----------------------------------------------------------------------

function xt = trapezEulerPC(tmin,tmax,delT,x0,ode)
% -----------------------------------------------------------------------
% Input:
% tmin: Starting Time
% tmax: Final time
% delT: Time step (fixed for all ODEs
% x0: Vector of initial x values
% ode: Function handle to the ODE to be solved
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
    % First step by Forward Euler
    xt(:,i) = xt(:,i-1) + ode(xt(:,i-1),t(i-1))*delT; 
    xt(:,i) = xt(:,i-1) + ((ode(xt(:,i-1),t(i-1))+ode(xt(:,i),t(i)))/2)*delT;
end
end
