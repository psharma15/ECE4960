% -----------------------------------------------------------------------
% This code solves ODE by Trapezoidal Euler Method
% This is an implicit method, using Newton's method to solve symbolic
% equation numerically
% -----------------------------------------------------------------------
% Pragya Sharma, April 2017
% ps847@cornell.edu
% -----------------------------------------------------------------------

function xt = trapezEulerNewton1(tmin,tmax,delT,xInit,ode,dOde)
% -----------------------------------------------------------------------
% Input:
% tmin: Starting Time
% tmax: Final time
% delT: Time step (fixed for all ODEs
% x0: Vector of initial x values
% ode: Function handle to the ODE to be solved
% dOde: Function handle to Jacobian needed for Newton method
%
% -----------------------------------------------------------------------
% Output:
% xt: Final x, columns are result x vectors at a time instant
% -----------------------------------------------------------------------

lengthT = uint32((tmax-tmin)/delT + 1);
t = linspace(tmin,tmax,lengthT);
lengthX = length(xInit);
xt = zeros(lengthX,lengthT);
xt(:,1) = xInit;
for i = 2:lengthT
    % Backward Euler to get xt_i+1
    x0 = xt(:,i-1);
    x1 = x0 - inv(eye(lengthX)-delT*dOde(x0,t(i-1),ode))*(x0-delT*ode(x0,t(i-1))-xt(:,i-1));
    while (norm(x1-x0)>1e-5)
        x0 = x1;
        x1 = x0 - inv(eye(lengthX)-delT*dOde(x0,t(i-1),ode))*(x0-delT*ode(x0,t(i-1))-xt(:,i-1));
    end
    xt(:,i) = x1;
    xt(:,i) = xt(:,i-1) + ((ode(xt(:,i-1),t(i-1))+ode(xt(:,i),t(i)))/2)*delT;
end
end
