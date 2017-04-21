% -----------------------------------------------------------------------
% This code solves ODE by Forward Euler Method
% -----------------------------------------------------------------------
% Pragya Sharma, April 2017
% ps847@cornell.edu
% -----------------------------------------------------------------------

function xt = forwardEuler(tmin,tmax,delT,x0,ode)
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
  xt(:,i) = xt(:,i-1) + ode(xt(:,i-1),t(i-1))*delT;
end
end
