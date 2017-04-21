% -----------------------------------------------------------------------
% This code solves ODE by RK34 method WITH time adaptivity at next step.
% It gives additional output of error Estimation and h for next time step.
% -----------------------------------------------------------------------
% Pragya Sharma, April 2017
% ps847@cornell.edu
% -----------------------------------------------------------------------

function [xt,t,h,errEstim] = rk34Adaptive(tmin,tmax,delT,x0,epsilonR,ode)
% -----------------------------------------------------------------------
% Input:
% tmin: Starting Time
% tmax: Final time
% delT: Time step (fixed for all ODEs
% x0: Vector of initial x values
% ode: Function handle to the ODE to be solved
% epsilonR: Normalized relative accuracy ~ 1e-2
% -----------------------------------------------------------------------
% Output:
% xt: Final x, columns are result x vectors at a time instant
% t: Time instants at which claculated
% h: Possible h at next time step
% errEstim: Error estimation performed using RK3 and RK4

%% -----------------------------------------------------------------------
% Approximating length of time vector and time at each step
lengthT = uint32((tmax-tmin)/delT + 1);
t = linspace(tmin,tmax,lengthT);
lengthX = length(x0);
xt = zeros(lengthX,lengthT);
errEstim = zeros(lengthX,lengthT);
h = zeros(lengthT,1);

% Initializing values 
xt(:,1) = x0;
h(1) = delT;
k4old = ode(xt(:,1),t(1));
maxIter = lengthT; % Approximate number of iterations

%% -----------------------------------------------------------------------
iter = 2;
while iter <= maxIter
    % ---------------- First construct a possible RK3 --------------------
    % using k1 from previous step
    k1 = k4old; 
    k2 = ode(xt(:,iter-1)+(k1*h(iter-1)/2),t(iter-1)+(h(iter-1)/2));
    k3 = ode(xt(:,iter-1)+(3*k2*h(iter-1)/4),t(iter-1)+(3*h(iter-1)/4));
    xt(:,iter) = xt(:,iter-1) + (1/9)*(2*k1+3*k2+4*k3)*h(iter-1);
    % ------- Now use this to construct RK4 with matching terms ----------
    k4 = ode(xt(:,iter),t(iter-1)+h(iter-1));
    errEstim(:,iter) = (1/72)*(-5*k1 + 6*k2 + 8*k3 - 9*k4)*h(iter-1);
    xt(:,iter) = xt(:,iter-1) + (1/24)*(7*k1 + 6*k2 + 8*k3 + 3*k4)*h(iter-1);
    normErr = norm(errEstim(:,iter))/norm(xt(:,iter));
    h(iter) = h(iter-1) * (epsilonR/normErr)^(1/3);
    t(iter) = t(iter-1) + h(iter);
    % updating old k4 as the new k4
    k4old = k4;
    if (t(iter) >= tmax) && (t(iter-1) < tmax)
        maxIter = iter + 1 ;
    end
    iter = iter + 1;
end
t = t(1:iter-1);
xt = xt(:,1:iter-1);
h = h(1:iter-1);
errEstim(:,1:iter-1);
end