% -----------------------------------------------------------------------
% This code solves ODE by RK34 method WITH time adaptivity. 
% This adapts BOTH current and forward time step.
% -----------------------------------------------------------------------
% Pragya Sharma, April 2017
% ps847@cornell.edu
% -----------------------------------------------------------------------

function [xt,t,h,errEstim] = rk34Adaptive2(tmin,tmax,delT,x0,epsilonR,ode)
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
maxIterTime = lengthT; % Approximating number of time points 
maxIter = 100; % Iterations for correcting error at one time step

%% -----------------------------------------------------------------------
iterTime = 2;
while iterTime <= maxIterTime
    iter = 0;
    while (iter < maxIter)
        % ---------------- First construct a possible RK3 --------------------
        % using k1 from previous step
        k1 = k4old; 
        k2 = ode(xt(:,iterTime-1)+(k1*h(iterTime-1)/2),t(iterTime-1)+(h(iterTime-1)/2));
        k3 = ode(xt(:,iterTime-1)+(3*k2*h(iterTime-1)/4),t(iterTime-1)+(3*h(iterTime-1)/4));
        xt(:,iterTime) = xt(:,iterTime-1) + (1/9)*(2*k1+3*k2+4*k3)*h(iterTime-1);
        % ------- Now use this to construct RK4 with matching terms ----------
        k4 = ode(xt(:,iterTime),t(iterTime-1)+h(iterTime-1));
        errEstim(:,iterTime) = (1/72)*(-5*k1 + 6*k2 + 8*k3 - 9*k4)*h(iterTime-1);
        xt(:,iterTime) = xt(:,iterTime-1) + (1/24)*(7*k1 + 6*k2 + 8*k3 + 3*k4)*h(iterTime-1);
        normErr = norm(errEstim(:,iterTime))/norm(xt(:,iterTime));
        % In some cases decreasing h is not decreasing error. For that,
        % compare normErrNew to normErrPrev. If it is decreasing, continue,
        % else take the solution with delT fixed, and exit. 
        if (normErr > epsilonR) 
            hTemp = h(iterTime-1);
            h(iterTime-1) = h(iterTime-1) * (epsilonR/normErr)^(1/3);
            iter = iter + 1;                
        else
            break;
        end
    end
    h(iterTime) = h(iterTime-1) * (epsilonR/normErr)^(1/3);
    t(iterTime) = t(iterTime-1) + h(iterTime);

    % updating old k4 as the new k4
    k4old = k4;
    if (t(iterTime) >= tmax) && (t(iterTime-1) < tmax)
        maxIterTime = iterTime + 1 ;
    end
    iterTime = iterTime + 1;
end
t = t(1:iterTime-1);
xt = xt(:,1:iterTime-1);
h = h(1:iterTime-1);
errEstim(:,1:iterTime-1);
end