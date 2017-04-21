% -----------------------------------------------------------------------
% This code solves ODE by RK34 method WITHOUT time adaptivity
% It gives additional output of error Estimation and possible h at next
% time step, but does not use that
% -----------------------------------------------------------------------
% Pragya Sharma, April 2017
% ps847@cornell.edu
% -----------------------------------------------------------------------

function [xt,h,errEstim] = rk34(tmin,tmax,delT,x0,epsilonR,ode)
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
% h: Possible h at next time step
% errEstim: Error estimation performed using RK3 and RK4

%% -----------------------------------------------------------------------

lengthT = uint32((tmax-tmin)/delT + 1);
t = linspace(tmin,tmax,lengthT);
lengthX = length(x0);
xt = zeros(lengthX,lengthT);
errEstim = zeros(lengthX,lengthT);
h = zeros(lengthT,1);

% Initializing values at first time step
xt(:,1) = x0;
h(1) = delT;
k4old = ode(xt(:,1),t(1));

%% -----------------------------------------------------------------------
for i = 2:lengthT
    % ---------------- First construct a possible RK3 --------------------
    % using k1 from previous step
    k1 = k4old; 
    k2 = ode(xt(:,i-1)+(k1.*(delT/2)),t(i-1)+(delT/2));
    k3 = ode(xt(:,i-1)+(k2.*(3*delT/4)),t(i-1)+(3*delT/4));
    xt(:,i) = xt(:,i-1) + (1/9)*(2*k1+3*k2+4*k3)*delT;
    % ------- Now use this to construct RK4 with matching terms ----------
    k4 = ode(xt(:,i),t(i-1)+delT);
    xt(:,i) = xt(:,i-1) + (1/24)*(7*k1 + 6*k2 + 8*k3 + 3*k4)*delT;
    errEstim(:,i) = (1/72)*(-5*k1 + 6*k2 + 8*k3 - 9*k4)*delT;
    h(i) = delT * (epsilonR/(norm(errEstim(:,i))/norm(xt(:,i))))^(1/3);
    % updating old k4 as the new k4
    k4old = k4;
end
end