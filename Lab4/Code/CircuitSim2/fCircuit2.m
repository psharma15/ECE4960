function dxdt = fCircuit2(x,t)
rg = 10e3;
rl = 10e3;
c1 = 1e-12;
c2 = 1e-12;
is = 5e-6;
kappa = 0.7;
vth = 1;
vt = 26e-3;
vdd = 5;

% -----------------------------------------------------------------------
% Current at time t
tCurrent = rem(t*1e9,20) * 1e-9; % Wrapping time to range 0-20ns
if tCurrent <= 1e-9
    i = (0.1e-3/1e-9)*tCurrent;
else if (tCurrent > 1e-9) && (tCurrent <= 10e-9)
        i = 0.1e-3;
    else if (tCurrent > 10e-9) && (tCurrent <= 11e-9)
            i = -(0.1e-3/1e-9)*tCurrent + 11*0.1*1e-3;
        else
            i = 0;
        end
    end
end

% -----------------------------------------------------------------------
dxdt = zeros(length(x),1);
dxdt(1) = (-1/(rg*c1))*x(1) + (1/(rg*c1))*i;
id = is*(log(1+exp((kappa*(x(1)-vth))/(2*vt))))^2- ... 
    is*(log(1+exp((kappa*(x(1)-vth)-x(2))/(2*vt))))^2;
dxdt(2) = (-1/c2)*id - (1/(rl*c2))*x(2) + (1/(rl*c2))*vdd;
end
