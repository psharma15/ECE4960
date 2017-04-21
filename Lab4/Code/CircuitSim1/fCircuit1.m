function dxdt = fCircuit1(x,t)
r1 = 10e3;
r2 = 10e3;
r3 = 10e3;
c1 = 1e-12;
c2 = 1e-12;

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
a = -((1/(r1*c1))+(1/(r2*c2)));
b = (1/(c1*r2));
c = (1/(c2*r2));
d = -((1/(c2*r2))+(1/(c2*r3)));
e = 1/c1;
dxdt(1) = a*x(1) + b*x(2) + e*i;
dxdt(2) = c*x(1) + d*x(2);
end
