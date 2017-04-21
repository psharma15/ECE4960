% Plot current
tSimulate = (0:0.2:100)*1e-9;
lengthT = length(tSimulate);
i = zeros(1,lengthT);

for iter = 1:lengthT
    t = tSimulate(iter);
    tCurrent = rem(t*1e9,20) * 1e-9; % Wrapping time to range 0-20ns
    if tCurrent <= 1e-9
        i(iter) = (0.1e-3/1e-9)*tCurrent;
    else if (tCurrent > 1e-9) && (tCurrent <= 10e-9)
            i(iter) = 0.1e-3;
        else if (tCurrent > 10e-9) && (tCurrent <= 11e-9)
                i(iter) = -(0.1e-3/1e-9)*tCurrent + 11*0.1*1e-3;
            else
                i(iter) = 0;
            end
        end
    end
end

figure
plot(tSimulate,i);
