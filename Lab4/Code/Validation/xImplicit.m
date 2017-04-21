function xNew = xImplicit(x,t,delT)
xNew = (x + 4*delT*exp(0.8*t))/(1+0.5*delT);
end