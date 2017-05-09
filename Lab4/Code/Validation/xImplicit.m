function xNew = xImplicit(x,tNew,delT)
xNew = (x + 4*delT*exp(0.8*tNew))/(1+0.5*delT);
end