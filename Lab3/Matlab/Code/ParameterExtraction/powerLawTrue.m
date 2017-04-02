% This function calculates Measurement Values of S, logy using power law
% y = c0*x^m is converted to log(y) = log(c0) + m*log(x)
% ------------------------------------------------------------------------
% Pragya Sharma, March 30 2017
% ------------------------------------------------------------------------

function sMeas = powerLawTrue(x,c0,m)
    c = log10(c0);
    nSample = length(x);
    sMeas = c.* ones(nSample,1) + m.*x; 
end
