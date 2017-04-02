% This function calculates Model Values of S, logy using power law
% y = c0*x^m is converted to log(y) = log(c0) + m*log(x)
% ------------------------------------------------------------------------
% Pragya Sharma, March 30 2017
% ------------------------------------------------------------------------

function sMod = powerLawModel(a,sMeas)
    sMod = a(1) + a(2).*sMeas(:,1);
end
