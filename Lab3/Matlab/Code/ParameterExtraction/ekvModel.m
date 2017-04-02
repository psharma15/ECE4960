% This function calculates Model Values of S, Id using ekv model
% ------------------------------------------------------------------------
% Pragya Sharma, March 30 2017
% ------------------------------------------------------------------------

function sMod = ekvModel(a,sMeas)
    VT = 26e-3;
    expr1 = (a(2)*(sMeas(:,1)-a(3)))/(2*VT);
    expr2 = ((a(2)*(sMeas(:,1)-a(3))) - sMeas(:,2))/(2*VT);
    sMod = a(1).* (log(1+exp(expr1))).^2 - a(1).* (log(1+exp(expr2))).^2; 
end
