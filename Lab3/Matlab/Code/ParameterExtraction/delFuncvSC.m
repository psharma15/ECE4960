% This function calculates gradient of V using Secant method
% ------------------------------------------------------------------------
% Pragya Sharma, March 30 2017
% ------------------------------------------------------------------------

function delv = delFuncvSC(a,sMeas,sModel,h)
    delv = (funcV(a+h,sMeas,sModel) - funcV(a-h,sMeas,sModel))./(2.*h)
end
