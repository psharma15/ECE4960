% This function calculates gradient of V using Quasi-Newton method
% ------------------------------------------------------------------------
% Pragya Sharma, March 30 2017
% ------------------------------------------------------------------------

function delv = delFuncvQN(a,sMeas,sModel,h)
    nPara = length(a);
    delv = zeros(nPara,1);
    
    for i = 1:nPara
        hx = zeros(nPara,1);
        hx(i) = h(i);
        delv(i) = (funcV(a+hx,sMeas,sModel) - funcV(a-hx,sMeas,sModel))/(2*hx(i));
    end
end
