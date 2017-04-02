% This function calculates hessian of V using Quasi-Newton method
% ------------------------------------------------------------------------
% Pragya Sharma, March 30 2017
% ------------------------------------------------------------------------

function hessV = hessFuncvSC(a,sMeas,sModel,h)
    nPara = length(a);
    hessV = zeros(nPara);
    hessNumerator = (funcV(a+2*h,sMeas,sModel) - 2*funcV(a,sMeas,sModel)+...
                    funcV(a-2*h,sMeas,sModel));
    for i = 1:nPara
        for j = 1:nPara
            hessV(i,j) = (hessNumerator)/(4*h(i)*h(j));
        end
    end
   
end