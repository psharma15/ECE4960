% This function calculates hessian of V using Quasi-Newton method
% ------------------------------------------------------------------------
% Pragya Sharma, March 30 2017
% ------------------------------------------------------------------------

function hessV = hessFuncvQN(a,sMeas,sModel,h)
    nPara = length(a);
    
    hessV = zeros(nPara);

    for i = 1:nPara
        hx = zeros(nPara,1);
        hx(i) = h(i);
        for j = 1:nPara
            hy = zeros(nPara,1);
            hy(j) = h(j);
            hessV(i,j) = (funcV(a+hx+hy,sMeas,sModel) - funcV(a-hx+hy,sMeas,sModel)-...
                funcV(a+hx-hy,sMeas,sModel) + funcV(a-hx-hy,sMeas,sModel))/(4*hx(i)*hy(j));
        end
    end
end