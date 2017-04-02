% This function formulates V that describes the least square fitting.
% ------------------------------------------------------------------------
% Pragya Sharma, March 30 2017
% ------------------------------------------------------------------------

function v = funcV(a,sMeas,sModel)
    [nSample,~] = size(sMeas);
    v = 0;
    sMod = sModel(a,sMeas);
    for i = 1:nSample
        v = v + (sMod(i) - sMeas(i,end))^2;
    end
end
