% Pragya Sharma, March 26, 2017
% Task 2, Lab 4, ECE 4960
% Parameter extraction validation using Power Law
function [sol,aold,dela,v,iter,delv] = paraExtractTask2()

    % Define x_i to generate S_i and S_i,measured.
    nSample = 10;
    x = logspace(0,3,nSample)';
    
    % Generate measurement data
    c0 = 10;
    m = -0.5;
    sMeas = funcx(x,c0,m);
    
    % Starting quasi-newton here
    a0 = [0;0];
    nPara = length(a0);
    maxIter = 100;
    t = [10, 1, 0.1]; % Line-search parameter
    iter = 1;
    tol = 1e-9;
    v = zeros(1,maxIter);
    delv = zeros(2,maxIter);   % First differentiation
    dela = zeros(2,maxIter);
    aold = zeros(2,maxIter+1);
    aold(:,1) = a0;
    flagIter = 1;
    sol = 0;
    while (iter <= maxIter) && (flagIter == 1)
        v(iter) = funcV(x,aold(:,iter),sMeas);
        delv(:,iter) = delFuncv(x,aold(:,iter),sMeas);
        dela(:,iter) = -inv(hessFuncv(x,aold(:,iter),sMeas))*delv(:,iter);
        vLineSearch = zeros(length(t),1);
        anewLineSearch = zeros(nPara,length(t));
        normvLineSearch = zeros(1,length(t));
        for iter2 = 1:length(t)
            anewLineSearch(:,iter2) = aold(:,iter) + t(iter2)*dela(:,iter);
            vLineSearch(iter2) = funcV(x,anewLineSearch(:,iter2),sMeas);
            normvLineSearch(iter2) = norm(vLineSearch(iter2));
        end
        [~,I] = min(normvLineSearch);
        fprintf('t chosen %d\n',t(I));
        anew = anewLineSearch(:,I);
        aold(:,iter+1) = anew;
        if norm(dela(:,iter)) < tol
            flagIter = 0;
            sol = [10^anew(1); anew(2)];
            fprintf('The true parameters are \nc0 = %f, \nm = %f.\nThe solution is \nc0 = %f, \nm = %f.\n',c0,m,sol(1),sol(2));
        end
        iter = iter + 1;
    end
    iter = iter - 1;
    v = v(1:iter);
    delv = delv(:,1:iter);   % First differentiation
    dela = dela(:,1:iter);
    aold = aold(:,1:iter+1);
end

function sMeas = funcx(x,c0,m)
    c = log10(c0);
    nSample = length(x);
    % Element-wise operation, x can be a column vector.
    sMeas = c.* ones(nSample,1) + m.*x; 
end

function v = funcV(x,a,sMeas)
    nSample = length(x);
    v = 0;
    for i = 1:nSample
        v = v + (a(1) + a(2)*x(i) - sMeas(i))^2;
    end
end

function delv = delFuncv(x,a,sMeas)
    delv = zeros(2,1);
    hx1 = [0.0001;0];
    hx2 = [0;0.0001];
    delv(1) = (funcV(x,a+hx1,sMeas)-funcV(x,a,sMeas))/(hx1(1));
    delv(2) = (funcV(x,a+hx2,sMeas)-funcV(x,a,sMeas))/(hx2(2));
end

function hessV = hessFuncv(x,a,sMeas)
hessV = zeros(2,2);
hx1 = [0.001; 0];
hx2 = [0; 0.001];
% d^2V/da1^2
hessV(1,1) = -(-funcV(x,a+hx1,sMeas)+2*funcV(x,a,sMeas)-funcV(x,a-hx1,sMeas))/(hx1(1)^2);
% d^2V/da2^2
hessV(2,2) = -(-funcV(x,a+hx2,sMeas)+2*funcV(x,a,sMeas)-funcV(x,a-hx2,sMeas))/(hx2(2)^2);
% d^2V/da1da2
hessV(2,1) = (funcV(x,a+hx1+hx2,sMeas)-funcV(x,a+hx1-hx2,sMeas)-...
    funcV(x,a-hx1+hx2,sMeas)+funcV(x,a-hx1-hx2,sMeas))/(4*hx1(1)*hx2(2));
% d^2V/da2da1
hessV(1,2) = (funcV(x,a+hx1+hx2,sMeas)-funcV(x,a-hx1+hx2,sMeas)-...
    funcV(x,a+hx1-hx2,sMeas)+funcV(x,a-hx1-hx2,sMeas))/(4*hx1(1)*hx2(2));
end