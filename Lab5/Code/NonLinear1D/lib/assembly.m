function fNew = assembly(f,n,he,testOn,tol)
%ASSEMBLY function 
nn = n + 1;                 % Number of nodes
k = zeros(nn);              % Initializing main LHS matrix  
r = zeros(nn,1);            % Initializing RHS vector
ds = [-1/he,  1/he];          
lmn = zeros(n,2);           % Connectivity matrix of nodes
for i = 1:n
    lmn(i,:) = [i,i+1];
end

y1 = @func1;
y2 = @func2;
y3 = @func3;
y4 = @func4;

for i = 1:n
    lm = lmn(i,:);
    k11 = -(ds'*ds)*he + (simpInt(y3,0,he,he,i)*f(lm(1))...
        *simpInt(y2,0,he,he,i)*f(lm(2))) - simpInt(y1,0,he,he,i);
    
    % --------------------------------------------------------------------
    if (testOn)
        syms xSym                     
        sSym = [1 - xSym/he, xSym/he]; 
        dsSym = diff(sSym,xSym);
        k11True = -int(dsSym'*dsSym,xSym,0,he) + (int(sSym'*dsSym*...
        sSym(1),xSym,0,he)*f(lm(1))*int(sSym'*dsSym*sSym(2),xSym,0,he)...
        *f(lm(2))) - int(sSym'*sSym,xSym,0,he);
        if(norm(double(k11True)-k11)<tol)
            fprintf('Integration module works correctly.\n');
        else
            fprintf('Integration module doesn''t converge to true solution.\n');
        end
        testOn = 0;
    end       

    % --------------------------------------------------------------------
    f1 = simpInt(y4,0,he,he,i);
    k(lm,lm) = k(lm,lm) + k11;
    r(lm) = r(lm) + f1;
end

% Boundary conditions
k(1,:) = 0;
k(nn,:) = 0;
k(1,1) = 1;
k(nn,nn) = 1;
r(1,1) = f(1);
r(nn,1) = f(nn);

% Final calculations
d = k\r;
fNew = d;
end