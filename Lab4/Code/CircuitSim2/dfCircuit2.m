function jac = dfCircuit2(x,t,ode)
% DFCIRCUIT2 differentiates ode for Circuit 2 to get Jacobian
%   It takes input as x,t and function handle for the ode

lengthX = length(x);
jac = zeros(lengthX);
stepSize = 1e-3;

for i = 1:lengthX
    h = zeros(lengthX,1);
    h(i) = stepSize;
    jac(:,i) = (ode(x+h,t)-ode(x-h,t))/(2*stepSize);
end

end