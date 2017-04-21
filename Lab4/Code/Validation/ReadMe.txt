validateODE.m: main script
fValidate.m: Function for ode of circuit
dfValidate.m: Jacobian (here, just differentiation wrt one variable) for Newton method. Jacobian is calculated symbolically here, could be easily converted to numerical calculation.
xImplicit.m: Returns x_(i+1) on LHS (calculated symbolically).
fValidateSolution.m: True solution of ODE for checking solvers' accuracy.