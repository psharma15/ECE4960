circuit1ODE.m: main script
fCircuit1.m: Function for ode of circuit
dfCircuit.m: Jacobian for Newton method. Jacobian is calculated symbolically here, could be easily converted to numerical calculation.
xImplicit.m: Returns x_(i+1) on LHS (calculated symbolically).
plotCurrent.m: Plotting the current generated, needed as input.