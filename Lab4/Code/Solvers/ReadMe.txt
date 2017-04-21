Solvers folder has all the solvers:
forwardEuler.m: Uses forward euler method.

backwardEulerPC.m: Predicts using forward Euler and then calculates backward euler, like Predictor-Corrector.

backwardEulerSym.m: Calculated true backward euler for linear ODEs, where everything is taken to RHS and solved for x_(i+1).

backwardEulerNewton1.m: Uses Newton method for solving non-linear equation in backward euler formulation. 

trapezEulerPC.m, trapezEulerSym.m, trapezEulerNewton1.m: Refer to explanation for respective backward euler functions.

rk34.m: Runge-Kutta method, error estimation by 3-4. 

rk34Adaptive.m: Adaptive RK34, but only for the next step.

rk34Adaptive2.m: Adaptive RK34 for both forward and backward step. 
