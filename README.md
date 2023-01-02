# R2_Geodesics

Solve and plot geodesics of given 2-dimensional metric.

Problems:
* Numpy basic ODE solve does not work reliably for complex metrics or some initial conditions
* Division by 0 or other exceptions due to initial conditions are not handled
* For some reason using numpy odeint solver chrashes the kernel from time to time, though this might depend on the hardware
* Can be used for solving geodesics on sphere or torus (or other 2d manifolds), but the visualisation is not very helpful for these cases
 
