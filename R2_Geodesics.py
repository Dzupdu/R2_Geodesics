
from sympy import Eq,solve_linear_system, Matrix
from numpy import linalg
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.integrate import odeint


from sympy.diffgeom.rn import R2
from sympy.diffgeom import metric_to_Christoffel_1st, metric_to_Christoffel_2nd, TensorProduct
TP = TensorProduct

def metric_to_geodesic_eq(g):
    Gammas2 = metric_to_Christoffel_2nd(g)
    
    t = sp.symbols('t') 
    x, y = sp.symbols('x y', cls=sp.Function) 
    x=x(t)
    y=y(t)
    xys = [x,y]
    # simplify, for some reason this is not default
    Coefs = [[[(Gammas2[k,j,i]).simplify().subs([(R2.x,x), (R2.y,y)]) for i in range(2)] for j in range(2)] for k in range(2)]
    eqx=sp.diff(x,t,t)
    eqy=sp.diff(y,t,t)

    for i in range(2):
        for j in range(2):
            eqx = eqx + Coefs[0][i][j]*sp.diff(xys[i],t)*sp.diff(xys[j],t)
            eqy = eqy + Coefs[1][i][j]*sp.diff(xys[i],t)*sp.diff(xys[j],t)
            
    return [eqx, eqy]

#def sym2_to_numeric4_eq(eqs):
    

def solve_geodesic_eq(eqs,S0, tmax):
    dxdt, dydt = sp.symbols('dxdt dydt', cls=sp.Function)
    t = sp.symbols('t') 
    x, y = sp.symbols('x y', cls=sp.Function) 
    x=x(t)
    y=y(t)
    dxdt = dxdt(t)
    dydt = dydt(t)
    
    eqx = eqs[0]
    eqy = eqs[1]
    
    eqx =eqx.subs([(sp.diff(x,t),dxdt),(sp.diff(y,t),dydt)])
    eqy =eqy.subs([(sp.diff(x,t),dxdt),(sp.diff(y,t),dydt)])

    numeqs = sp.Matrix([dxdt,dydt, sp.solve(eqx,sp.diff(dxdt,t))[0], sp.solve(eqy,sp.diff(dydt,t))[0]])
    
    # cant make this work without repetition for some reason
    #df = sp.lambdify([t,x,y,dxdt,dydt], numeqs)
    #def dFs(S,t2):
    #    return df(t2, S[0],S[1],S[2],S[3])

    dS1 = sp.lambdify([t,x,y,dxdt,dydt], numeqs[0])
    dS2 = sp.lambdify([t,x,y,dxdt,dydt], numeqs[1])
    dS3 = sp.lambdify([t,x,y,dxdt,dydt], numeqs[2])
    dS4 = sp.lambdify([t,x,y,dxdt,dydt], numeqs[3])

    def dFs(S,t2):
        return [dS1(t2, S[0],S[1],S[2],S[3]),dS2(t2, S[0],S[1],S[2],S[3]),dS3(t2, S[0],S[1],S[2],S[3]),dS4(t2, S[0],S[1],S[2],S[3])]
    
    t_range =np.linspace(0,tmax,round(tmax*1000))
    
    sol = odeint(dFs, y0=S0, t=t_range, full_output =0)
    
    return sol.T

def plot_geodesics(eqs, x0, n_angles = 8,xlim = 1,ylim = 1, tmax = 1):

    fig, ax = plt.subplots(figsize=(6, 6))
    angles = np.linspace(0,2*np.pi,n_angles + 1)
    
    for a in angles[:-1]:
        S0=[x0[0],x0[1],np.cos(a),np.sin(a)]
        
        sol = solve_geodesic_eq(eqs,S0, tmax=tmax)
        ax.plot(sol[0], sol[1])

    ax.set(xlim=(-xlim, xlim), ylim=(-ylim, ylim))
    plt.show()