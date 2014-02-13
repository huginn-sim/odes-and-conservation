# -*- coding: utf-8 -*-
"""
.. module:: sho
   :synopsis: Describes both discrete and analytic models that 

.. moduleauthor:: Huginn
"""
#~ Modules
import  sys
from pylab import *
#/~ Modules

#~ Custom Modules
sys.path.append("C:\Users\Evin\Documents\Github")
from viz.display.plot import configure
#/~ Custom Modules

#~ ODE Solvers
def euler(f,t,dt,x):
    return x + f(t,x)*dt

def euler_richardson(f,t,dt,x):
    return x + f(t + dt/2., x + f(t,x)*dt/2.)*dt

def rk4(f,t,dt,x):
    k1 = f(t        , x        )*dt
    k2 = f(t + dt/2., x + k1/2.)*dt
    k3 = f(t + dt/2., x + k2/2.)*dt
    k4 = f(t + dt   , x + k3   )*dt

    return x + (1./6.)*(k1 + 2*k2 + 2*k3 + k4)

def predict_correct(f,t,dt,x):
    t[0:0] = [t[0] - dt]
    # Roll back by one time-step.
    pc_state = [rk4(f,t[0],-dt,x[0])]
    pc_state.append(x[0])
    # Roll forward.
    for t in times:
        xp = pc_state[-2] + f(t,pc_state[-1])*2*dt
        xc = pc_state[-1] + 0.5*(f(t,xp) + f(t,pc_state[-1]))*dt
        pc_state.append(xc)

    return times, array(pc_state)
#/~ ODE Solvers

#~ ODEs
an = lambda t,xn: -k/m*xn
discrete_spring_ode = lambda t,x: array([x[1], an(t,x[0]), -k/m])
analytic_spring_ode = lambda t,x: array([x[0]*cos(sqrt(k/m)*t), -x[0]*sqrt(k/m)*sin(sqrt(k/m)*t)]).T

spring_energy = lambda t,x: .5*k*x[:,0]**2 + .5*m*x[:,1]**2
#/~ ODEs

#~ Plot Functions
def plot_pvt():
    astate = analytic_spring_ode(times, state[-1])

    dstate = list(state)
    for t in times[1:]:
        dstate.append(rk4(discrete_spring_ode, t, dt, dstate[-1]))
    dstate = array(dstate)

    agreement = astate[:,0] - dstate[:,0]

    fig, ax = subplots()
    agmark, = ax.plot(times, agreement, 'k--', lw=1, alpha=.7)
    amark, = ax.plot(times, astate[:,0], 'b-', lw=6)
    dmark, = ax.plot(times, dstate[:,0], 'r-', lw=2)

    ax.legend(  [amark, dmark, agmark],
                [r'Analytic $\left(\mathcal{A}_t\right)$', r'Discrete $\left(\mathcal{D}_t\right)$', r'Agreement $\left(\mathcal{A}_t-\mathcal{D}_t\right)$'],
                numpoints=1,
                loc="lower right")

    configure(  ax=ax,
                title=r"$k=1.0$ $\frac{1}{seconds^2}$; $m=1.0$ $gram(s)$; $v_0=0.0$ $\frac{meter(s)}{second}$",
                xlabel="Time (seconds)\n$\\Delta t=" + str(dt) + "$",
                ylabel="Spatial Displacement (meters)\n$x_0="+str(x0)+"$",
                xbounds=(t0, tf), ybounds=(-1.,1.))   

    fig.suptitle("Simple Harmonic Oscillator", size=30)
    fig.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.08)

def plot_evt(dt, integral, suptitle='Simple Harmonic Oscillator'):
    # Do this once.
    a_state = analytic_spring_ode(times, state[-1])
    a_energy = spring_energy(times, a_state)
    
    d_state = list(state)
    if integral != predict_correct:
        for t in times[1:]:
            d_state.append(integral(discrete_spring_ode, t, dt, d_state[-1]))
    else:
        d_times, d_state = predict_correct(discrete_spring_ode, list(times), dt, d_state)

    d_state = array(d_state)
    d_energy = spring_energy(times, d_state)

    fig, ax = subplots()
    
    amark, = ax.plot(times, a_energy, 'b-', lw=6)

    if integral == predict_correct:
        print len(times); print len(d_energy)
        dmark, = ax.plot(d_times, d_energy[:-2], 'r-', lw=2)
    else:
        print len(times); print len(d_energy)
        dmark, = ax.plot(times, d_energy, 'r-', lw=2)

    ax.legend(  [amark, dmark],
                [r'Analytic $\left(\mathcal{A}_t\right)$', r'Discrete $\left(\mathcal{D}_t\right)$'],
                numpoints=1,
                loc="lower right")

    configure(  ax=ax,
                title=r"$k=1.0$ $\frac{1}{seconds^2}$; $m=1.0$ $gram(s)$; $v_0=0.0$ $\frac{meter(s)}{second}$",
                xlabel="Time (seconds)\n$\\Delta t=" + str(dt) + "$",
                ylabel="Energy (joules)\n$E_0="+str(E0)+"$",
                xbounds=None, ybounds=(0.49,0.51))   

    error = abs(E0 - d_energy[-1]) / E0
    print error
    fig.suptitle(suptitle, size=30)
    fig.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.08)

def plot_pc():
    a_state = analytic_spring_ode(times, state[-1])
    pc_times, pc_state = predict_correct(discrete_spring_ode,times,dt,state)

    fig, ax = subplots()
    amark, = ax.plot(times, a_state[:,0], 'b-', lw=6)
    pcmark, = ax.plot(pc_times, pc_state[1:-1,0], 'r-', lw=2)

    ax.legend(  [amark, pcmark],
                [r'Analytic $\left(\mathcal{A}_t\right)$', r'Predictor-Corrector $\left(\mathcal{D}_t\right)$'],
                numpoints=1,
                loc="lower right")

    configure(  ax=ax,
                title=r"$k=1.0$ $\frac{1}{seconds^2}$; $m=1.0$ $gram(s)$; $v_0=0.0$ $\frac{meter(s)}{second}$",
                xlabel="Time (seconds)\n$\\Delta t=" + str(dt) + "$",
                ylabel="Position (meters)\n$x_0="+str(x0)+"$",
                xbounds=(t0,tf), ybounds=None)   

    fig.suptitle("Predictor Corrector", size=30)
    fig.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.08)
#/~ Plot Functions

#~ Error Analysis
def find_error(interval=(0,2), step=.001, accuracy=.01):
    a_state = analytic_spring_ode(times, state[-1])
    a_energy = spring_energy(times, a_state)

    integrals = [euler, euler_richardson, rk4, predict_correct]
    errors = {}
    dts = arange(interval[0], interval[1], step)[::-1]
    for dt in dts:
        for integral in integrals:
            if integral.__name__ not in errors:
                errors[integral.__name__] = None

            if errors[integral.__name__] == None:
                d_state = list(state)
                if integral != predict_correct:
                    for t in times:
                        if integral != predict_correct:
                            d_state.append(integral(discrete_spring_ode, t, dt, d_state[-1]))
                
                else:
                    d_state = integral(discrete_spring_ode, list(times), dt, d_state)[1]
                
                d_state = array(d_state)
                d_energy = spring_energy(times, d_state)

                error = abs(E0 - d_energy[-1]) / E0

                if error <= accuracy:
                    errors[integral.__name__] = (error, dt)
                    break

    return errors
#/~

#~ Global Variables
t0 = 0; tf = 10*pi; dt = .1
times = np.arange(t0, tf, dt)

k = 1.; m = 1.;
x0 = 1.; v0 = 0.; a0 = an(t0,x0); E0 = 0.5;
state = [np.array([x0,v0,a0])]
#/~ Global Variables

if __name__ == "__main__":
    plot_evt(dt=.005, integral=euler, suptitle="Euler Method")
    plot_evt(dt=.106, integral=euler_richardson, suptitle="Euler-Richardson")
    plot_evt(dt=.3634, integral=rk4, suptitle="Runge Kutta (RK4) Method")
    plot_evt(dt=.089, integral=predict_correct, suptitle="Predictor Corrector")

    show()
    sys.exit()
# Results of Error Analysis ####################################
    #errors = find_error()
    rk4_dt = .3634
    rk4_e = .009958

    er_dt = .106
    er_e = .009991

    pc_dt = .089
    pc_e = .009639

    e_dt = .005
    e_e = .007906