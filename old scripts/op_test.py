import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
plt.style.use('dark_background')

from sys import path
path.append('C:\\Users\\pasca\\OneDrive\\Documents\\GitHub\\orbitalmechanics\\python_tools')
from OrbitPropagator import OrbitPropagator as OP
import planetary_data as pd
cb=pd.mars

if __name__ == '__main__':
    #initial conditions of orbit parameters
    r_mag = cb['radius']+100.0 # km mag -> magnitude
    v_mag = np.sqrt(cb['mu']/r_mag)-1 # km/s
    
    #initial position and velocity vectors
    #r0=np.array([r_mag,0,0])
    #v0=np.array([0,v_mag,0])
    r0=[r_mag,r_mag*0.01,r_mag*-0.7]
    v0=[0,v_mag,v_mag*0.1]
    
    #timespan how long is simulation
    tspan= 6*3600*24.0 # s
    
    #timestep
    dt = 100.0
    
    op=OP(r0,v0,tspan,dt,cb=cb)
    op.propagate_orbit()
    op.plot_3d(show_plot=True)
    