import planetary_data as pd
import numpy as np
import tools2 as t_
from orbprop3d import OrbitPropagator as OP
import matplotlib.pyplot as plt
import baseline as ba

# time parameterss
tspan = 10000 # seconds
dt=0.1
cb = pd.mars
sim = ba.sim
craft = ba.craft
#Physical Parameters
Cd=2.2
A=(1e-3)**2/4.0 # km^2
fpa = sim['fpa']

if __name__ == '__main__':
    
    r_mag = cb['radius']+sim['entry_altitude']
    v_mag = sim['velocity']
    r0 = np.array([0,r_mag,0])
    v0 = np.array([v_mag*np.cos(fpa),v_mag*np.sin(fpa),0])
    state0 = np.concatenate((r0,v0),axis=0)
    
    op=OP(state0,tspan,dt)
    op.plot_state(show_plot=True)
    #op.plot_3d(show_plot=True)
    
