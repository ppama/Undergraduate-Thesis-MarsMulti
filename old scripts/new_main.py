import numpy as np
from math import sqrt #faster than numpy

import planetary_data as pd
import tools as t
from OrbitPropagator import OrbitPropagator as OP

#time parameters
tspan=3600*24*1.0 # seconds
dt = 10.0
cb=pd.mars

if __name__ == '__main__':
    r_mag=cb['radius']+400
    v_mag=np.sqrt(cb['mu']/r_mag) #circular orbit
    r0=np.array([r_mag,0,0])
    v0=np.array([0,v_mag,0])
    
    r_mag=cb['radius']+1000
    v_mag=np.sqrt(cb['mu']/r_mag)*1.3 #circular orbit
    r00=np.array([r_mag,0,0.6])
    v00=np.array([0,v_mag,0.3])
    
    op0=OP(r0,v0,tspan,dt)
    op00=OP(r00,v00,tspan,dt)
    
    op0.propagate_orbit(coes=False)
    op00.propagate_orbit()
    
    t.plot_n_orbits([op0.rs,op00.rs],labels=['orbit1','orbit2'],show_plot=True)