import numpy as np
from math import sqrt #faster than numpy

import planetary_data as pd
import tools as t
from OrbitPropagator import OrbitPropagator as OP

#time parameters
tspan=10000 # seconds
dt = 1.0
cb=pd.mars

if __name__ == '__main__':
    r_mag=cb['radius']+400
    v_mag=np.sqrt(cb['mu']/r_mag) #circular orbit
    r0=np.array([r_mag,0,0])
    v0=np.array([0,v_mag,0])
    state0=np.concatenate((r0,v0),axis=None)
    
    op0=OP(state0,tspan,dt)
    
    op0.propagate_orbit()