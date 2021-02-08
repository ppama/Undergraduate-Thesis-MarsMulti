import numpy as np
import tools as t_
from orbprop3d import OrbitPropagator as OP
import matplotlib.pyplot as plt
import sys
import json
import os



if __name__ == '__main__':
    
    if len(sys.argv) != 2:
        print('usage: %s parameter_file' % sys.argv[0])
        exit(1)

    with open(sys.argv[1], 'r') as f:
        params = json.load(f)

    # loading input parameters
    craft = params['craft']
    planet = params['planet']
    sim = params['sim']
    
    atm_data = np.loadtxt(planet['atm_data']) # [altitude,temp,pressure,density]
    fpa = craft['fpa']
    
    r_mag = planet['radius']+craft['entry_altitude']
    v_mag = craft['velocity']
    r0 = np.array([0,r_mag,0])
    v0 = np.array([v_mag*np.cos(fpa),v_mag*np.sin(fpa),0])
    state0 = np.concatenate((r0,v0),axis=0)
    
    
    op=OP(state0,planet,craft,sim)
    #maxrange=op.rs[op.step-1,0]
    
    op.plot_state(show_plot=True)
    op.plot_3d(show_plot=True)
    
