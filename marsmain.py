import numpy as np
import tools as t_
from orbprop3d import OrbitPropagator as OP
import matplotlib.pyplot as plt
import sys
import json
import os

class craftParams(object):
        def __init__(self,fpa,entry_altitude,beta,v_init):
            self.fpa = craft['fpa']
            self.entry_altitude = craft['entry_altitude']
            self.beta = craft['ballistic_coef']
            self.v_init = craft['velocity']


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
    zs=atm_data[:,0], # km
    temp=atm_data[:,1], # atmospheric temps [K]
    p=atm_data[:,2], #atmospheric pressure [mb]=100Pa
    rhos=atm_data[:,3], # density [g/cm^3]
    fpa = np.radians(craft['fpa'])
    entry_altitude = craft['entry_altitude'] 
    beta = craft['ballistic_coef']
    r_mag = planet['radius']+craft['entry_altitude']
    v_init = craft['velocity'] 
    #r0 = np.array([0,r_mag,0])
    #v0 = np.array([v_init*np.cos(fpa),v_init*np.sin(fpa),0])
    #state0 = np.concatenate((r0,v0),axis=0)
    
    finalRange = [[],[]]
    fpaList = finalRange[0]
    rangeList = finalRange[1]
    escape = False
    fpa = np.radians(-20)
    op=OP(planet,craft,sim,fpa,entry_altitude,beta,v_init)
    #while not escape:
    #    op=OP(planet,craft,sim,fpa,entry_altitude,beta,v_init)
        
    #    fpa = fpa + 1
    #    fpaList.append(fpa)
    #    rangeList.append(op.range[op.step-1])
    #    escape = op.escape
        
    #print(finalRange)
        
        
    #maxrange=op.rs[op.step-1,0]
    
    op.plot_state(show_plot=True)
    #op2.plot_state(show_plot=True)
    #op.plot_3d(show_plot=True)
    
