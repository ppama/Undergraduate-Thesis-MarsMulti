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
            
def LandingPrecision(planet,craft,sim,state0,beta,fpa,deltav_mag,test_type):
        finalRange = [] #[fpa],[range]
        count = 0
        if test_type == 'fpa':
            while fpa <=np.radians(0):
                op=OP(planet,craft,sim,state0,beta,fpa)
                if op.escape:
                    break
                finalRange.append([np.rad2deg(fpa),np.abs(op.range[op.step-1][0]-baseRange)]) # FPA in degrees, range precision as absolute difference between range and base range
                fpa = fpa + np.radians(0.1)
                v0 = np.array([(v_mag+deltav_mag)*np.cos(fpa),(v_mag+deltav_mag)*np.sin(fpa),0])
                r0 = np.array([0,r_mag,0])
                state0 = np.concatenate((r0,v0),axis=0)
                count+=1
            print(f"done in {count} iterations.")
            
        elif test_type == 'deltav':
            while deltav_mag <= -0.01:
                op=OP(planet,craft,sim,state0,beta,fpa)
                if op.escape:
                    break
                finalRange.append([deltav_mag,np.abs(op.range[op.step-1][0]-baseRange)]) # deltav, range precision as absolute difference between range and base range
                deltav_mag += 0.0001
                v0 = np.array([(v_mag+deltav_mag)*np.cos(fpa),(v_mag+deltav_mag)*np.sin(fpa),0])
                r0 = np.array([0,r_mag,0])
                state0 = np.concatenate((r0,v0),axis=0)
                count+=1
            print(f"done in {count} iterations.")
            op.plot_state(show_plot=True)
            op.plot_3d(show_plot=True)
        
        elif test_type == 'beta':
            while beta < 100:
                op=OP(planet,craft,sim,state0,beta,fpa)
                if op.escape:
                    break
                finalRange.append([beta,np.abs(op.range[op.step-1][0]-baseRange)]) # beta, range precision as absolute difference between range and base range
                beta +=0.1
                count+=1
            print(f"done in {count} iterations.")
            op.plot_state(show_plot=True)
            op.plot_3d(show_plot=True)
        return(finalRange)

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
    
    #Simulate Baseline profile
    fpa = np.radians(craft['fpa']) # entry path angle
    entry_altitude = craft['entry_altitude'] 
    beta = craft['ballistic_coef'] * 10**6 # ballistic coeff. [kg/m^2] -> [kg/km^2]
    r_mag = planet['radius']+ entry_altitude
    v_mag = np.sqrt(planet['mu']/(planet['radius']+entry_altitude))
    deltav_mag = -0.01
    fpa = np.radians(0)
    v0 = np.array([(v_mag+deltav_mag)*np.cos(fpa),(v_mag+deltav_mag)*np.sin(fpa),0])
    r0 = np.array([0,r_mag,0])
    state0 = np.concatenate((r0,v0),axis=0)
    baseline=OP(planet,craft,sim,state0,beta,fpa)
    baseRange=baseline.range[baseline.step-1]
    baseline.plot_state(show_plot=True)
    baseline.plot_3d(show_plot=True)
    
    beta = 0.1 * 10**6
    v0 = np.array([(v_mag+deltav_mag)*np.cos(fpa),(v_mag+deltav_mag)*np.sin(fpa),0])
    r0 = np.array([0,r_mag,0])
    state0 = np.concatenate((r0,v0),axis=0)
    #op1=OP(planet,craft,sim,state0,beta,fpa)
    
    beta = 10
    v0 = np.array([(v_mag+deltav_mag)*np.cos(fpa),(v_mag+deltav_mag)*np.sin(fpa),0])
    r0 = np.array([0,r_mag,0])
    state0 = np.concatenate((r0,v0),axis=0)
    #op2=OP(planet,craft,sim,state0,beta,fpa)
    
    beta = 100
    v0 = np.array([(v_mag+deltav_mag)*np.cos(fpa),(v_mag+deltav_mag)*np.sin(fpa),0])
    r0 = np.array([0,r_mag,0])
    state0 = np.concatenate((r0,v0),axis=0)
    #op3=OP(planet,craft,sim,state0,beta,fpa)
    
    #t_.plot_n_orbits([op1.rs,baseline.rs,op2.rs,op3.rs],xrange=[op1.range[op1.step-1],baseline.range[baseline.step-1],op2.range[op2.step-1]],cb=planet,labels=['beta=0.1','beta=1','beta=10','beta=100'],show_plot=True)
    
    #fpa = np.radians(-90)
    #v0 = np.array([(v_mag+deltav_mag)*np.cos(fpa),(v_mag+deltav_mag)*np.sin(fpa),0])
    #r0 = np.array([0,r_mag,0])
    #state0 = np.concatenate((r0,v0),axis=0)
    #op=OP(planet,craft,sim,state0,beta,fpa)
    
    #deltav_range=LandingPrecision(planet,craft,sim,state0,beta,fpa,deltav_mag,test_type='deltav')
    #fpa_range=LandingPrecision(planet,craft,sim,state0,beta,fpa,deltav_mag,test_type='fpa')
    #beta_range=LandingPrecision(planet,craft,sim,state0,beta,fpa,deltav_mag,test_type='beta')
    
    #f=open("deltavrange.csv","w")
    #for elem in deltav_range:
    #    f.write(f"{elem[0]},{elem[1][0]}\n") #2nd element saves as 1 element list so make address [1][0] ¯\_(ツ)_/¯
    #f.close()
    
    #t_.atm_plot()
    
    
    
    
    