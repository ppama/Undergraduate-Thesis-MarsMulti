import planetary_data as pd
import numpy as np
import tools2 as t_
from marsOrbitPropagator import marsOrbitPropagator as OP
from marsOrbitPropagator import null_perts
import matplotlib.pyplot as plt

# time parameterss
tspan = 10000 # hour
dt=0.1
cb = pd.mars
#Physical Parameters
Cd=2.2
A=(1e-3)**2/4.0 # km^2

if __name__ == '__main__':
    
    perts=null_perts()
    perts['Cd']=2.2
    perts['A']=(1e-3)**2/4.0 # km^2
    
    center2interface=(cb['atmos_interface']+cb['radius'])*1000
    entry_apo_alt=500
    entry_peri_alt=50
    mass0=100 # kg
    
    # get velocity and flight path angle at entry interface given entry orbit
    #v_interface,gamma_interface=t_.planetEntryConditions(entry_apo_alt,entry_peri_alt)
    #vx_interface=v_interface*np.cos(gamma_interface)
    #vy_interface=v_interface*np.sin(gamma_interface)
    # v=omega*r, atmospheric velocity at interface (m/s) 
    #vrot_interface=2*np.pi*cb['radius']/(cb['period']*3600) 
    #vx_init=vx_interface+vrot_interface
    #vy_init=vy_interface
    
    # Entry velocity (m/s) (relative to planet frame)
    v_init=5600 # m/s
    # Entry FPA [rad] relative to planet frame 
    gamma_init=6.5 # deg
    vx_init=v_init*np.cos(gamma_init) # m/s
    vy_init=v_init*np.sin(gamma_init)
    h_init=125*1000 # entry altitude
    y_init=h_init+cb['radius']*1000
    s_init=0 # entry range [m]
    theta_init=s_init*np.cos(gamma_init)/y_init
    
    # state0: rx, ry, gamma, theta, vx,vy
    state0=[s_init,y_init,gamma_init,theta_init,vx_init,vy_init]
    
    op=OP(state0,tspan,dt,mass0=mass0,perts=perts)
    
    #t_.atm_plot()
    
    