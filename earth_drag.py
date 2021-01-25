import planetary_data as pd
import tools as t
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts

#Cd coefficient of drag, 2.2 spacecraft
#A area

cb = pd.earth

# time params
tspan = 3600*24*3
dt=100.0

if __name__ == '__main__':
    #perturbations dictionary, arbitrary
    
    perts=null_perts()
    perts['aero']=True
    perts['Cd']=2.2
    perts['A']=(1e-3)**2/4.0 # km^2

    # initial mass of spacecraft
    mass0=10.0 # kg

    # perigee and apogee
    rp=215+cb['radius'] # km
    ra=300+cb['radius'] # km

    # orbital element angles
    raan=340.0
    i=65.2
    aop=58.0
    ta=332.0

    # calculate other orbital elements
    a=(rp+ra)/2.0 # km
    e=(ra-rp)/(ra+rp)

    # calculate initial state vector
    state0=[a,e,i,ta,aop,raan]

    op=OP(state0,tspan,dt,deg=True,coes=True,mass0=mass0,perts=perts)
    op.plot_alts(show_plot=True,hours=True)
    op.plot_3d(show_plot=True)
    op.calculate_coes()
    op.plot_coes(show_plot=True,hours=True)
    op.calculate_apoapse_periapse()
    op.plot_apoapse_periapse(show_plot=True,hours=True)
    
    