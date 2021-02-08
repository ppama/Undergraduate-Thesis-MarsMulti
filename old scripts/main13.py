import planetary_data as pd
import tools as t_
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts

#Cd coefficient of drag, 2.2 spacecraft
#A area
cb = pd.mars

# time parameterss
tspan = 3600*6 # hour
dt=100.0

if __name__ == '__main__':
    #perturbations dictionary, arbitrary

    # initial mass of spacecraft
    mass0=50.0 # kg

    # perigee and apogee
    rp=300+cb['radius'] # km
    ra=300+cb['radius'] # km

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
    
    