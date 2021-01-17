import planetary_data as pd
import tools as t
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts

#Cd coefficient of drag, 2.2 spacecraft
#A area
if __name__ == '__main__':
    perts=null_perts()
    perts['aero']=True
    perts['Cd']=2.2
    perts['A']=(1e-3)**2/4.0 # km^2


    # initial mass of spacecraft
    mass0=10.0 # kg

    # perigee and apogee
    rp=215+cb['radius'] # km
    rp=300+cb['radius'] # km

    # orbital element angles
    raan=340.0
    i=65.2
    aop=58.0
    ta=332.0

    # calculate other orbital elements
    a=(rp+ra)/2.0 # km
    e=(ra-rp)/(ra+rp)

    # calculate initial state vector
    state0=[a,e,i,ta,aop.raan]

    op=OP