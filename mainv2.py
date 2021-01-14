import planetary_data as pd
import tools as t
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts

#time parameters
tspan=3600*24*1.0 # seconds
dt = 10.0

# central body
cb=pd.earth

if __name__ == '__main__':
    perts=null_perts()
    perts['J2']=True
    #,e,i,ta,aop,raan=coes
    
    # ISS
    c0=[cb['radius']+414.0,0.0000291,51.6425,0.0,285.3988,22.3967]
    
    op=OP(c0,tspan,dt,coes=True,perts=perts)
    op.propagate_orbit()
    op.plot_3d(show_plot=True)