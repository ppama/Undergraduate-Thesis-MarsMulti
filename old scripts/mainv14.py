import planetary_data as pd
import tools as t_
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts

# time parameterss
tspan = 3600*6 # hour
dt=100.0

cb = pd.earth

if __name__ == '__main__':
    #perturbations dictionary, arbitrary
    
    perts=null_perts()
    perts['thrust']=16.627
    perts['thrust_direction']=-1 # 1 prograde -1 retrograde
    perts['isp']=4300
    perts['aero']=True
    perts['Cd']=2.2
    perts['A']=(1e-3)**2/4.0 # km^2
    
    # define stop condition dictionary
    sc={'min_alt':0.0}
    
    mass0=50.0 # kg
    
    # calculate initial state vector
    state0=[cb['radius']+800,0.03,10.0,0.0,0.0,0.0]
    
    op=OP(state0,tspan,dt,deg=True,coes=True,mass0=mass0,perts=perts,sc=sc)
    op.plot_alts(show_plot=True,hours=True)
    op.plot_3d(show_plot=True)
    op.calculate_coes()
    op.plot_coes(show_plot=True,hours=True)
    op.calculate_apoapse_periapse()
    op.plot_apoapse_periapse(show_plot=True,hours=True)