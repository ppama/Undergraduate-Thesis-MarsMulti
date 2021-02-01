import numpy as np
import matplotlib.pyplot as plt
#plt.style.use('dark_background')
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D


import planetary_data as pd
import tools as t_
def null_perts(): #dictionary of perturbations
    return {
            'aero':False,
            'thrust':0,
            'orbit_direction':0
    }
class marsOrbitPropagator:
    def __init__(self,state0,tspan,dt,perts=null_perts(),cb=pd.mars,mass0=0,propagator='lsoda'):
        # state0: position x0=s0,y0=r0,  gamma0,theta0,vx,vy
        self.r0=state0[:2]
        self.v0=state0[4:]
        self.gamma=state0[2] # FPA [rad]
        self.theta=state0[3] # angle of orbit subtended
        
        
        self.tspan=tspan
        self.dt=dt
        self.cb=cb
        self.mass0=mass0
        # define perturbation dictionary
        self.perts=perts
        
        #total number of steps
        self.n_steps=int(np.ceil(self.tspan/self.dt))+1
        
        #initialize variables
        self.y=np.zeros((self.n_steps+1,6))#6 states, rx,ry,gamma,theta,vx,vy
        self.ts=np.zeros((self.n_steps+1,1))
        self.alts=np.zeros((self.n_steps+1))
        self.propagator=propagator
        self.step=0
        
        #initial conditions
        print(self.r0.tolist())
        print(self.v0.tolist())
        self.y[0,:]=self.r0.tolist()+[self.gamma]+[self.theta]+self.v0.tolist()
        self.alts[0]=t_.norm(self.r0)-self.cb['radius']
        
        #initiate solver
        self.solver=ode(self.diffy_q)
        self.solver.set_integrator(self.propagator)
        self.solver.set_initial_value(self.y[0,:],0)
        
        self.propagate_orbit()
        
    def propagate_orbit(self):
        print('Propagating orbit...')
        
        #propagate orbit, check for max time and stop conditions at each step
        while self.solver.successful() and self.step<self.n_steps: # and self.check_stop_conditions():
            #integrate step
            self.solver.integrate(self.solver.t+self.dt)
            self.step+=1
            
            #extract values from solver instance
            self.ts[self.step]=self.solver.t
            self.y[self.step]=self.solver.y
            
            # calcualte altitude at this time step
            self.alts[self.step]=t_.norm(self.solver.y[:2])-self.cb['radius']
        
        # extract arrays at the step where the propagation stopped
        self.ts=self.ts[:self.step]
        self.rs=self.y[:self.step,:2]
        self.vs=self.y[:self.step,4:5]
        self.gamma=self.y[:self.step,2]
        self.theta=self.y[:self.step,3]
        self.alts=self.alts[:self.step]
    
    def diffy_q(self,t,y):
        # unpack state
        rx,ry,gamma,theta,vx,vy=y
        r=np.array([rx,ry])
        v=np.array([vx,vy])
        # norm of the radius vector
        norm_r=np.linalg.norm(r)
        norm_v=np.linalg.norm(v)
        
        # two body acceleration
        a=-r*self.cb['mu']/norm_r**3*np.sin(gamma)
        
        # calculate angles
        dgamma=(self.cb['mu']/(norm_v*norm_r**2)-norm_v/norm_r)*np.cos(gamma)
        dtheta=norm_v*np.cos(gamma)/norm_r
        # calculate radius of orbit
        r=norm_v*np.sin(gamma)
        
        
        # calculate altitude and air density
        h=norm_r-self.cb['radius']
        rho=t_.calc_atmospheric_density(h)
        
        # calculate motion of spacecraft w/ respect to a rotating atmosphere
        v_rel=v-np.cross(self.cb['atm_rot_vector'],r)
        # calculate drag
        drag=-v_rel*0.5*rho*(t_.norm(v_rel))**2*self.perts['Cd']*self.perts['A']/self.mass0
        a+=drag
        
        return [vx,vy,dgamma,dtheta,a[0],a[1]]