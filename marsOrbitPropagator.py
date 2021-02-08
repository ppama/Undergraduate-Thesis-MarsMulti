import numpy as np
import matplotlib.pyplot as plt
#plt.style.use('dark_background')
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
import baseline as ba
import planetary_data as pd
import tools as t_
from rungekutta import RK4

class marsOrbitPropagator:
    """ Calculates orbital trajectory for a craft reentering the Martian atmosphere
    Assumes circular orbit, with origin at the planet center, uniform surface terrain, 
    ODEs used:
    gravitational acceleration = (mu/r^2)*(r/r_mag)-> towards planet center, surface 
    drag = 0.5*rho(z)*v^2/beta -> opposite velocity vector
    dv/dt = g - D
    dgamma/dt = ((mu/v*r^2)-(v/r))*cos(gamma) 
    dr/dt=v*singamma
    
    """
    
    def __init__(self,cb=pd.mars,sim=ba.sim,mass0=0,craft=ba.craft,propagator='lsoda'):
        self.cb = cb
        self.sim = sim
        self.tspan = sim['timespan']
        
        self.beta = craft['ballistic_coef']
        self.mass = craft['mass']
        self.n_steps=0
    
    def sim_run(self,cb,sim,craft,tspan,beta,mass):
        self.r=np.array([0,self.sim['entry_altitude'] + self.cb['radius']])
        self.x = np.zeros(self.tspan)
        self.y = np.zeros(self.tspan)
        self.gamma = np.zeros(self.tspan)
        self.gamma[0] = np.radians(sim['fpa'])
        self.theta = np.zeros(self.tspan) # angle of orbit subtended
        
        self.v = sim['velocity'] * np.array([np.cos(self.gamma[0]),np.sin(self.gamma[0])])
        self.vx = np.zeros(tspan)
        self.vy = np.zeros(tspan)
        
        self.a = np.array([0,0])
        self.ax = np.zeros(tspan)
        self.ay = np.zeros(tspan)
        
        self.dt = sim['delta_t']
        self.t = np.arange(0,tspan*self.dt,self.dt)
        
        def g_acc
        #total number of steps
        self.n_steps=int(np.ceil(self.tspan/self.dt))+1
        k = 0
        for _ in range(0,self.n_steps):
            self.r = self.r + self.v*self.dt
            self.x[k],self.y[k] = self.r
            
            # norm of the position and velocity vectors
            r_mag=np.linalg.norm(self.r)
            v_mag=np.linalg.norm(self.v)
            
            # calculate altitude and air density
            h=r_mag-self.cb['radius'] #center to craft - planet radius = altitude
            rho=t_.calc_atmospheric_density(h)
            
            # gravitational acceleration
            g_acc=-(self.cb['mu']/r_mag**2)*(self.r/r_mag) # negative because towards center, 
            # calculate motion of spacecraft w/ respect to a rotating atmosphere
            v_rel = np.subtract(self.v,[2*np.pi*r_mag*self.cb['atm_rot'],0]) # velocity relative to planet surface
            # calculate drag
            drag = - 0.5*rho*(v_rel)*np.linalg.norm(v_rel)/self.beta
            
            # calculate angles
            dgamma=(self.cb['mu']/(v_mag*r_mag**2)-v_mag/r_mag)*np.cos(self.gamma)
            dtheta=v_mag*np.cos(self.gamma)/r_mag
            
            self.a = g_acc + drag
            self.ax[k],self.ay[k] = self.a
            
            self.v = self.v + self.a * self.dt
            self.vx[k],self.vy[k] = self.v
            
            self.gamma[k] = self.gamma[k] + dgamma*self.dt
            self.theta[k] = self.theta[k] + dtheta*self.dt
            
            k +=1
            
            if r_mag - cb['radius'] <= sim['stop_alt']:
                print('done in %d iterations' % k)
                break
            
            return(
                np.resize(self.x,k),
                np.resize(self.y,k),
                np.resize(self.vx,k),
                np.resize(self.vy,k),
                np.resize(self.ax,k),
                np.resize(self.ay,k),
                np.resize(self.t,k)
            )
            
            
            