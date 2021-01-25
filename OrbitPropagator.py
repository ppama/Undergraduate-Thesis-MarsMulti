import numpy as np
import matplotlib.pyplot as plt
#plt.style.use('dark_background')
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D

import planetary_data as pd
import tools as t_

def null_perts(): #dictionary of perturbations
    return {
            'oblateness':False, #no need
            'aero':False,
            'moon_grav':False, #no need
            'solar_grav':False, #no need
            'thrust':0,
            'thrust_direction':0
    }
class OrbitPropagator:
    def __init__(self,state0,tspan,dt,coes=False,deg=True,mass0=0,cb=pd.earth,perts=null_perts(),propagator='lsoda'):
        if coes:
            self.r0,self.v0=t_.coes2rv(state0,deg=deg,mu=cb['mu'])
        else:
            self.r0=state0[:3]
            self.v0=state0[3:]
        
        self.y0=self.r0.tolist()+self.v0.tolist()
        self.tspan=tspan
        self.dt=dt
        self.cb=cb
        self.mass0=mass0
        
        #total number of steps
        self.n_steps=int(np.ceil(self.tspan/self.dt))+1
        
        #initialize variables
        self.y=np.zeros((self.n_steps,7))#7 states, xyz, for position and v
        self.ts=np.zeros((self.n_steps,1))
        self.propagator=propagator
        self.step=1
        
        #initial conditions
        self.y[0,:]=self.r0.tolist()+self.v0.tolist()+[self.mass0]
        
        #initiate solver
        self.solver=ode(self.diffy_q)
        self.solver.set_integrator(self.propagator)
        self.solver.set_initial_value(self.y[0,:],0)
        
        # define perturbations dictionary
        self.perts=perts

        self.propagate_orbit()
        
    def propagate_orbit(self):
        print('Propagating orbit...')
        
        #propagate orbit, check for max time and stop conditions at each step
        while self.solver.successful() and self.step<self.n_steps:
            #integrate step
            self.solver.integrate(self.solver.t+self.dt)
            
            #extract values from solver instance
            self.ts[self.step]=self.solver.t
            self.y[self.step]=self.solver.y
            self.step+=1
        
        # extract arrays at the step where the propagation stopped
        self.ts=self.ts[:self.step]
        self.rs=self.y[:self.step,:3]
        self.vs=self.y[:self.step,3:6]
        self.masses=self.y[:self.step,-1]
        self.alts=(np.linalg.norm(self.rs,axis=1)-self.cb['radius']).reshape((self.step,1))
        
    def diffy_q(self,t,y):
        # unpack state
        rx,ry,rz,vx,vy,vz,mass0=y
        r=np.array([rx,ry,rz])
        v=np.array([vx,vy,vz])
        # norm of the radius vector
        norm_r=np.linalg.norm(r)
        
        # two body acceleration
        a=-r*self.cb['mu']/norm_r**3
        
        # J2 perturbation oblateness
        if self.perts['oblateness']:
            z2=r[2]**2
            r2=norm_r**2
            tx=r[0]/norm_r*(5*z2/r2-1)
            ty=r[1]/norm_r*(5*z2/r2-1)
            tz=r[2]/norm_r*(5*z2/r2-3)
            a+=1.5*self.cb['J2']*self.cb['mu']*self.cb['radius']**2/norm_r**4*np.array([tx,ty,tz])
        
        # aerodynamic drag
        if self.perts['aero']:
            # calculate altitude and air density
            z=norm_r-self.cb['radius']
            rho=t_.calc_atmospheric_density(z)
            
            # calculate motion of spacecraft w/ respect to a rotating atmosphere
            v_rel=v-np.cross(self.cb['atm_rot_vector'],r)
            
            drag=-v_rel*0.5*rho*t_.norm(v_rel)*self.perts['Cd']*self.perts['A']/self.mass0
            
            a+=drag
        
        if self.perts['thrust']:
            # thrust vector
            a+=self.perts['thrust_direction']*t_.normed(v)*self.perts['thrust']/mass0/1000.0
            
            # derivative of total mass
            dmdt=-self.perts['thrust']/self.perts['isp']/9.81 #decreasing mass as fuel is spent
            
        return [vx,vy,vz,a[0],a[1],a[2],dmdt]
    
    def calculate_coes(self,degrees=True):
        print('Calculating COEs...')
        
        self.coes=np.zeros((self.n_steps,6))
        
        for n in range(self.n_steps):
            self.coes[n,:]=t_.rv2coes(self.rs[n,:],self.vs[n,:],mu=self.cb['mu'],degrees=degrees)
        
    def plot_coes(self,hours=False,days=False,show_plot=False,save_plot=False,title='COEs',figsize=(16,8),dpi=500):
        print('Plotting COES...')
            
        #create figure and axes instances
        fig,axs=plt.subplots(nrows=2,ncols=3,figsize=figsize)
        
        # figure title
        fig.suptitle(title,fontsize=20)
        
        # x axis
        if hours:
            ts=self.ts/3600.0
            xlabel='Time elapsed (hours)'
        elif days:
            ts=self.ts/(3600.0*24.0)
            xlabel='Time elapsed (days)'
        else:
            ts=self.ts
            xlabel='Time elapsed (seconds)'
            
        # plot true anomaly
        axs[0,0].plot(self.ts,self.coes[:,3])
        axs[0,0].set_title('True Anomaly vs. Time')
        axs[0,0].grid(True)
        axs[0,0].set_ylabel('Angle (degrees)')
        axs[0,0].set_xlabel(xlabel)
        
        # plot semi major axis
        axs[1,0].plot(self.ts,self.coes[:,0])
        axs[1,0].set_title('Semi-Major Axis vs. Time')
        axs[1,0].grid(True)
        axs[1,0].set_ylabel('Semi-Major Axis (km)')
        axs[1,0].set_xlabel(xlabel)
        
        # plot eccentricity
        axs[0,1].plot(self.ts,self.coes[:,1])
        axs[0,1].set_title('Eccentricity vs. Time')
        axs[0,1].grid(True)
        
        # plot argument of periapse
        axs[0,2].plot(self.ts,self.coes[:,4])
        axs[0,2].set_title('Argument of Periapse vs. Time')
        axs[0,2].grid(True)
        
        # plot inclination
        axs[1,1].plot(self.ts,self.coes[:,2])
        axs[1,1].set_title('Inclination vs. Time')
        axs[1,1].grid(True)
        axs[1,1].set_ylabel('Angle (degrees)')
        axs[1,1].set_xlabel(xlabel)
        
        # plot RAAN
        axs[1,2].plot(self.ts,self.coes[:,5])
        axs[1,2].set_title('RAAN vs. Time')
        axs[1,2].grid(True)
        axs[1,2].set_xlabel(xlabel)
        
        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title+'.png',dpi=300)
            
    def calculate_apoapse_periapse(self):
        #define empty arrays
        print('calculating Apoapse and Periapse...')
        self.apoapses=self.coes[:,0]*(1+self.coes[:,1]) #TA = 0, semimajor axis * 1+eccentricity
        self.periapses=self.coes[:,0]*(1-self.coes[:,1]) #TA = 180
        
    def plot_apoapse_periapse(self,hours=False,days=False,show_plot=False,save_plot=False,title='Apoapse and Periapse',dpi=500):
        #create figure
        plt.figure(figsize=(20,10))
        
        if hours:
            ts=self.ts/3600.0
            x_unit='Hours'
        elif days:
            ts=self.ts/(3600.0*24.0)
            x_unit='Days'
        else:
            ts=self.ts
            x_unit='Seconds'
            
        #plot each
        plt.plot(ts,self.apoapses,'b',label='Apoapse')
        plt.plot(ts,self.periapses,'r',label='Periapse')
        
        #labels
        plt.xlabel('Time (%s)' % x_unit)
        plt.ylabel('Altitude (km)')
        
        #other param
        plt.grid(True)
        plt.title(title)
        plt.legend()
        
        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title+'.png',dpi=300)
        
    # plot altitude over time
    def plot_alts(self,show_plot=False,save_plot=False,hours=False,days=False,title='Radial Distance vs. Time',figsize=(16,8),dpi=500):
        if hours:
            ts=self.ts/3600.0
            x_unit='Hours'
        elif days:
            ts=self.ts/(3600.0*24.0)
            x_unit='Days'
        else:
            ts=self.ts
            x_unit='Seconds'
            
        plt.figure(figsize=figsize)
        plt.plot(ts,self.alts,'b')
        plt.grid(True)
        plt.xlabel('Time (%s)' % x_unit)
        plt.ylabel('Altitude (km)')
        plt.title(title)
        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title+'.png',dpi=300)
            
    def plot_3d(self,show_plot=False,save_plot=False,title='3D orbit'):
        # 3D plot
        fig = plt.figure(figsize=(16,8))
        ax = fig.add_subplot(111,projection='3d')
        
        # plot trajectory and starting point 
        ax.plot(self.rs[:,0],self.rs[:,1],self.rs[:,2],'k',label='Trajectory',zorder=10)
        ax.plot([self.rs[0,0]],[self.rs[0,1]],[self.rs[0,2]],'wo',label='Initial Position',zorder=10)
        
        # plot central body
        _u,_v=np.mgrid[0:2*np.pi:20j,0:np.pi:10j]
        _x=self.cb['radius']*np.cos(_u)*np.sin(_v)
        _y=self.cb['radius']*np.sin(_u)*np.sin(_v)
        _z=self.cb['radius']*np.cos(_v)
        ax.plot_surface(_x,_y,_z,cmap='Blues',zorder=1)
        
        # plot the x,y,z vectors
        l = self.cb['radius']*2.0
        x,y,z=[[0,0,0],[0,0,0],[0,0,0]]
        u,v,w=[[l,0,0],[0,l,0],[0,0,l]]
        ax.quiver(x,y,z,u,v,w,color='k')
        
        # check for custom axe limits
        max_val=np.max(np.abs(self.rs))
        
        # set labels and title
        ax.set_xlim([-max_val,max_val])
        ax.set_ylim([-max_val,max_val])
        ax.set_zlim([-max_val,max_val])
        ax.set_xlabel('X (km)'); ax.set_ylabel('Y (km)'); ax.set_zlabel('Z (km)')
        #ax.set_aspect('equal')
        ax.set_title(title)
        plt.legend(['Trajectory', 'Starting Position'])
        
        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title+'.png',dpi=300)