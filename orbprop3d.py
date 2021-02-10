import numpy as np
import matplotlib.pyplot as plt
#plt.style.use('dark_background')
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
import math

import tools as t_
from matplotlib.offsetbox import AnchoredText

class OrbitPropagator(object):
    def __init__(self,planet,craft,sim,state0,beta,fpa):
        self.r0=state0[:3]
        self.v0=state0[3:]
        self.gamma0=[fpa]
        #r_mag = planet['radius']+entry_altitude
        #v_init = v_init
        #self.r0 = np.array([0,r_mag,0])
        #self.v0 = np.array([v_init*np.cos(fpa),v_init*np.sin(fpa),0])
        
        self.y0=self.r0.tolist()+self.v0.tolist()+self.gamma0
        self.tspan=sim['timespan']
        
        self.cb=planet
        self.beta = beta
        self.craft = craft
        self.sim = sim
        self.dt=sim['delta_t']
        
        #total number of steps
        self.n_steps=int(np.ceil(self.tspan/self.dt))
        
        #initialize variables
        self.ys=np.zeros((self.n_steps,7))#7 states, xyz, for position and v then gamma
        self.ts=np.zeros((self.n_steps,1))
        self.ts[0]=0
        self.ys[0,:]=self.y0
        self.step=1
        self.alts=np.zeros((self.n_steps,1))
        self.alts[0]=t_.norm(self.r0)-self.cb['radius']
        self.vmag=np.zeros((self.n_steps,1))
        self.vmag[0]=t_.norm(self.v0)
        self.range=np.zeros((self.n_steps,1))
        self.theta=np.zeros((self.n_steps,1))
        self.gamma=np.zeros((self.n_steps,1))
        self.gamma[0]=np.rad2deg(self.gamma0)
        self.escape = False
        self.nrev = 0 # revolutions around planet during entry
        
        #initiate solver
        self.solver=ode(self.diffy_q)
        self.solver.set_integrator('lsoda')
        self.solver.set_initial_value(self.y0,0)

        self.propagate_orbit()
        
    def propagate_orbit(self):
        """
        Calculates the trajectory of an object entering a planetary atmosphere
        Iteratively solves the ODEs using scipy integrate, with a break condition once the stop altitude is reached
        """
        print('Propagating orbit...')
        #propagate orbit, check for max time and stop conditions at each step
        while self.solver.successful() and self.step<self.n_steps:
            #integrate step
            self.solver.integrate(self.solver.t+self.dt)
            
            #extract values from solver instance
            self.ts[self.step]=self.solver.t
            self.ys[self.step]=self.solver.y
            
            # calculate altitude,gamma,theta at this time step
            self.vmag[self.step]=t_.norm(self.solver.y[3:6])
            self.alts[self.step]=t_.norm(self.solver.y[:3])-self.cb['radius']
            self.gamma[self.step]=np.rad2deg(self.solver.y[6])
            #calculate orbital angle subtended
            if self.solver.y[0]<0:
                self.theta[self.step]=np.arctan2(self.solver.y[0],self.solver.y[1])+2*np.pi
            else:
                self.theta[self.step]=np.arctan2(self.solver.y[0],self.solver.y[1])
            #calculate range based on angle subtended
            #if self.solver.y[0] >= 0.0:
            #    self.range[self.step]=self.cb['radius']*(self.nrev*2*np.pi+np.pi/2-self.theta[self.step])
            #elif self.solver.y[0] < 0.0:
            #    self.range[self.step]=self.cb['radius']*(self.nrev*2*np.pi+(5*np.pi/2)-self.theta[self.step])
            self.range[self.step]=(2*np.pi*self.nrev+self.theta[self.step])*self.cb['radius'] # convert x tangent to surface distance
            if np.abs(self.theta[self.step]-self.theta[self.step-1]) > 6.28: 
                self.range[self.step] = 2*np.pi*self.cb['radius'] #range at 1 revolution = 2pi*r
                self.nrev +=1 # revolution counter

            altcheck=self.alts[self.step]
            self.step+=1
            
            # altitude break condition
            if altcheck <= self.craft['stop_alt']:
                break
            if self.alts[self.n_steps-1] > self.craft['stop_alt']:
                print("Craft was unable to complete EDL.")
                self.escape=True
                break
                
        
        # extract arrays at the step where the propagation stopped
        self.ts=self.ts[:self.step]
        self.rs=self.ys[:self.step,:3]
        self.vs=self.ys[:self.step,3:6]
        self.alts=self.alts[:self.step]
        self.vmag=self.vmag[:self.step]
        self.range=self.range[:self.step]
        self.gamma=self.gamma[:self.step]
        self.theta=self.theta[:self.step]
        
    def diffy_q(self,t,y):
        """
        Function containing the system of ODEs representing behavior of an orbital equation
        The orbitting object is considered as a point function
        Differential equations:
        gravitational acceleration g = -r*mu/r^3 (negative as acceleration is towards the center of the celestial body)
        relative velocity = v - cross product between r vector and planetary rotation vector 
        (the latter is clockwise to assume reentry is in the direction of rotation, +x-axis)
        rho = atmospheric density at altitude z
        drag = 0.5*v_rel*rho*v_rel*beta
        dv/dt = g - drag
        returns V and A vectors
        """
        # unpack state
        rx,ry,rz,vx,vy,vz,gamma=y
        r=np.array([rx,ry,rz])
        v=np.array([vx,vy,vz])
        
        # norm of the radius vector
        norm_r=np.linalg.norm(r)
        
        # two body acceleration
        a=-r*self.cb['mu']/norm_r**3
        #a[0] = -r[0]*self.cb['mu']/norm_r**3
        #a[1] = -r[1]*self.cb['mu']/norm_r**3*np.sin(gamma)
        #a[2] = -r[2]*self.cb['mu']/norm_r**3
        # calculate altitude and air density
        z=norm_r-self.cb['radius']
        rho=t_.calc_atmospheric_density(z)
        
        # calculate motion of spacecraft w/ respect to a rotating atmosphere
        v_rel=v-np.cross(self.cb['atm_rot_vector'],r)
        drag=-v_rel*0.5*rho*t_.norm(v_rel)/self.beta
        gamma = (self.cb['mu']/(t_.norm(v_rel)*norm_r**2)-t_.norm(v_rel)/norm_r)*np.cos(gamma)
        #gamma = np.arctan(vy/vx)
        a+=drag
        return [vx,vy,vz,a[0],a[1],a[2],gamma]
    
   
    # plot state over time
    def plot_state(self,hours=False,days=False,show_plot=False,save_plot=False,figsize=(16,8),dpi=500):
        print('Plotting Trajectory Profile...')
        #create figure and axes instances
        fig,axs=plt.subplots(nrows=1,ncols=4,figsize=figsize)
        
        # figure title
        title = f'Trajectory Profile (beta={self.craft["ballistic_coef"]})'
        fig.suptitle(title,fontsize=20)
        
        # plot alts
        axs[0].plot(self.ts,self.alts)
        axs[0].set_title('Altitude vs Time')
        axs[0].grid(True)
        axs[0].set_ylabel('Altitude (km)')
        axs[0].set_xlabel('Time (s)')
        anchored_text = AnchoredText(f'Time of Flight: {self.truncate(self.ts[self.step-1,0],2)} s', loc=1)
        axs[0].add_artist(anchored_text)
    
        # plot velocity
        axs[1].plot(self.vmag,self.alts)
        axs[1].set_title('Altitude vs. Velocity')
        axs[1].grid(True)
        axs[1].set_ylabel('Altitude (km)')
        axs[1].set_xlabel('Velocity (km/s)')
        axs[1].set_xlim([0,12])
        anchored_text = AnchoredText(f'Velocity at impact: {self.truncate(self.vmag[self.step-1,0],2)} km/s', loc=1)
        axs[1].add_artist(anchored_text)
        
        # plot range
        axs[2].plot(self.range,self.alts)
        #axs[2].plot(self.rs[:self.step,0],self.alts,'b') # print pure x value or surface arc?
        axs[2].set_title('Altitude vs Range')
        axs[2].grid(True)
        axs[2].set_ylabel('Altitude (km)')
        axs[2].set_xlabel('Range (km)')
        anchored_text = AnchoredText(f'Maximum Range: {self.truncate(self.range[self.step-1,0],2)} km', loc=1)
        axs[2].add_artist(anchored_text)
        
        # plot FPA
        axs[3].plot(self.ts,self.gamma)
        axs[3].set_title('Flight Path Angle vs Time')
        axs[3].grid(True)
        axs[3].set_ylabel('Flight Path Angle (degrees)')
        axs[3].set_xlabel('Time (s)')
        anchored_text = AnchoredText(f'Final FPA: {self.truncate(self.gamma[self.step-1,0],2)} degrees', loc=1)
        axs[3].add_artist(anchored_text)
        
        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title+'.png',dpi=300)
    
    def plot_alts(self,show_plot=False,save_plot=False,hours=False,days=False,title='Altitude vs. Time',figsize=(16,8),dpi=500):
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
            
    def plot_velo(self,show_plot=False,save_plot=False,hours=False,days=False,title='Altitude vs velocity',figsize=(16,8),dpi=500):
        plt.figure(figsize=figsize)
        plt.plot(self.vmag,self.alts,'b')
        plt.grid(True)
        plt.xlabel('Velocity (km)')
        plt.ylabel('Altitude (km)')
        plt.title(title)
        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title+'.png',dpi=300)
        
    def plot_range(self,show_plot=False,save_plot=False,hours=False,title='Altitude vs. Range',figsize=(16,8),dpi=500):
        f,ax=plt.subplots(1,1)
        plt.figure(figsize=figsize)
        ax.plot(self.rs[:self.step,0],self.alts,'b')
        ax.grid(True)
        ax.set_xlabel('Range (km)')
        ax.set_ylabel('Altitude (km)')
        ax.set_title(title)
        anchored_text = AnchoredText("Test", loc=1)
        ax.add_artist(anchored_text)
        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title+'.png',dpi=300)
        
        print('Maximum range reached: %.1f km' % self.rs[self.step-1,0])
    def plot_3d(self,show_plot=False,save_plot=False,title='Ballistic Entry to Martian Surface'):
        # 3D plot
        fig = plt.figure(figsize=(16,8))
        ax = fig.add_subplot(111,projection='3d')
        
        # plot trajectory and starting point 
        ax.plot(self.rs[:,0],self.rs[:,1],self.rs[:,2],'k',label='Trajectory',zorder=10)
        ax.plot([self.rs[0,0]],[self.rs[0,1]],[self.rs[0,2]],'bo',label='Initial Position',zorder=10)
        
        # plot ending point
        ax.plot([self.rs[len(self.rs)-1,0]],[self.rs[len(self.rs)-1,1]],[self.rs[len(self.rs)-1,2]],'go',label='Ending Position',zorder=10)
        
        # plot central body
        _u,_v=np.mgrid[0:2*np.pi:20j,0:np.pi:10j]
        _x=self.cb['radius']*np.cos(_u)*np.sin(_v)
        _y=self.cb['radius']*np.sin(_u)*np.sin(_v)
        _z=self.cb['radius']*np.cos(_v)
        ax.plot_surface(_x,_y,_z,cmap='Reds',zorder=1)
        
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
        plt.legend(['Trajectory', 'Starting Position','Ending Position'])
        
        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title+'.png',dpi=300)
            
    def truncate(self,number, decimals=0):
        """
        Returns a value truncated to a specific number of decimal places.
        """
        if not isinstance(decimals, int):
            raise TypeError("decimal places must be an integer.")
        elif decimals < 0:
            raise ValueError("decimal places has to be 0 or more.")
        elif decimals == 0:
            return math.trunc(number)

        factor = 10.0 ** decimals
        return math.trunc(number * factor) / factor