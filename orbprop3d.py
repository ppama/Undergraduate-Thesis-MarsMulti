import numpy as np
import matplotlib.pyplot as plt
#plt.style.use('dark_background')
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D

import baseline as ba
import planetary_data as pd
import tools as t_
from matplotlib.offsetbox import AnchoredText

class OrbitPropagator:
    def __init__(self,state0,tspan,dt,deg=True,mass0=0,cb=pd.mars,craft=ba.craft,sim=ba.sim):
        self.r0=state0[:3]
        self.v0=state0[3:]
        
        self.y0=self.r0.tolist()+self.v0.tolist()
        self.tspan=tspan
        self.dt=dt
        self.cb=cb
        self.mass=mass0
        self.craft = craft
        self.sim = sim
        
        #total number of steps
        self.n_steps=int(np.ceil(self.tspan/self.dt))
        
        #initialize variables
        self.ys=np.zeros((self.n_steps,6))#6 states, xyz, for position and v
        self.ts=np.zeros((self.n_steps,1))
        self.ts[0]=0
        self.ys[0,:]=self.y0
        self.step=1
        self.alts=np.zeros((self.n_steps,1))
        self.alts[0]=t_.norm(self.r0)-self.cb['radius']
        self.vmag=np.zeros((self.n_steps,1))
        self.vmag[0]=t_.norm(self.v0)
        self.range=np.zeros((self.n_steps,1))
        
        #initiate solver
        self.solver=ode(self.diffy_q)
        self.solver.set_integrator('lsoda')
        self.solver.set_initial_value(self.y0,0)

        self.propagate_orbit()
        
    def propagate_orbit(self):
        print('Propagating orbit...')
        
        #propagate orbit, check for max time and stop conditions at each step
        while self.solver.successful() and self.step<self.n_steps:
            #integrate step
            self.solver.integrate(self.solver.t+self.dt)
            
            #extract values from solver instance
            self.ts[self.step]=self.solver.t
            self.ys[self.step]=self.solver.y
            
            # calcualte altitude,gamma,theta at this time step
            self.vmag[self.step]=t_.norm(self.solver.y[3:])
            self.alts[self.step]=t_.norm(self.solver.y[:3])-self.cb['radius']
            self.range[self.step]=np.arctan(self.solver.y[0]/self.solver.y[1])*0.01745*self.cb['radius'] # convert x tangent to surface distance
            altcheck=self.alts[self.step]
            #print(self.step)
            #print(self.sim['stop_alt'])
            self.step+=1
            
            # altitude break condition
            if altcheck <= self.sim['stop_alt']:
                #print('Spacecraft reached target altitude after %.1f seconds' % self.dt*self.step)
                break
                
        
        # extract arrays at the step where the propagation stopped
        self.ts=self.ts[:self.step]
        self.rs=self.ys[:self.step,:3]
        self.vs=self.ys[:self.step,3:]
        self.alts=self.alts[:self.step]
        self.vmag=self.vmag[:self.step]
        self.range=self.range[:self.step]
        
        
        
    def diffy_q(self,t,y):
        # unpack state
        rx,ry,rz,vx,vy,vz=y
        r=np.array([rx,ry,rz])
        v=np.array([vx,vy,vz])
        
        # norm of the radius vector
        norm_r=np.linalg.norm(r)
        
        # two body acceleration
        a=-r*self.cb['mu']/norm_r**3
        
        # calculate altitude and air density
        z=norm_r-self.cb['radius']
        rho=t_.calc_atmospheric_density(z)
        # calculate motion of spacecraft w/ respect to a rotating atmosphere
        v_rel=v-np.cross(self.cb['atm_rot_vector'],r)
        drag=-v_rel*0.5*rho*t_.norm(v_rel)*self.craft['ballistic_coef']
        a+=drag
        return [vx,vy,vz,a[0],a[1],a[2]]
    
   
    # plot state over time
    def plot_state(self,hours=False,days=False,show_plot=False,save_plot=False,figsize=(16,8),dpi=500):
        print('Plotting Trajectory Profile...')
        #create figure and axes instances
        fig,axs=plt.subplots(nrows=1,ncols=3,figsize=figsize)
        
        # figure title
        title = f'Trajectory Profile (beta={self.craft["ballistic_coef"]})'
        fig.suptitle(title,fontsize=20)
        
        # plot alts
        axs[0].plot(self.ts,self.alts)
        axs[0].set_title('Altitude vs Time')
        axs[0].grid(True)
        axs[0].set_ylabel('Altitude (km)')
        axs[0].set_xlabel('Time (s)')
    
        # plot velocity
        axs[1].plot(self.vmag,self.alts)
        axs[1].set_title('Altitude vs. Velocity')
        axs[1].grid(True)
        axs[1].set_ylabel('Altitude (km)')
        axs[1].set_xlabel('Velocity (km/s)')
        
        # plot range
        #axs[2].plot(self.range,self.alts)
        axs[2].plot(self.rs[:self.step,0],self.alts,'b')
        axs[2].set_title('Altitude vs Range')
        axs[2].grid(True)
        axs[2].set_ylabel('Altitude (km)')
        axs[2].set_xlabel('Range (km)')
        anchored_text = AnchoredText("Test", loc=1)
        axs[2].add_artist(anchored_text)
        
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
    def plot_3d(self,show_plot=False,save_plot=False,title='3D orbit'):
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
        plt.legend(['Trajectory', 'Starting Position','Ending Position'])
        
        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title+'.png',dpi=300)