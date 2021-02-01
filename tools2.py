import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math as m
import datetime
import planetary_data as pd

d2r=np.pi/180.0
r2d=180.0/np.pi
G_meters=6.67408e-11
G=G_meters*10**-9

cb = pd.mars

def norm(v):
    return np.linalg.norm(v)

def normed(v):
    return np.array(v)/norm(v)

def sqrt(a):
    return np.sqrt(a)

def calc_atmospheric_density(z):
    rhos,zs=find_rho_z(z)
    if rhos[0]==0: return 0.0
    
    Hi=-(zs[1]-zs[0])/m.log(rhos[1]/rhos[0]) # interpolate between 2 data points using exponential fxn
    
    return rhos[0]*m.exp(-(z-zs[0])/Hi)

# find endpoints of altitude and density surrounding input altitude
def find_rho_z(z,zs=cb['zs'],rhos=cb['rhos']):
    if not 1.0<z<1000.0: # no atmosphere beyond 1km altitude
        return [[0.0,0.0],[0.0,0.0]]
    
    # find the two points surrounding the given input altitude
    for n in range(len(rhos)-1):
        if zs[n]<z<zs[n+1]: #within which 2 altitude values
            return [[rhos[n],rhos[n+1]],[zs[n],zs[n+1]]]
        
    #if out of range return zeros
    return [[0.0,0.0],[0.0,0.0]]

def planetEntryConditions(h_apo,h_peri):
    center2interface=(cb['radius']+cb['atmos_interface'])*1000 # radius at martian atmosphere interface [m]
    
    #apoapsis and periapsis calculation
    r_apo=h_apo+cb['radius']*1000
    r_peri=h_peri+cb['radius']*1000
    
    #semimajor axis
    a=(r_apo+r_peri)/2.0
    
    # entry conditions calculations
    v_entry=np.sqrt(G_meters*cb['mass']*(2.0/center2interface-1.0/a)) # velocity at planet atmospheric interface
    print(2.0/center2interface)
    v_apo=np.sqrt(G_meters*(2.0/r_apo-1.0/a)) # velocity at orbit apoapsis
    H=v_apo*r_apo # angular momentum of body in orbit
    gamma_entry=-1*np.arccos(H/(v_entry*center2interface)) # angle at planet atmospheric interface (negative since negative solution gives entry)
    
    return(v_entry,gamma_entry)

def atm_plot(show_plot=False,save_plot=False,title='Martian Atmospheric Model',figsize=(16,8),dpi=500):
    
    fig,axs=plt.subplots(nrows=1,ncols=3,figsize=figsize)
    fig.suptitle(title,fontsize=20)
    
    zs=np.arange(-10,410,10)
    rho=np.zeros(len(zs))
    for z in range(len(zs)-1):
        rho[z]=calc_atmospheric_density(z)
        
    # plot density
    axs[0].plot(cb['rhos'],cb['zs'])
    axs[0].set_title('Altitude vs Density')
    axs[0].grid(True)
    axs[0].set_ylabel('Altitude (km)')
    axs[0].set_xlabel('Density (kg/m^3)')
    axs[0].set_xscale('log')
    
    # plot temperature
    axs[1].plot(cb['temp'],cb['zs'])
    axs[1].set_title('Temperature vs Density')
    axs[1].grid(True)
    axs[1].set_ylabel('Altitude (km)')
    axs[1].set_xlabel('Tempearature [K]')
    
    # plot pressure
    axs[2].plot(cb['p'],cb['zs'])
    axs[2].set_title('Altitude vs Pressure')
    axs[2].grid(True)
    axs[2].set_ylabel('Altitude (km)')
    axs[2].set_xlabel('Pressure[Pa]')
    axs[2].set_xscale('log')
    
    plt.show()
    