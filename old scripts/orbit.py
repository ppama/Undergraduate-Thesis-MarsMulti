import math
import planetary_data as pd
import tools2 as t_

cb=pd.mars

total=10000
dt=0.1
v_init=5.6 # km/s
h_init=125 # km

j = 0
step=0

def V(g,v,r,a): # dv/dt
    V=-gravity(r)*math.sin(g)-Drag(v,a,r)/m
    return V

def gamma(g,v,r): #dgamma/dt
    g1 = (gravity(r)/v * math.cos(g)-v/r)*math.cos(g)
    return g1

def cr(v,g): # dr/dt
    cr=v*math.sin(g)
    return cr

def theta(v,g,r): #dtheta/dt
    theta = v*math.cos(g)/r
    return theta

def gravity(r): # gravitational force from orbital radius r and acceleration g [km]
    g = cb['mu']/r**2
    return g

def Drag(v,a,r): # drag force from atmosphere
    # calculate altitude and air density
    h=norm_r-self.cb['radius']
    rho=t_.calc_atmospheric_density(h)
        
    # calculate motion of spacecraft w/ respect to a rotating atmosphere
    v_rel=v-np.cross(self.cb['atm_rot_vector'],r)
    # calculate drag
    drag=-v_rel*0.5*rho*(t_.norm(v_rel))**2*self.perts['Cd']*self.perts['A']/self.mass0

while step <= total: # h -> change in speed, k -> change in gamma
    t = step * dt
    h = [0,0,0,0]
    k = [0,0,0,0]
    n = [0,0,0,0]
    s = [0,0,0,0]

    # initial conditions
    h[0] = dt * V(g,v,r,c)
    k[0] = dt * gamma(g,v,r,c)
    n[0] = dt * cr(v,g)
    s[0] = dt * theta(v,g,r)
      
    h[1] = dt * V(g + k[0]*0.5, v + h[0]*0.5 ,r + n[0]*0.5,c)
    k[1] = dt * gamma(g + k[0]*0.5,v+h[0]*0.5 , r +n[0]*0.5,c)
    n[1] = dt * cr(v + h[0]*0.5 ,g +k [0]*0.5 )
    s[1] = dt * theta(v + h[0]*0.5,g + k[0]*0.5,r + n[0]*0.5)

    h[2] = dt * V( g + k[1]*0.5,v+h[1]*0.5 , r + n[1]*0.5,c)
    k[2] = dt * gamma( g + k[1]*0.5,v+h[1]*0.5,r + n[1]*0.5,c)
    n[2] = dt * cr(v + h[1]*0.5 ,g +k [1]*0.5 )
    s[2] = dt * theta(v + h[1]*0.5,g + k[1]*0.5,r + n[1]*0.5)
      
    h[3] = dt * V( g + k[2],v+h[2], r + n[2],c)
    k[3] = dt * gamma( g + k[2],v+h[2], r + n[2],c)
    n[3] = dt * cr(v + h[2] ,g +k [2] )
    s[3] = dt * theta(v + h[2],g + k[2],r + n[2])
      
    v = v + (h[0] + 2*h[1] + 2*h[2] + h[3]) / 6
    g = g + (k[0] + 2*k[1] + 2*k[2] + k[3]) / 6
    r = r + (n[0] + 2*n[1] + 2*n[2] + n[3]) / 6
    b = b + (s[0] + 2*s[1] + 2*s[2] + s[3]) / 6
