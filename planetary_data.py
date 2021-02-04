import numpy as np

G_meters=6.67408e-11

G=G_meters*10**-9


sun={
    'name':'Sun',
    'mass':1.989e30,
    'mu':1.32712e11,
    'radius':695700.0
}

#altitude data points 63, 251, 1000 km, -> density
atm=np.array([[63.096,2.059e-4],[251.189,5.909e-11],[1000.0,3.561e-15]])
earth={
    'name':'Earth',
    'mass':5.972e24,
    'mu':5.972e24*G,
    'radius':6378.0,
    'J2':1.082635854e-3,
    'zs':atm[:,0], # km
    'rhos':atm[:,1]*10**8, # kg / km^3
    'atm_rot_vector':np.array([0.0,0.0,72.9211e-6]), # 1/s
    'deorbit_altitude':10.0 # km
}

data = np.loadtxt("./mars.dat") # [altitude,temp,pressure,density]
mars={
    'name':'Mars',
    'mass':6.4171e23,
    'mu':6.4171e23*G,
    'radius':3389.5, # km
    'zs':data[:,0], # km
    'temp':data[:,1], # atmospheric temps [K]
    'p':data[:,2], #atmospheric pressure [mb]=100Pa
    'rhos':data[:,3], # density [g/cm^3]
    'k':1.3, # ratio of specific heats Cp/CV
    'Rgas':192, # gas constant [J/kg/K]
    'atmos_interface': 3522.2, # atmospheric interface radius [km]
    'period':24.6229, #sidereal rotation period [hr]
    'atm_rot':70.88235e-6, # 1/s planetary rotation freq
    'atm_rot_vector':np.array([0.0,70.88235e-6]) # 1/s planetary rotation time
}

