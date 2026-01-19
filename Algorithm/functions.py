import numpy as np
from scipy.interpolate import interp1d
import os


current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
############################################################################
temperature_file_path = os.path.join(project_root, 'resources', 'temperature_profile_0_to_600km.txt')
temp_file = np.loadtxt(temperature_file_path)
alt_ref = temp_file[:, 0]
temp_ref = temp_file[:, 1]
temp_interp = interp1d(alt_ref, temp_ref, kind = 'linear', fill_value='extrapolate')
############################################################################


############################################################################
cd_file_path = os.path.join(project_root, 'resources', 'CD-Mach_relation.txt')
drag_file = np.loadtxt(cd_file_path)
Mach_ref = drag_file[:, 0]
Cd_ref = drag_file[:, 1]
drag_interp = interp1d(Mach_ref, Cd_ref, kind='linear')
############################################################################


############################################################################
def temperature_by_altitude(Rplanet, x, y, z):
    altitude = np.sqrt(x**2 + y**2 + z**2) - Rplanet
    Local_Temperature = temp_interp(altitude)
    
    return Local_Temperature
############################################################################


############################################################################
def gravity(Rplanet, Mplanet, G, x, y, z):
    
    r = np.sqrt(x**2 + y**2 + z**2)

    if r < Rplanet:
        accelx = 0.0
        accely = 0.0
        accelz = 0.0
    else:
        accelx = G * Mplanet / (r**3) * x
        accely = G * Mplanet / (r**3) * y
        accelz = G * Mplanet / (r**3) * z
    
    return np.array([accelx, accely, accelz]) 
############################################################################


############################################################################
# def propulsion():

############################################################################


############################################################################
def density(x, y, z, Rplanet):
    """Air density on current altitude"""
    rho0 = 1.225 # density of the air at sea level [kg/m^3]
    Hm = 8500 # [m]
    alt = np.sqrt(x**2 + y**2 + z**2) - Rplanet
    rho = rho0 * np.exp(- alt / Hm)
    return rho
############################################################################
    

############################################################################
def Derivatives(Rplanet, Mplanet, G, state):
    # State vector
    x = state[0]
    y = state[1]
    z = state[2]
    velx = state[3]
    vely = state[4]
    velz = state[5]
    mass = state[6]

    # compute xdot, ydot and zdot
    xdot = velx
    ydot = vely
    zdot = velz
    
    # Compute the forces
    
    # Gravity force
    gravityF = - gravity(Rplanet, Mplanet, G, x, y, z) * mass
    

    # Aerodynamic force
    rho_alt = density(x, y, z, Rplanet)
    aeroF = np.asarray([0.0, 0.0])

    # Thrust
    thrustF, mdot = propulsion()


    Forces = gravityF + aeroF + thrustF

    #--------------------------------------------    
    # Compute the resulting acceleration
    if mass > 0:
        vdot = Forces / mass
    else:
        vdot = 0.0
        mdot = 0.0
    

    statedot = np.array([xdot, zdot, vdot[0], vdot[1], mdot])
    return statedot
############################################################################