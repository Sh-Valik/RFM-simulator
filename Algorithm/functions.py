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
    """ Temperature computation on current altitude """
    altitude = np.sqrt(x**2 + y**2 + z**2) - Rplanet
    Local_Temperature = temp_interp(altitude)
    
    return Local_Temperature
############################################################################


############################################################################
def gravity(Rplanet, g0, x, y, z):
    """ Gravitational acceleration computation """
    
    r = np.sqrt(x**2 + y**2 + z**2)

    if r < Rplanet:
        accelx = 0.0
        accely = 0.0
        accelz = 0.0
    else:
        accelx = g0 * ((Rplanet**2) / (r**3)) * x
        accely = g0 * ((Rplanet**2) / (r**3)) * y
        accelz = g0 * ((Rplanet**2) / (r**3)) * z
    
    return np.array([accelx, accely, accelz])
############################################################################


############################################################################
def propulsion(t, t_burn, T_mag, m_flow, theta):
    """ Thrust and mass flow rate computation """
    if t <= t_burn:
        thrust_x = T_mag * np.cos(theta)
        thrust_y = 0.0 # placeholder
        thrust_z = T_mag * np.sin(theta)
        mdot = - m_flow  # mass flow rate [kg/s]
    else:
        thrust_x = 0.0
        thrust_y = 0.0
        thrust_z = 0.0
        mdot = 0.0

    return np.array([thrust_x, thrust_y, thrust_z]), mdot


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
def Derivatives(state, t, t_burn, T_mag, m_flow, theta, Rplanet, g0, Area):
    """ Computes the state derivatives """
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
    gravityF = - gravity(Rplanet, g0, x, y, z) * mass
    

    # Aerodynamic force
    V = np.sqrt(velx**2 + vely**2 + velz**2)
    rho_alt = density(x, y, z, Rplanet)
    Temp_local = temperature_by_altitude(Rplanet, x, y, z)
    loc_sound_speed = np.sqrt(1.4 * 287.05 * Temp_local)  # speed of sound [m/s]
    Mach = V / loc_sound_speed
    Cd = drag_interp(Mach)
    aeroF = -0.5 * min(rho_alt, 1.293) * Area * abs(V) * Cd * np.array([velx, vely, velz]) # Aerodynamic drag force

    # Thrust
    thrustF, mdot = propulsion(t, t_burn, T_mag, m_flow, theta)


    Forces = gravityF + aeroF + thrustF

    #--------------------------------------------    
    # Compute the resulting acceleration
    if mass > 0:
        vdot = Forces / mass
    else:
        vdot = 0.0
        mdot = 0.0
    

    statedot = np.array([xdot, ydot, zdot, vdot[0], vdot[1], vdot[2], mdot])
    return statedot
############################################################################


############################################################################
def mass_of_stages_eps_mode(data_list):
    #placeholder for future implementation
    pass
############################################################################
