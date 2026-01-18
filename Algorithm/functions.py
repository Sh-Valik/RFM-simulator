import numpy as np

def gravity(z):
    r = np.sqrt(z**2)

    if r < Rplanet:
        accelx = 0.0
        accelz = 0.0
    else:
        accelx = 0.0
        accelz = g0 * ((Rplanet**2) / (r**2))
    
    return np.array([accelx, accelz]) 

def propulsion():
    if t < t_burn:
        T = Psi * mass * g0
        mdot = - T / Ve
        if gravity_turn == 'YES':
            if gt_velx < 1e-6:
                theta_rad_gt = np.deg2rad(90) - np.deg2rad(0.629) # gamma0 - kick angle
            else:
                theta_rad_gt = np.arctan(gt_velz / gt_velx)
            thrust_x = T * np.cos(theta_rad_gt)
            thrust_z = T * np.sin(theta_rad_gt)
        else:
            thrust_x = T * np.cos(theta_rad_const)
            thrust_z = T * np.sin(theta_rad_const)
    else:
        thrust_x = 0.0
        thrust_z = 0.0
        mdot = 0.0
    return np.array([thrust_x, thrust_z]), mdot

def Derivatives(state, t, theta_rad, t_burn, gravity_turn='YES'):

    
    # State vector
    x = state[0]
    z = state[1]
    velx = state[2]
    velz = state[3]
    mass = state[4]

    # compute xdot and zdot
    xdot = velx
    zdot = velz

    # Compute the forces
    
    # Gravity force
    gravityF = - gravity(z) * mass
    

    # Aerodynamic force
    aeroF = np.asarray([0.0, 0.0])

    # Thrust
    thrustF, mdot = propulsion(t, theta_rad, t_burn, mass, z, velx, velz, gravity_turn=gravity_turn)


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