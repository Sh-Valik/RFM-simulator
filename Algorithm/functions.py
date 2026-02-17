import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import brentq
from scipy.integrate import solve_ivp, odeint
import os
import streamlit as st
import plotly.express as px
import pandas as pd
import plotly.graph_objects as go

G = 6.6742 * 10**-11  # gravitational constant [N.m^2/kg^2]
g0 = 9.80665  # standard gravitational acceleration [m/s^2]
Rplanet = 6371000  # mean radius of the Earth [m]
Mplanet = 5.97219 * 10**24  # mass of the Earth [kg]


current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
############################################################################
temperature_file_path = os.path.join(project_root, 'resources', 'temperature_profile_0_to_600km.txt')
temp_file = np.loadtxt(temperature_file_path)
alt_ref = temp_file[:, 0]
temp_ref = temp_file[:, 1]
temp_interp = interp1d(alt_ref, temp_ref, kind = 'linear', fill_value='extrapolate')

new_alt = np.linspace(min(alt_ref), max(alt_ref), 200)
new_temp = temp_interp(new_alt)
temp_data = pd.DataFrame({'Temperature (K)': new_temp, 'Altitude (km)': new_alt/1000})
############################################################################


############################################################################
cd_file_path = os.path.join(project_root, 'resources', 'CD-Mach_relation.txt')
drag_file = np.loadtxt(cd_file_path)
Mach_ref = drag_file[:, 0]
Cd_ref = drag_file[:, 1]
# drag_interp = interp1d(Mach_ref, Cd_ref, kind='linear')
drag_interp = interp1d(
    Mach_ref,
    Cd_ref,
    kind='linear',
    bounds_error=False,
    fill_value=(Cd_ref[0], Cd_ref[-1])
)

new_Mach = np.linspace(min(Mach_ref), max(Mach_ref), 200)
new_Cd = drag_interp(new_Mach)
cd_data = pd.DataFrame({'Mach': new_Mach, 'Cd': new_Cd})
############################################################################


############################################################################
def temperature_by_altitude(x, y, z):
    """ Temperature computation on current altitude """
    altitude = np.sqrt(x**2 + y**2 + z**2) - Rplanet
    Local_Temperature = temp_interp(altitude)
    
    return Local_Temperature
############################################################################


############################################################################
def gravity(x, y, z):
    """ Gravitational acceleration computation with J2 perturbation """
    J2 = 0.00108263
    r = np.sqrt(x**2 + y**2 + z**2)

    if r < Rplanet:
        accelx = 0.0
        accely = 0.0
        accelz = 0.0
    else:
        accelx = - (((G * Mplanet) / r**3) * x) + (3 * J2 * ((G * Mplanet) / 2) * Rplanet**2 * (x / (r**5)) * (1 - ((5 * (z**2)) / (r**2))))
        accely = - (((G * Mplanet) / r**3) * y) + (3 * J2 * ((G * Mplanet) / 2) * Rplanet**2 * (y / (r**5)) * (1 - ((5 * (z**2)) / (r**2))))
        accelz = - (((G * Mplanet) / r**3) * z) + (3 * J2 * ((G * Mplanet) / 2) * Rplanet**2 * (z / (r**5)) * (3 - ((5 * (z**2)) / (r**2))))
    
    return np.array([accelx, accely, accelz])
############################################################################



############################################################################
def density(x, y, z):
    """Air density on current altitude"""
    rho0 = 1.293 # density of the air at sea level [kg/m^3]
    Hm = 8432.56 # [m]
    alt = np.sqrt(x**2 + y**2 + z**2) - Rplanet
    rho = rho0 * np.exp(- alt / Hm)
    return rho
############################################################################

############################################################################
def q_dynamic_pressure(rho, V):
    """Dynamic pressure computation"""
    q = 0.5 * rho * abs(V)**2
    return q
############################################################################

############################################################################
def Derivatives_with_boosters(state, t, stages_info, boosters_info, Area_pf, Area_bf, Cd_of_crosflow_cylinder, t_vertical, Az_rad, stage_index, kick_angle_deg):
    """ Computes the state derivatives for stages onhly"""
    kick_angle = np.deg2rad(kick_angle_deg)
    # Unpack stages_info
    T_mag_stages = stages_info[1]
    mass_flow_stages = stages_info[2]
    
    
    # Unpack boosters_info
    T_mag_boosters = boosters_info[1]
    mass_flow_boosters = boosters_info[2]
    
    
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

    # Aerodynamic block
    V = np.sqrt(velx**2 + vely**2 + velz**2)
    rho_alt = density(x, y, z)
    Temp_local = temperature_by_altitude(x, y, z)
    loc_sound_speed = np.sqrt(1.4 * 287.05 * Temp_local)  # speed of sound [m/s] = sqrt(R * gamma * local temperature)
    Mach = V / loc_sound_speed
    
    # Area
    r_vec = np.array([x, y, z])
    v_vec = np.array([velx, vely, velz])
    v_radial = np.dot(r_vec, v_vec) / np.linalg.norm(r_vec)
    apogee_reached = v_radial < 0

    if apogee_reached:
        Area = Area_bf
        Cd = Cd_of_crosflow_cylinder
    else:
        Area = Area_pf
        Cd = float(drag_interp(Mach))
    
    # Compute the forces
    
    # Gravity force
    gravityF = gravity(x, y, z) * mass
    

    # Aerodynamic force
    
    if rho_alt < 1e-6:
        aeroF = np.zeros(3)
    else:
        aeroF = -0.5 * min(rho_alt, 1.293) * Area * abs(V) * Cd * np.array([velx, vely, velz]) # Aerodynamic drag force



    # Thrust magnitude and mass flow

    T_mag = T_mag_stages[stage_index] + sum(T_mag_boosters)
    mdot =- (mass_flow_stages[stage_index] + sum(mass_flow_boosters))

    # Compute thrust vector in ECEF using local reference frame
    if T_mag > 0:
        r_pos = np.sqrt(x**2 + y**2 + z**2)
        r_hat = np.array([x, y, z]) / r_pos

        if t <= t_vertical:
            # Vertical flight — thrust along local vertical (radial direction)
            thrustF = T_mag * r_hat
        # elif t > t_vertical and t <= t_vertical + 1e-3:
        #     k_hat = np.array([0.0, 0.0, 1.0])  # Earth's rotation axis
        #     east_raw = np.cross(k_hat, r_hat)
        #     east_norm = np.linalg.norm(east_raw)
        #     if east_norm > 1e-10:
        #         east_hat = east_raw / east_norm
        #     else:
        #         east_hat = np.array([1.0, 0.0, 0.0])
        #     north_hat = np.cross(r_hat, east_hat)

        #     # Kick: small angle from vertical in the azimuth direction
        #     thrust_dir = (np.cos(kick_angle) * r_hat +
        #                   np.sin(kick_angle) * (
        #                       np.sin(Az_rad) * east_hat +
        #                       np.cos(Az_rad) * north_hat))
        #     thrustF = T_mag * thrust_dir
        

        else:
            # Phase 3: Pitch program below atmosphere, gravity turn above
            altitude = r_pos - 6371000.0  # approximate altitude
            k_hat = np.array([0.0, 0.0, 1.0])
            east_raw = np.cross(k_hat, r_hat)
            east_norm = np.linalg.norm(east_raw)
            if east_norm > 1e-10:
                east_hat = east_raw / east_norm
            else:
                east_hat = np.array([1.0, 0.0, 0.0])
            north_hat = np.cross(r_hat, east_hat)

            if altitude < 80000.0:
                # Below atmosphere: fixed pitch at kick angle from local vertical
                thrust_dir = (np.cos(kick_angle) * r_hat +
                              np.sin(kick_angle) * (
                                  np.sin(Az_rad) * east_hat +
                                  np.cos(Az_rad) * north_hat))
            else:
                # Above atmosphere: true gravity turn (follow velocity vector)
                if V > 1e-6:
                    thrust_dir = v_vec / V
                else:
                    thrust_dir = r_hat
            thrustF = T_mag * thrust_dir
    else:
        thrustF = np.zeros(3)


    Forces = gravityF + aeroF + thrustF

    #--------------------------------------------    
    # Compute the resulting acceleration
    if mass > 0:
        vdot = Forces / mass
    else:
        vdot = np.zeros(3)
        mdot = 0.0
    

    statedot = np.array([xdot, ydot, zdot, vdot[0], vdot[1], vdot[2], mdot])
    return statedot

def Derivatives_propelled(state, t, stages_info, boosters_info, Area_pf, Area_bf, Cd_of_crosflow_cylinder, t_vertical, Az_rad, stage_index, kick_angle_deg):
    """ Computes the state derivatives for stages onhly"""
    kick_angle = np.deg2rad(kick_angle_deg)
    # Unpack stages_info
    t_burn_stages = stages_info[0]
    T_mag_stages = stages_info[1]
    mass_flow_stages = stages_info[2]
    m_construction_stages = stages_info[3]
    
    
    # Unpack boosters_info
    t_burn_boosters = boosters_info[0]
    T_mag_boosters = boosters_info[1]
    mass_flow_boosters = boosters_info[2]
    
    
    # State vector
    x = state[0]
    y = state[1]
    z = state[2]
    velx = state[3]
    vely = state[4]
    velz = state[5]
    mass = state[6]
    if mass <= m_construction_stages[stage_index]:
        mass = m_construction_stages[stage_index]
    
    # compute xdot, ydot and zdot
    xdot = velx
    ydot = vely
    zdot = velz

    # Aerodynamic block
    V = np.sqrt(velx**2 + vely**2 + velz**2)
    rho_alt = density(x, y, z)
    Temp_local = temperature_by_altitude(x, y, z)
    loc_sound_speed = np.sqrt(1.4 * 287.05 * Temp_local)  # speed of sound [m/s] = sqrt(R * gamma * local temperature)
    Mach = V / loc_sound_speed
    
    # Area
    r_vec = np.array([x, y, z])
    v_vec = np.array([velx, vely, velz])
    v_radial = np.dot(r_vec, v_vec) / np.linalg.norm(r_vec)
    apogee_reached = v_radial < 0

    if apogee_reached:
        Area = Area_bf
        Cd = Cd_of_crosflow_cylinder
    else:
        Area = Area_pf
        Cd = float(drag_interp(Mach))
    
    # Compute the forces
    
    # Gravity force
    gravityF = gravity(x, y, z) * mass
    

    # Aerodynamic force
    
    if rho_alt < 1e-6:
        aeroF = np.zeros(3)
    else:
        aeroF = -0.5 * min(rho_alt, 1.293) * Area * abs(V) * Cd * np.array([velx, vely, velz]) # Aerodynamic drag force



    # Thrust magnitude and mass flow
    T_mag = T_mag_stages[stage_index]
    mdot = - mass_flow_stages[stage_index]
    if mass <= m_construction_stages[stage_index]:
        T_mag = 0.0
        mdot = 0.0

    # Compute thrust vector in ECEF using local reference frame
    if T_mag > 0:
        r_pos = np.sqrt(x**2 + y**2 + z**2)
        r_hat = np.array([x, y, z]) / r_pos

        if t <= t_vertical:
            # Vertical flight — thrust along local vertical (radial direction)
            thrustF = T_mag * r_hat
        # elif t > t_vertical and t <= t_vertical + 1e-3:
        #     k_hat = np.array([0.0, 0.0, 1.0])  # Earth's rotation axis
        #     east_raw = np.cross(k_hat, r_hat)
        #     east_norm = np.linalg.norm(east_raw)
        #     if east_norm > 1e-10:
        #         east_hat = east_raw / east_norm
        #     else:
        #         east_hat = np.array([1.0, 0.0, 0.0])
        #     north_hat = np.cross(r_hat, east_hat)

        #     # Kick: small angle from vertical in the azimuth direction
        #     thrust_dir = (np.cos(kick_angle) * r_hat +
        #                   np.sin(kick_angle) * (
        #                       np.sin(Az_rad) * east_hat +
        #                       np.cos(Az_rad) * north_hat))
        #     thrustF = T_mag * thrust_dir
        

        else:
            # Phase 3: Pitch program below atmosphere, gravity turn above
            altitude = r_pos - 6371000.0
            k_hat_local = np.array([0.0, 0.0, 1.0])
            east_raw = np.cross(k_hat_local, r_hat)
            east_norm = np.linalg.norm(east_raw)
            if east_norm > 1e-10:
                east_hat = east_raw / east_norm
            else:
                east_hat = np.array([1.0, 0.0, 0.0])
            north_hat = np.cross(r_hat, east_hat)

            if altitude < 95000.0:
                thrust_dir = (np.cos(kick_angle) * r_hat +
                              np.sin(kick_angle) * (
                                  np.sin(Az_rad) * east_hat +
                                  np.cos(Az_rad) * north_hat))
            else:
                if V > 1e-6:
                    thrust_dir = v_vec / V
                else:
                    thrust_dir = r_hat
            thrustF = T_mag * thrust_dir
    else:
        thrustF = np.zeros(3)


    Forces = gravityF + aeroF + thrustF

    # q = q_dynamic_pressure(min(rho_alt, 1.293), V)
    #--------------------------------------------    
    # Compute the resulting acceleration
    if mass > 0:
        vdot = Forces / mass
    else:
        vdot = np.zeros(3)
        mdot = 0.0
    

    statedot = np.array([xdot, ydot, zdot, vdot[0], vdot[1], vdot[2], mdot])
    return statedot

def Derivatives_balistic(state, t, Area_pf, Area_bf, Cd_of_crosflow_cylinder):
    """ Computes the state derivatives for stages onhly"""
     
    
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

    # Aerodynamic block
    V = np.sqrt(velx**2 + vely**2 + velz**2)
    rho_alt = density(x, y, z)
    Temp_local = temperature_by_altitude(x, y, z)
    loc_sound_speed = np.sqrt(1.4 * 287.05 * Temp_local)  # speed of sound [m/s] = sqrt(R * gamma * local temperature)
    Mach = V / loc_sound_speed
    
    # Area
    r_vec = np.array([x, y, z])
    v_vec = np.array([velx, vely, velz])
    v_radial = np.dot(r_vec, v_vec) / np.linalg.norm(r_vec)
    apogee_reached = v_radial < 0

    if apogee_reached:
        Area = Area_bf
        Cd = Cd_of_crosflow_cylinder
    else:
        Area = Area_pf
        Cd = float(drag_interp(Mach))
    
    # Compute the forces
    
    # Gravity force
    gravityF = gravity(x, y, z) * mass
    

    # Aerodynamic force
    
    if rho_alt < 1e-6:
        aeroF = np.zeros(3)
    else:
        aeroF = -0.5 * min(rho_alt, 1.293) * Area * abs(V) * Cd * np.array([velx, vely, velz]) # Aerodynamic drag force



    Forces = gravityF + aeroF
    mdot = 0.0

    
    #--------------------------------------------    
    # Compute the resulting acceleration
    if mass > 0:
        vdot = Forces / mass
    else:
        vdot = np.zeros(3)
        mdot = 0.0
    

    statedot = np.array([xdot, ydot, zdot, vdot[0], vdot[1], vdot[2], mdot])
    return statedot



############################################################################
def Derivatives_boosters(state, t, t_burn_boosters, T_mag_boosters, mass_flow_boosters, m_construction_each_boosters, Area_pf, Area_bf, Cd_of_crosflow_cylinder, t_vertical, Az_rad, kick_angle_deg):
    
    kick_angle = np.deg2rad(kick_angle_deg)
    x = state[0]
    y = state[1]
    z = state[2]
    velx = state[3]
    vely = state[4]
    velz = state[5]
    if t <= t_burn_boosters:
        mass = state[6]
    else:
        mass = m_construction_each_boosters

    # compute xdot, ydot and zdot
    xdot = velx
    ydot = vely
    zdot = velz

    # Aerodynamic block
    V = np.sqrt(velx**2 + vely**2 + velz**2)
    rho_alt = density(x, y, z)
    Temp_local = temperature_by_altitude(x, y, z)
    loc_sound_speed = np.sqrt(1.4 * 287.05 * Temp_local)  # speed of sound [m/s] = sqrt(R * gamma * local temperature)
    Mach = V / loc_sound_speed
    
    # Area
    r_vec = np.array([x, y, z])
    v_vec = np.array([velx, vely, velz])
    v_radial = np.dot(r_vec, v_vec) / np.linalg.norm(r_vec)
    apogee_reached = v_radial < 0

    if apogee_reached:
        Area = Area_bf
        Cd = Cd_of_crosflow_cylinder
    else:
        Area = Area_pf
        Cd = float(drag_interp(Mach))
    
    # Compute the forces
    
    # Gravity force
    gravityF = gravity(x, y, z) * mass
    

    # Aerodynamic force
    
    if rho_alt < 1e-6:
        aeroF = np.zeros(3)
    else:
        aeroF = -0.5 * min(rho_alt, 1.293) * Area * abs(V) * Cd * np.array([velx, vely, velz]) # Aerodynamic drag force



    # Thrust magnitude and mass flow
    T_mag = 0.0
    mdot = 0.0

    if t <= t_burn_boosters:
        T_mag = T_mag_boosters
        mdot = -mass_flow_boosters
    else:
        T_mag = 0.0
        mdot = 0.0

    # Compute thrust vector in ECEF using local reference frame
    if T_mag > 0:
        r_pos = np.sqrt(x**2 + y**2 + z**2)
        r_hat = np.array([x, y, z]) / r_pos

        if t <= t_vertical:
            # Vertical flight — thrust along local vertical (radial direction)
            thrustF = T_mag * r_hat
        elif t > t_vertical and t <= t_vertical + 1e-3:
            k_hat = np.array([0.0, 0.0, 1.0])  # Earth's rotation axis
            east_raw = np.cross(k_hat, r_hat)
            east_norm = np.linalg.norm(east_raw)
            if east_norm > 1e-10:
                east_hat = east_raw / east_norm
            else:
                east_hat = np.array([1.0, 0.0, 0.0])
            north_hat = np.cross(r_hat, east_hat)

            # Kick: small angle from vertical in the azimuth direction
            thrust_dir = (np.cos(kick_angle) * r_hat +
                          np.sin(kick_angle) * (
                              np.sin(Az_rad) * east_hat +
                              np.cos(Az_rad) * north_hat))
            thrustF = T_mag * thrust_dir
        

        else:
            # Phase 3: Pitch program below atmosphere, gravity turn above
            altitude = r_pos - 6371000.0
            k_hat_local = np.array([0.0, 0.0, 1.0])
            east_raw = np.cross(k_hat_local, r_hat)
            east_norm = np.linalg.norm(east_raw)
            if east_norm > 1e-10:
                east_hat = east_raw / east_norm
            else:
                east_hat = np.array([1.0, 0.0, 0.0])
            north_hat = np.cross(r_hat, east_hat)

            if altitude < 80000.0:
                thrust_dir = (np.cos(kick_angle) * r_hat +
                              np.sin(kick_angle) * (
                                  np.sin(Az_rad) * east_hat +
                                  np.cos(Az_rad) * north_hat))
            else:
                if V > 1e-6:
                    thrust_dir = v_vec / V
                else:
                    thrust_dir = r_hat
            thrustF = T_mag * thrust_dir
    else:
        thrustF = np.zeros(3)


    Forces = gravityF + aeroF + thrustF

    if mass > 0:
        vdot = Forces / mass
    else:
        vdot = np.zeros(3)
        mdot = 0.0
    

    statedot = np.array([xdot, ydot, zdot, vdot[0], vdot[1], vdot[2], mdot])
    return statedot


def extract_results(stateout):
    x = stateout[:, 0]
    y = stateout[:, 1]
    z = stateout[:, 2]
    velx = stateout[:, 3]
    vely = stateout[:, 4]
    velz = stateout[:, 5]
    mass = stateout[:, 6]

    return x, y, z, velx, vely, velz, mass

############################################################################

def integration_stages(stateinitial, tout, stages_info, boosters_info, Area_pf, Area_bf, Cd_of_crosflow_cylinder, t_vertical, Az_rad, stage_index, kick_angle_deg, stages_count, rocket_has_boosters):
    t_burn_stages = stages_info[0]
    t_burn_boosters = boosters_info[0]
    simulation_time = 200000
    m_construction_each_boosters = boosters_info[3]

    if stage_index == 0 and rocket_has_boosters:
        time_with_boosters = np.linspace(0, t_burn_boosters, 10000)
        time_1st_stage_without_boosters = np.linspace(t_burn_boosters, t_burn_stages[stage_index], 10000)
        tout_propelled = np.concatenate((time_with_boosters, time_1st_stage_without_boosters))

        stateout_with_boosters = odeint(Derivatives_with_boosters, stateinitial, time_with_boosters, args=(stages_info, boosters_info, Area_pf, Area_bf, Cd_of_crosflow_cylinder, t_vertical, Az_rad, stage_index, kick_angle_deg,))
        state_initial_without_boosters = stateout_with_boosters[-1].copy()
        mass_after_booster_separation = state_initial_without_boosters[6] - sum(m_construction_each_boosters)
        state_initial_without_boosters[6] = mass_after_booster_separation
        stateout_without_boosters = odeint(Derivatives_propelled, state_initial_without_boosters, time_1st_stage_without_boosters, args=(stages_info, boosters_info, Area_pf, Area_bf, Cd_of_crosflow_cylinder, t_vertical, Az_rad, stage_index, kick_angle_deg,))
        stateout_propelled = np.concatenate((stateout_with_boosters, stateout_without_boosters))

        time_1st_stage_balistic = np.linspace(t_burn_stages[stage_index], simulation_time, 10000)
        tout = np.concatenate((tout_propelled, time_1st_stage_balistic))

        state_initial_1st_stage_balistic = stateout_propelled[-1].copy()
        
        stateout_balistic = odeint(Derivatives_balistic, state_initial_1st_stage_balistic, time_1st_stage_balistic, args=(Area_pf, Area_bf, Cd_of_crosflow_cylinder,))
        stateout = np.concatenate((stateout_propelled, stateout_balistic))
    
    elif stage_index == 0 and not rocket_has_boosters:
        tout_propelled = np.linspace(0, t_burn_stages[stage_index], 10000)
        time_stage_balistic = np.linspace(t_burn_stages[stage_index], simulation_time, 10000)
        tout = np.concatenate((tout_propelled, time_stage_balistic))

        stateout_propelled = odeint(Derivatives_propelled, stateinitial, tout_propelled, args=(stages_info, boosters_info, Area_pf, Area_bf, Cd_of_crosflow_cylinder, t_vertical, Az_rad, stage_index, kick_angle_deg,))
        state_initial_stage_balistic = stateout_propelled[-1].copy()
        
        stateout_balistic = odeint(Derivatives_balistic, state_initial_stage_balistic, time_stage_balistic, args=(Area_pf, Area_bf, Cd_of_crosflow_cylinder,))
        stateout = np.concatenate((stateout_propelled, stateout_balistic))

    elif stage_index == stages_count - 1:
        # Last stage: partial ascent burn, coast to apogee, circularize
        fuel_reserve_fraction = 0.06534  # 6.534% of fuel reserved for circularization
        
        t_burn_stages_info = stages_info[0]
        mass_flow_stages = stages_info[2]
        m_construction_stages = stages_info[3]
        
        m_fuel_total = stateinitial[6] - m_construction_stages[stage_index]
        m_fuel_ascent = m_fuel_total * (1.0 - fuel_reserve_fraction)
        m_fuel_circ = m_fuel_total * fuel_reserve_fraction
        
        # Ascent burn time (only burn the ascent fuel fraction)
        t_burn_ascent = m_fuel_ascent / mass_flow_stages[stage_index]
        
        t_start = sum(t_burn_stages_info[:stage_index])
        t_end_ascent = t_start + t_burn_ascent
        t_end_full = t_start + t_burn_stages_info[stage_index]
        
        # --- Phase 1: Ascent burn (gravity turn) ---
        # Temporarily increase construction mass to stop burn early
        original_m_construction = m_construction_stages[stage_index]
        m_construction_stages[stage_index] = original_m_construction + m_fuel_circ
        
        tout_ascent = np.linspace(t_start, t_end_ascent, 10000)
        stateout_ascent = odeint(Derivatives_propelled, stateinitial, tout_ascent,
                                 args=(stages_info, boosters_info, Area_pf, Area_bf,
                                       Cd_of_crosflow_cylinder, t_vertical, Az_rad,
                                       stage_index, kick_angle_deg,))
        
        # Restore original construction mass for circularization
        m_construction_stages[stage_index] = original_m_construction
        
        # --- Phase 2: Coast to apogee ---
        state_after_ascent = stateout_ascent[-1].copy()
        
        # Coast for a long time to find apogee
        coast_duration = 10000.0  # seconds max coast
        tout_coast = np.linspace(t_end_ascent, t_end_ascent + coast_duration, 20000)
        stateout_coast = odeint(Derivatives_balistic, state_after_ascent, tout_coast,
                                args=(Area_pf, Area_bf, Cd_of_crosflow_cylinder,))
        
        # Find apogee: where radial velocity changes sign (positive to negative)
        r_dot = np.array([np.dot(stateout_coast[j, :3], stateout_coast[j, 3:6]) /
                          np.linalg.norm(stateout_coast[j, :3])
                          for j in range(len(stateout_coast))])
        
        # Find first zero-crossing (positive to negative = apogee)
        apogee_idx = None
        for j in range(1, len(r_dot)):
            if r_dot[j-1] > 0 and r_dot[j] <= 0:
                apogee_idx = j
                break
        
        if apogee_idx is None:
            # No apogee found (escape trajectory or already descending)
            # Just use end of coast
            apogee_idx = len(stateout_coast) - 1
        
        # Trim coast to apogee
        tout_coast_to_apogee = tout_coast[:apogee_idx+1]
        stateout_coast_to_apogee = stateout_coast[:apogee_idx+1]
        
        # --- Phase 3: Circularization burn at apogee ---
        state_at_apogee = stateout_coast_to_apogee[-1].copy()
        t_apogee = tout_coast_to_apogee[-1]
        
        # Burn time for remaining fuel
        t_burn_circ = m_fuel_circ / mass_flow_stages[stage_index]
        tout_circ = np.linspace(t_apogee, t_apogee + t_burn_circ, 5000)
        
        # Use Derivatives_propelled for circularization — above 80km it thrusts prograde
        stateout_circ = odeint(Derivatives_propelled, state_at_apogee, tout_circ,
                               args=(stages_info, boosters_info, Area_pf, Area_bf,
                                     Cd_of_crosflow_cylinder, t_vertical, Az_rad,
                                     stage_index, kick_angle_deg,))
        
        # --- Phase 4: Coast after circularization ---
        state_after_circ = stateout_circ[-1].copy()
        t_end_circ = tout_circ[-1]
        tout_final_coast = np.linspace(t_end_circ, simulation_time, 5000)
        stateout_final_coast = odeint(Derivatives_balistic, state_after_circ, tout_final_coast,
                                       args=(Area_pf, Area_bf, Cd_of_crosflow_cylinder,))
        
        # Concatenate all phases
        tout_propelled = np.concatenate((tout_ascent, tout_circ))
        stateout_propelled = np.concatenate((stateout_ascent, stateout_circ))
        
        tout = np.concatenate((tout_ascent, tout_coast_to_apogee, tout_circ, tout_final_coast))
        stateout = np.concatenate((stateout_ascent, stateout_coast_to_apogee, stateout_circ, stateout_final_coast))
    else:
        t_start = sum(t_burn_stages[:stage_index])
        t_end = t_start + t_burn_stages[stage_index]
        tout_propelled = np.linspace(t_start, t_end, 10000)
        time_stage_balistic = np.linspace(t_end, simulation_time, 10000)
        tout = np.concatenate((tout_propelled, time_stage_balistic))

        stateout_propelled = odeint(Derivatives_propelled, stateinitial, tout_propelled, args=(stages_info, boosters_info, Area_pf, Area_bf, Cd_of_crosflow_cylinder, t_vertical, Az_rad, stage_index, kick_angle_deg,))
        state_initial_balistic = stateout_propelled[-1].copy()
        
        stateout_balistic = odeint(Derivatives_balistic, state_initial_balistic, time_stage_balistic, args=(Area_pf, Area_bf, Cd_of_crosflow_cylinder,))
        stateout = np.concatenate((stateout_propelled, stateout_balistic))


    

    return tout, stateout, tout_propelled, stateout_propelled
##############################################################################

def integration_boosters(stateinitial, tout, t_burn_boosters, T_mag_boosters, mass_flow_boosters, m_construction_each_boosters, Area_pf, Area_bf, Cd_of_crosflow_cylinder, t_vertical, Az_rad, kick_angle_deg):
    simulation_time = 3000
    tout = np.linspace(0, simulation_time, 10000)
    stateout = odeint(Derivatives_boosters, stateinitial, tout, args=(t_burn_boosters, T_mag_boosters, mass_flow_boosters, m_construction_each_boosters, Area_pf, Area_bf, Cd_of_crosflow_cylinder, t_vertical, Az_rad, kick_angle_deg))

    tout_burn = np.linspace(0, t_burn_boosters, 10000)
    stateout_burn = odeint(Derivatives_boosters, stateinitial, tout_burn, args=(t_burn_boosters, T_mag_boosters, mass_flow_boosters, m_construction_each_boosters, Area_pf, Area_bf, Cd_of_crosflow_cylinder, t_vertical, Az_rad, kick_angle_deg))

    return tout, stateout, tout_burn, stateout_burn

##############################################################################

def compute_corrected_Azimuth(launch_lat, t_o_i, t_o_a):
    t_o_a = t_o_a * 1000  # km -> m
    t_day = 24*60*60  # [s]

    V_orb = np.sqrt(G * Mplanet / t_o_a)

    Vpad = ((2 * np.pi * Rplanet) / t_day) * np.cos(np.deg2rad(launch_lat))
    Azimuth = np.arcsin(np.cos(np.deg2rad(t_o_i)) / np.cos(np.deg2rad(launch_lat)))

    VL_north = V_orb * np.cos(Azimuth)
    VL_east = V_orb * np.sin(Azimuth) - Vpad
    corrected_Azimuth = np.arctan2(VL_east, VL_north)

    return corrected_Azimuth, Vpad

def utc_to_julian_date(launch_date, launch_time):
    date_parts = launch_date.split("-")
    time_parts = launch_time.split(":")
    year = int(date_parts[0])
    month = int(date_parts[1])
    day = int(date_parts[2])
    hour = int(time_parts[0])
    minute = int(time_parts[1])
    second = int(time_parts[2])

    Df = day + (hour / 24) + (minute / 1440) + (second / 86400)

    if month <= 2:
        year_prim = year - 1
        month_prim = month + 12
    else:
        year_prim = year
        month_prim = month
    
    A = int(year_prim / 100)
    B = 2 - A + int(A / 4)

    # Julian Date
    JD = int(365.25 * (year_prim + 4716)) + int(30.6001 * (month_prim + 1)) + Df + B - 1524.5

    return JD

def compute_gmst(jd_utc):
    
    JD_J2000 = 2451545.0  # Julian Date of J2000 epoch
    
    # Centuries since J2000
    T = (jd_utc - JD_J2000) / 36525.0
    
    # GMST in degrees
    theta_GMST_deg = (280.46061837 +
                      360.98564736629 * (jd_utc - JD_J2000) +
                      0.000387933 * T**2 -
                      T**3 / 38710000.0)
    
    # Reduce to [0, 360)
    theta_GMST_deg = theta_GMST_deg % 360.0
    
    return np.deg2rad(theta_GMST_deg)

def eci_greenwich_to_j2000(state, gmst_rad):
    cos_g = np.cos(gmst_rad)
    sin_g = np.sin(gmst_rad)
    
    # Rotation matrix Rz(-GMST): from Greenwich ECI to J2000
    # r_J2000 = Rz(-GMST) * r_greenwich
    x_j2000 = state[0] * cos_g + state[1] * sin_g
    y_j2000 = -state[0] * sin_g + state[1] * cos_g
    z_j2000 = state[2]
    
    vx_j2000 = state[3] * cos_g + state[4] * sin_g
    vy_j2000 = -state[3] * sin_g + state[4] * cos_g
    vz_j2000 = state[5]
    
    result = np.array([x_j2000, y_j2000, z_j2000, vx_j2000, vy_j2000, vz_j2000])
    
    # Preserve mass if present
    if len(state) > 6:
        result = np.append(result, state[6:])
    
    return result

def cartesian_to_keplerian(state):
    """
    Convert Cartesian state vector to Keplerian orbital elements.
    """
    mu=G * Mplanet
    r_vec = state[0:3]
    v_vec = state[3:6]
    
    r = np.linalg.norm(r_vec)
    v = np.linalg.norm(v_vec)
    
    # Specific angular momentum
    h_vec = np.cross(r_vec, v_vec)
    h = np.linalg.norm(h_vec)
    
    # Node vector
    k_hat = np.array([0.0, 0.0, 1.0])
    n_vec = np.cross(k_hat, h_vec)
    n = np.linalg.norm(n_vec)
    
    # Eccentricity vector
    # e_vec = (np.cross(v_vec, h_vec) / mu) - (r_vec / r)
    e_vec = 1 / mu * (np.cross(v_vec, h_vec) - mu * r_vec / r)
    e = np.linalg.norm(e_vec)
    
    # Specific orbital energy
    epsilon = v**2 / 2.0 - mu / r
    
    # Semi-major axis
    if abs(1 - e) > 1e-10:  # not parabolic
        # a = (-mu / (2 * epsilon)) / 1000
        a = ((2 / r) - (v**2 / mu))**(-1) / 1000
    else:
        a = np.inf
    
    # Inclination
    i = np.arccos(np.clip(h_vec[2] / h, -1, 1))
    i = np.rad2deg(i)
    
    # RAAN
    if n > 1e-10:
        RAAN = np.arccos(np.clip(n_vec[0] / n, -1, 1))
        if n_vec[1] < 0:
            RAAN = 2 * np.pi - RAAN
    else:
        RAAN = 0.0
    
    # Argument of periapsis
    if n > 1e-10 and e > 1e-10:
        omega = np.arccos(np.clip(np.dot(n_vec, e_vec) / (n * e), -1, 1))
        if e_vec[2] < 0:
            omega = 2 * np.pi - omega
    else:
        omega = 0.0
    
    # True anomaly
    if e > 1e-10:
        nu = np.arccos(np.clip(np.dot(e_vec, r_vec) / (e * r), -1, 1))
        if np.dot(r_vec, v_vec) < 0:
            nu = 2 * np.pi - nu
    else:
        # Circular orbit: use argument of latitude
        if n > 1e-10:
            nu = np.arccos(np.clip(np.dot(n_vec, r_vec) / (n * r), -1, 1))
            if r_vec[2] < 0:
                nu = 2 * np.pi - nu
        else:
            nu = 0.0
    
    return {
        'a': a,
        'e': e,
        'i': i,
        'RAAN': RAAN,
        'omega': omega,
        'nu': nu
    }
############################################################################
def geodetic_to_cartesian_WGS84(latitude_deg, longitude_deg, altitude):
    latitude = np.deg2rad(latitude_deg)
    longitude = np.deg2rad(longitude_deg)
    
    a = 6378137 # [m]
    e2 = 0.00669437999014
    N = a / (np.sqrt(1 - e2 * np.sin(latitude)**2))

    x = (N + altitude) * np.cos(latitude) * np.cos(longitude)
    y = (N + altitude) * np.cos(latitude) * np.sin(longitude)
    z = (N*(1 - e2) + altitude) * np.sin(latitude)

    return x, y, z
############################################################################

############################################################################
def parameters_of_stages(input_mode, data_list, m_payload_without_boosters, payload_mass_ratio_total, rocket_type, stages_count):
    """Return Ve, mass_flow, m0, m_prop, Vf_id, Lambda for each stage depending on input mode and rocket type"""
    
    if input_mode == "Start mass & Propellant":
        m0 = [stage["Start Mass (kg)"] for stage in data_list]
        m_prop = [stage["Propellant (kg)"] for stage in data_list]
        Ve = [stage["Ve (m/s)"] for stage in data_list]
        mass_flow = [stage["Mass flow (kg/s)"] for stage in data_list]

        Lambda = [m0[i] / (m0[i] - m_prop[i]) for i in range(len(m0))]
        Vf_id = [Ve[i] * np.log(Lambda[i]) for i in range(len(m0))]
        

    else:
        eps = [stage["EPS"] for stage in data_list]
        mass_flow = [stage["Mass_flow (kg/s)"] for stage in data_list]
        Ve = [stage["Ve (m/s)"] for stage in data_list]
        diameter_stages = [stage["Diameter (m)"] for stage in data_list]
        height_stages = [stage["Height (m)"] for stage in data_list]

        if rocket_type == "Optimal":
            m0, m_prop, Vf_id, Lambda = optimal_rocket_parameters(eps, payload_mass_ratio_total, Ve, m_payload_without_boosters, stages_count)
        else:
            m0, m_prop, Vf_id, Lambda = non_optimal_rocket_parameters(eps, payload_mass_ratio_total, Ve, m_payload_without_boosters, stages_count)

    
    return Ve, mass_flow, diameter_stages, height_stages, m0, m_prop, Vf_id, Lambda
############################################################################


############################################################################
def parameters_of_boosters(input_mode, data_list, m_payload_without_boosters, m_payload_with_boosters, Vf_id_stages, m0_stages, m_prop_stages, stage_count, booster_count, Ve_stages, mass_flow_stages, t_burn_ratio):
    if input_mode == "Start mass & Propellant":
        # placeholder
        pass       

    else:
        eps = [booster["EPS"] for booster in data_list]
        mass_flow_boosters = [booster["Mass_flow (kg/s)"] for booster in data_list]
        Ve_boosters = [booster["Ve (m/s)"] for booster in data_list]
        diameter_boosters = [booster["Diameter (m)"] for booster in data_list]
        height_boosters = [booster["Height (m)"] for booster in data_list]

        delta_m_payload = m_payload_with_boosters - m_payload_without_boosters

        new_m0_stages = [m0_stages[i] + delta_m_payload for i in range(stage_count)]
        phi_with_boosters = [None] * (stage_count - 1)
        Lambda_with_boosters = [None] * (stage_count - 1)

        for i in range(stage_count -1, 0, -1):
            phi_with_boosters[i-1] = m_prop_stages[i] / new_m0_stages[i]
            Lambda_with_boosters[i-1] = 1 / (1 - phi_with_boosters[i-1])

        Vf_id_with_boosters = [Ve_stages[i] * np.log(Lambda_with_boosters[i]) for i in range(stage_count - 1)]

        Vf_id_first_stage_with_boosters = sum(Vf_id_stages) - sum(Vf_id_with_boosters)

        phi_first_stage_with_boosters = m_prop_stages[0] / (m0_stages[0] + delta_m_payload)

        Lambda_double_prim = (1 - phi_first_stage_with_boosters * t_burn_ratio) / (1 - phi_first_stage_with_boosters)
        deltaV_double_prim = Ve_stages[0] * np.log(Lambda_double_prim)
        deltaV_prim = Vf_id_first_stage_with_boosters - deltaV_double_prim
        
        mass_flow_boosters_total = sum(mass_flow_boosters)
        Ve_equivalent = Ve_stages[0] - ((1 / (1 + mass_flow_stages[0] / mass_flow_boosters_total)) * (Ve_stages[0] - Ve_boosters[0]))
        
        Lambda_prim = np.exp(deltaV_prim / Ve_equivalent)
        
        y = 1 - (1/Lambda_prim)
        m_prop_boosters = ((y * (m0_stages[0] + delta_m_payload)) - m_prop_stages[0] * t_burn_ratio) / (1 - (y / (1 - eps[0])))
        m_construction_boosters = (eps[0] * m_prop_boosters) / (1 - eps[0])

        m_prop_each_boosters = [m_prop_boosters / booster_count for i in range(booster_count)]
        m_construction_each_boosters = [m_construction_boosters / booster_count for i in range(booster_count)]
        m0_each_boosters = [m_prop_each_boosters[i] + m_construction_each_boosters[i] for i in range(booster_count)]


        new_m0_stages[0] = m0_stages[0] + delta_m_payload + sum(m0_each_boosters)

        return Ve_boosters, mass_flow_boosters, diameter_boosters, height_boosters, m0_each_boosters, m_construction_each_boosters, m_prop_each_boosters, new_m0_stages


############################################################################


############################################################################
def optimal_rocket_parameters(eps, payload_mass_ratio_total, Ve, m_payload_without_boosters, stages_count):
    """Calculate parameters for an optimal rocket configuration."""
    mu = calculate_optimal_mu(payload_mass_ratio_total, eps, Ve)

    lambda_optimal = [None] * stages_count
    for i in range(stages_count):
        lambda_optimal[i] = (mu * eps[i]) / ((Ve[i] - mu) * (1 - eps[i]))
    
    phi = [None] * stages_count
    for i in range(stages_count):
        phi[i] = (1 - eps[i]) * (1 - lambda_optimal[i])
   
    Lambda = [None] * stages_count
    for i in range(stages_count):
        Lambda[i] = 1 / (1 - phi[i])
    
    Vf_id = [None] * stages_count
    for i in range(stages_count):
        Vf_id[i] = Ve[i] * np.log(Lambda[i])

    m0 = [None] * stages_count
    m_prop = [None] * stages_count
    for i in range(stages_count - 1, -1, -1):
        if i == stages_count - 1: # last section
            m0[i] = m_payload_without_boosters / lambda_optimal[i] # last section mass
            m_prop[i] = m0[i] * phi[i]
        else:
            m0[i] = m0[i+1] / lambda_optimal[i]
            m_prop[i] = m0[i] * phi[i]
    
    return m0, m_prop, Vf_id, Lambda
############################################################################
      
############################################################################
def non_optimal_rocket_parameters(eps, payload_mass_ratio_total, Ve, m_payload_without_boosters, stages_count):
    lambda_of_non_optimal_rocket = payload_mass_ratio_total**(1 / stages_count)
    
    phi = [None] * stages_count
    for i in range(stages_count):
        phi[i] = (1 - eps[i]) * (1 - lambda_of_non_optimal_rocket)
   
    Lambda = [None] * stages_count
    for i in range(stages_count):
        Lambda[i] = 1 / (1 - phi[i])
    
    Vf_id = [None] * stages_count
    for i in range(stages_count):
        Vf_id[i] = Ve[i] * np.log(Lambda[i])

    m0 = [None] * stages_count
    m_prop = [None] * stages_count

    for i in range(stages_count - 1, -1, -1):
        if i == stages_count - 1: # last section
            m0[i] = m_payload_without_boosters / lambda_of_non_optimal_rocket # last section mass
            m_prop[i] = m0[i] * phi[-1]
        else:
            m0[i] = m0[i+1] / lambda_of_non_optimal_rocket
            m_prop[i] = m0[i] * phi[i]
    
    return m0, m_prop, Vf_id, Lambda
############################################################################

############################################################################
def calculate_optimal_mu(lambda_total, eps_list, Ve_list):
    """Calculate the optimal mu for given lambda_total, eps_list, and Ve_list."""
    
    if len(eps_list) != len(Ve_list):
        raise ValueError("EPS and Ve must have the same length.")
    
    min_Ve = min(Ve_list)

    target_log = np.log(lambda_total)
    
    
    def equation(mu):
        current_log_sum = 0
        
        for eps, Ve in zip(eps_list, Ve_list):
            term = np.log(mu) + np.log(eps) - np.log(Ve - mu) - np.log(1 - eps)
            current_log_sum += term
        return current_log_sum - target_log
    
    try:
        mu_optimal = brentq(equation, 1e-9, min_Ve - 1e-5)
        return mu_optimal
    except ValueError:
        raise ValueError("No solution found for the given parameters.")
############################################################################



############################################################################
############################################################################
# Functions for Result Page
def test_plot(stages_count, velmag_stages, tout_stages):
    st.markdown("### Select stages")

    stage_visibility = []
    for i in range(stages_count):
        checked = st.checkbox(f"Stage {i + 1}", value=True, key=f"stage_{i}")
        stage_visibility.append(checked)

    # ---------- Plot ----------
    fig = go.Figure()

    for i in range(stages_count):
        if stage_visibility[i]:
            fig.add_trace(
                go.Scatter(
                    x=tout_stages[i],
                    y=velmag_stages[i],
                    mode="lines",
                    name=f"Stage {i + 1}"
                )
            )

    fig.update_layout(
        xaxis_title="Time (s)",
        yaxis_title="Velocity magnitude (m/s)",
        legend_title="Stages",
        template="plotly_white",
        height=600
    )

    st.plotly_chart(fig, use_container_width=True)
def plot_3d_orbit(data):
    """Function to plot 3D Orbit"""
    pass
    

def plot_velocity_vs_time(data):
    """Function to plot velocity vs time"""
    pass # placeholder for velocity vs time plotting code

def plot_altitude_vs_time(data):
    """Function to plot altitude vs time"""
    pass # placeholder for altitude vs time plotting code

def plot_mass_vs_time(data):
    """Function to plot mass vs time"""
    pass # placeholder for mass vs time plotting code

def plot_density_vs_altitude(data):
    """Function to plot density vs altitude"""
    pass # placeholder for density vs altitude plotting code

def plot_temperature_profile():
    """Function to plot temperature profile"""
    st.write("Plotting Temperature Profile...")
    fig = px.line(
        temp_data, 
        x='Temperature (K)', 
        y='Altitude (km)',
    )

    st.plotly_chart(fig)

def plot_drag_coefficient_vs_mach():
    """Function to plot drag coefficient vs Mach"""
    st.write("Plotting Drag Coefficient vs Mach")
    fig = px.line(
        cd_data,
        x='Mach',
        y='Cd',
    )
    st.plotly_chart(fig)
############################################################################
############################################################################
