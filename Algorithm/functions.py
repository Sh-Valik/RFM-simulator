import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import brentq
import os
import streamlit as st
import plotly.express as px
import pandas as pd

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
drag_interp = interp1d(Mach_ref, Cd_ref, kind='linear')

new_Mach = np.linspace(min(Mach_ref), max(Mach_ref), 200)
new_Cd = drag_interp(new_Mach)
cd_data = pd.DataFrame({'Mach': new_Mach, 'Cd': new_Cd})
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

        if rocket_type == "Optimal":
            m0, m_prop, Vf_id, Lambda = optimal_rocket_parameters(eps, payload_mass_ratio_total, Ve, m_payload_without_boosters, stages_count)
        else:
            m0, m_prop, Vf_id, Lambda = non_optimal_rocket_parameters(eps, payload_mass_ratio_total, Ve, m_payload_without_boosters, stages_count)

    
    return Ve, mass_flow, m0, m_prop, Vf_id, Lambda
############################################################################


############################################################################
def parameters_of_boosters(input_mode, data_list, m_payload_without_boosters, m_payload_with_boosters, Vf_id_stages, m0_stages, m_prop_stages, stage_count, booster_count, Ve_stages, t_burn_ratio):
    if input_mode == "Start mass & Propellant":
        # placeholder
        pass       

    else:
        eps = [booster["EPS"] for booster in data_list]
        mass_flow_boosters = [booster["Mass_flow (kg/s)"] for booster in data_list]
        Ve_boosters = [booster["Ve (m/s)"] for booster in data_list]

        delta_m_payload = m_payload_with_boosters - m_payload_without_boosters

        new_m0_stages = [m0_stages[i] + delta_m_payload for i in range(stage_count -1, 0, -1)]
        phi_with_boosters = [None] * (stage_count - 1)
        Lambda_with_boosters = [None] * (stage_count - 1)
        for i in range(stage_count -1, 0, -1):
            phi_with_boosters[i] = m_prop_stages[i] / new_m0_stages[i]
            Lambda_with_boosters[i] = 1 / (1 - phi_with_boosters[i])

        Vf_id_with_boosters = [Ve_stages[i] * np.log(Lambda_with_boosters[i]) for i in range(stage_count - 1)]

        Vf_id_first_stage_with_boosters = sum(Vf_id_stages) - sum(Vf_id_with_boosters)

        phi_first_stage_with_boosters = m_prop_stages[0] / (m0_stages[0] + delta_m_payload)

        Lambda_double_prim = (1 - phi_first_stage_with_boosters * t_burn_ratio) / (1 - phi_first_stage_with_boosters)
        deltaV_double_prim = Ve_stages[0] * np.log(Lambda_double_prim)
        deltaV_prim = Vf_id_first_stage_with_boosters - deltaV_double_prim
        
        Lambda_prim = [np.exp(deltaV_prim / Ve_boosters[i]) for i in range(booster_count)]
        mps_over_mp1 = [mps_over_mp1(Lambda_prim[i], phi_first_stage_with_boosters, eps[i], t_burn_ratio) for i in range(booster_count)]
        m_prop_boosters = [mps_over_mp1[i] * m0_stages[0] for i in range(booster_count)]

        m_construction_boosters = [(eps[i] * m_prop_boosters[i]) / (1 - eps[i]) for i in range(booster_count)]
        m0_boosters = [m_prop_boosters[i] + m_construction_boosters[i] for i in range(booster_count)]


        new_m0_stages = np.insert(new_m0_stages, 0, m0_stages[0] + delta_m_payload + sum(m0_boosters))

        return Ve_boosters, mass_flow_boosters, m0_boosters, m_prop_boosters, m_construction_boosters, Lambda_with_boosters, new_m0_stages


############################################################################
def mps_over_mp1(Lambda, phi, eps, tbs_over_tb1):
    A = (1 - 1 / Lambda_) * phi / (1 - epsilon)
    B = (1 - 1 / Lambda_) - phi * tbs_over_tb1
    C = -phi

    D = B**2 - 4 * A * C
    x1 = (-B + math.sqrt(D)) / (2 * A)
    x2 = (-B - math.sqrt(D)) / (2 * A)

    if x1 > 0 and x2 > 0:
        return min(x1, x2)
    elif x1 > 0:
        return x1
    elif x2 > 0:
        return x2
    else:
        return None

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
def plot_3d_orbit(data):
    """Function to plot 3D Orbit"""
    pass # placeholder for 3D orbit plotting code

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
