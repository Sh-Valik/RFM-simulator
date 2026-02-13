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
def propulsion(t, t_burn, T_mag_stages, m_flow_stages, t_burn_boosters, T_mag_boosters, m_flow_boosters, T_mag_, m_flow_, has_boosters=True):
    """ Thrust magnitude and mass flow rate computation """
    if has_boosters:
        if t <= t_burn_boosters:
            T_mag = T_mag_stages[0] + sum(T_mag_boosters)
            m_flow = m_flow_stages[0] + sum(m_flow_boosters)
        else:
            T_mag = T_mag_stages[0]
            m_flow = m_flow_stages[0]
    else:
        T_mag = T_mag_
        m_flow = m_flow_

    if t <= t_burn:
        mdot = - m_flow  # mass flow rate [kg/s]
    else:
        T_mag = 0.0
        mdot = 0.0

    return T_mag, mdot


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
def Derivatives(state, t, stages_info, boosters_info, m_construction, Area_pf, Area_bf, Cd_of_crosflow_cylinder, t_vertical, theta_rad, Az_rad, T_mag_, m_flow_, t_burn, has_boosters=True):
    """ Computes the state derivatives """

    # Unpack stages_info
    
    T_mag_stages = stages_info[1]
    mass_flow_stages = stages_info[2]
    
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
    if t <= t_burn:
        mass = state[6]
    else:
        mass = m_construction

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
    apogee_reached = False
    if velz < 0:
        apogee_reached = True
    
    if apogee_reached:
        Area = Area_bf
        Cd = Cd_of_crosflow_cylinder
    else:
        Area = Area_pf
        Cd = drag_interp(Mach)
    
    # Compute the forces
    
    # Gravity force
    gravityF = gravity(x, y, z) * mass
    

    # Aerodynamic force
    
    if rho_alt < 1e-6:
        aeroF = np.zeros(3)
    else:
        aeroF = -0.5 * min(rho_alt, 1.293) * Area * abs(V) * Cd * np.array([velx, vely, velz]) # Aerodynamic drag force



    # Thrust magnitude and mass flow
    T_mag, mdot = propulsion(t, t_burn, T_mag_stages, mass_flow_stages, t_burn_boosters, T_mag_boosters, mass_flow_boosters, T_mag_, m_flow_, has_boosters)

    # Compute thrust vector in ECEF using local reference frame
    if T_mag > 0:
        r_pos = np.sqrt(x**2 + y**2 + z**2)
        r_hat = np.array([x, y, z]) / r_pos

        if t <= t_vertical:
            # Vertical flight — thrust along local vertical (radial direction)
            thrustF = T_mag * r_hat
        else:
            # Gravity turn — thrust at angle theta from local vertical
            # Build local ENU (East-North-Up) frame
            k_hat = np.array([0.0, 0.0, 1.0])  # Earth's rotation axis
            east_raw = np.cross(k_hat, r_hat)
            east_norm = np.linalg.norm(east_raw)
            if east_norm > 1e-10:
                east_hat = east_raw / east_norm
            else:
                east_hat = np.array([1.0, 0.0, 0.0])
            north_hat = np.cross(r_hat, east_hat)

            thrust_dir = (np.cos(theta_rad) * r_hat +
                          np.sin(theta_rad) * (np.sin(Az_rad) * east_hat + np.cos(Az_rad) * north_hat))
            thrustF = T_mag * thrust_dir
    else:
        thrustF = np.zeros(3)


    Forces = gravityF + aeroF + thrustF

    q = q_dynamic_pressure(min(rho_alt, 1.293), V)
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

############################################################################
def integration_stop_event(t, state):
    """
    Event function for solve_ivp. 
    Returns 0 when the rocket hits the ground (Altitude = 0).
    """
    x, y, z = state[0], state[1], state[2]
    # Текущее расстояние от центра планеты
    r_current = np.sqrt(x**2 + y**2 + z**2)
    
    # Событие срабатывает, когда высота становится 0
    return r_current - Rplanet
    
############################################################################

############################################################################
# def integration(stateinitial, stages_info, boosters_info, m_construction, Area_pf, Area_bf, Cd_of_crosflow_cylinder, t_vertical, theta_rad, Az_rad, t_burn, T_mag_, m_flow_, has_boosters):
#     max_simulation_time = 3000
#     t_span = (0, max_simulation_time)


#     fun = lambda t, y: Derivatives(y, t, stages_info, boosters_info, m_construction, Area_pf, Area_bf, Cd_of_crosflow_cylinder, t_vertical, theta_rad, Az_rad, T_mag_, m_flow_, t_burn, has_boosters)


#     stop_event = lambda t, y: integration_stop_event(t, y)
#     stop_event.terminal = True
#     stop_event.direction = -1

#     sol = solve_ivp(fun, t_span, stateinitial, events=stop_event, method='RK45', rtol=1e-6, max_step=0.5)

#     stateout = sol.y.T
#     tout = sol.t


#     tout_burn = np.linspace(0, t_burn, 10000)
#     stateout_burn = odeint(Derivatives, stateinitial, tout_burn, args=(stages_info, boosters_info, m_construction, Area_pf, Area_bf, Cd_of_crosflow_cylinder, t_vertical, theta_rad, Az_rad, T_mag_, m_flow_, t_burn, has_boosters,))

#     return tout, stateout, tout_burn, stateout_burn


def integration(stateinitial, tout, tout_burn, stages_info, boosters_info, m_construction, Area_pf, Area_bf, Cd_of_crosflow_cylinder, t_vertical, theta_rad, Az_rad, t_burn, T_mag_, m_flow_, has_boosters, stage_number):
    t_burn_stages = stages_info[0]
    t_burn_boosters = boosters_info[0]
    m_construction_each_boosters = boosters_info[3]
    simulation_time = 3000

    if stage_number == 0:
        time_with_boosters = np.linspace(0, t_burn_boosters, 10000)
        time_1st_stage_without_boosters = np.linspace(t_burn_boosters, simulation_time, 10000)
        
        tout[stage_number] = np.concatenate((time_with_boosters, time_1st_stage_without_boosters))

        tout_burn_without_boosters = np.linspace(t_burn_boosters, t_burn_stages[stage_number], 10000)
        tout_burn[stage_number] = np.concatenate((time_with_boosters, tout_burn_without_boosters))

        stateout_with_boosters = odeint(Derivatives, stateinitial, time_with_boosters, args=(stages_info, boosters_info, m_construction, Area_pf, Area_bf, Cd_of_crosflow_cylinder, t_vertical, theta_rad, Az_rad, T_mag_, m_flow_, t_burn, has_boosters,))
        
        last_state_with_boosters = stateout_with_boosters[-1]
        mass_after_separation = last_state_with_boosters[6] - sum(m_construction_each_boosters)

        state_initial_without_boosters = np.array([
            last_state_with_boosters[0], last_state_with_boosters[1], last_state_with_boosters[2], # x, y, z
            last_state_with_boosters[3], last_state_with_boosters[4], last_state_with_boosters[5], # velx, vely, velz
            mass_after_separation
        ])
        stateout_without_boosters = odeint(Derivatives, state_initial_without_boosters, time_1st_stage_without_boosters, args=(stages_info, boosters_info, m_construction, Area_pf, Area_bf, Cd_of_crosflow_cylinder, t_vertical, theta_rad, Az_rad, T_mag_, m_flow_, t_burn, has_boosters,))
        stateout = np.concatenate((stateout_with_boosters, stateout_without_boosters))
        
        stateout_burn_without_boosters = odeint(Derivatives, state_initial_without_boosters, tout_burn_without_boosters, args=(stages_info, boosters_info, m_construction, Area_pf, Area_bf, Cd_of_crosflow_cylinder, t_vertical, theta_rad, Az_rad, T_mag_, m_flow_, t_burn, has_boosters,))
        stateout_burn = np.concatenate((stateout_with_boosters, stateout_burn_without_boosters))
    else:
        # tout[stage_number] = np.linspace(t_burn_stages[stage_number-1], simulation_time, 10000)
        # tout_burn[stage_number] = np.linspace(t_burn_stages[stage_number-1], t_burn_stages[stage_number], 10000)
        tout[stage_number] = np.linspace(0, simulation_time, 10000)
        tout_burn[stage_number] = np.linspace(0, t_burn, 10000)

        stateout = odeint(Derivatives, stateinitial, tout[stage_number], args=(stages_info, boosters_info, m_construction, Area_pf, Area_bf, Cd_of_crosflow_cylinder, t_vertical, theta_rad, Az_rad, T_mag_, m_flow_, t_burn, has_boosters,))
        stateout_burn = odeint(Derivatives, stateinitial, tout_burn[stage_number], args=(stages_info, boosters_info, m_construction, Area_pf, Area_bf, Cd_of_crosflow_cylinder, t_vertical, theta_rad, Az_rad, T_mag_, m_flow_, t_burn, has_boosters,))

    return tout, stateout, tout_burn, stateout_burn
############################################################################


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
