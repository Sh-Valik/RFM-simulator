import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional
from scipy.interpolate import interp1d
from scipy.optimize import brentq, minimize_scalar
from scipy.integrate import solve_ivp
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
# ======================= DATA CLASSES =======================
############################################################################

@dataclass
class RocketStage:
    """Encapsulates all parameters for one rocket stage or booster."""
    T_mag: float          # thrust magnitude [N]
    m_flow: float         # mass flow rate [kg/s]
    t_burn: float         # burn time [s]
    m_construction: float # dry mass (structure) [kg]
    Area_pf: float        # cross-section area in propelled flight [m^2]
    Area_bf: float        # cross-section area in ballistic flight [m^2]


@dataclass
class FlightConfig:
    """Complete flight configuration for the simulation."""
    stages: List[RocketStage]
    boosters: Optional[List[RocketStage]]
    t_vertical: float              # duration of vertical flight [s]
    kick_angle_rad: float          # kick angle at end of vertical flight [rad] (~0.6°)
    Az_rad: float                  # launch azimuth [rad]
    Cd_crossflow: float = 1.25     # Cd for ballistic (crossflow) flight
    t_burn_boosters: float = 0.0   # booster burnout time [s]
    kick_duration: float = 2.0     # duration of kick maneuver [s]


############################################################################
# ======================= DERIVATIVES (EVENT-DRIVEN) =======================
############################################################################

def Derivatives(t, state, config: FlightConfig, stage_index: int, boosters_active: bool):
    """
    Computes the state derivatives for solve_ivp.

    Parameters
    ----------
    t : float
        Current time [s]
    state : np.ndarray
        State vector [x, y, z, vx, vy, vz, mass]
    config : FlightConfig
        Flight configuration with all stages/boosters
    stage_index : int
        Index of the current stage being simulated
    boosters_active : bool
        Whether boosters are currently attached and firing

    Returns
    -------
    np.ndarray
        State derivative [xdot, ydot, zdot, ax, ay, az, mdot]
    """
    stage = config.stages[stage_index]

    # Unpack state
    x, y, z = state[0], state[1], state[2]
    velx, vely, velz = state[3], state[4], state[5]
    mass = state[6]

    # Velocity derivatives = current velocity
    xdot = velx
    ydot = vely
    zdot = velz

    # ---- Aerodynamics ----
    V = np.sqrt(velx**2 + vely**2 + velz**2)
    rho_alt = density(x, y, z)
    Temp_local = temperature_by_altitude(x, y, z)
    loc_sound_speed = np.sqrt(1.4 * 287.05 * Temp_local)
    Mach = V / loc_sound_speed

    # Determine if apogee has been reached (radial velocity < 0)
    r_vec = np.array([x, y, z])
    v_vec = np.array([velx, vely, velz])
    v_radial = np.dot(r_vec, v_vec) / np.linalg.norm(r_vec)
    apogee_reached = v_radial < 0

    if apogee_reached:
        Area = stage.Area_bf
        Cd = config.Cd_crossflow
    else:
        Area = stage.Area_pf
        Cd = float(drag_interp(Mach))

    # Gravity force
    gravityF = gravity(x, y, z) * mass

    # Aerodynamic drag force
    if rho_alt < 1e-6:
        aeroF = np.zeros(3)
    else:
        aeroF = -0.5 * min(rho_alt, 1.293) * Area * abs(V) * Cd * np.array([velx, vely, velz])

    # ---- Propulsion (event-driven) ----
    # Determine total thrust and mass flow based on flight phase
    T_mag_total = 0.0
    mdot = 0.0

    # Stage engine: active if within burn time
    if t <= stage.t_burn:
        T_mag_total += stage.T_mag
        mdot -= stage.m_flow

    # Booster engines: active if boosters are attached and within booster burn time
    if boosters_active and config.boosters and t <= config.t_burn_boosters:
        for booster in config.boosters:
            T_mag_total += booster.T_mag
            mdot -= booster.m_flow

    # ---- Thrust direction ----
    if T_mag_total > 0:
        r_pos = np.sqrt(x**2 + y**2 + z**2)
        r_hat = np.array([x, y, z]) / r_pos

        if stage_index == 0 and t <= config.t_vertical:
            # Phase 1: Vertical flight — thrust along local vertical
            thrustF = T_mag_total * r_hat
        elif stage_index == 0 and t <= config.t_vertical + config.kick_duration:
            # Phase 2: Kick — apply kick angle for kick_duration after vertical flight ends
            # This initiates the gravity turn by tilting slightly from vertical
            k_hat = np.array([0.0, 0.0, 1.0])  # Earth's rotation axis
            east_raw = np.cross(k_hat, r_hat)
            east_norm = np.linalg.norm(east_raw)
            if east_norm > 1e-10:
                east_hat = east_raw / east_norm
            else:
                east_hat = np.array([1.0, 0.0, 0.0])
            north_hat = np.cross(r_hat, east_hat)

            # Kick: small angle from vertical in the azimuth direction
            thrust_dir = (np.cos(config.kick_angle_rad) * r_hat +
                          np.sin(config.kick_angle_rad) * (
                              np.sin(config.Az_rad) * east_hat +
                              np.cos(config.Az_rad) * north_hat))
            thrustF = T_mag_total * thrust_dir
        else:
            # Phase 3: True gravity turn — thrust along velocity vector
            thrust_dir = v_vec / V
            thrustF = T_mag_total * thrust_dir
    else:
        thrustF = np.zeros(3)

    # ---- Total forces & acceleration ----
    Forces = gravityF + aeroF + thrustF

    if mass > 0:
        vdot = Forces / mass
    else:
        vdot = np.zeros(3)
        mdot = 0.0

    return np.array([xdot, ydot, zdot, vdot[0], vdot[1], vdot[2], mdot])


############################################################################
# ======================= EVENT FUNCTIONS =======================
############################################################################

def make_ground_hit_event():
    """Creates a terminal event that fires when the rocket hits the ground."""
    def ground_hit(t, state, config, stage_index, boosters_active):
        r = np.sqrt(state[0]**2 + state[1]**2 + state[2]**2)
        return r - Rplanet
    ground_hit.terminal = True
    ground_hit.direction = -1
    return ground_hit


def make_booster_separation_event(t_burn_boosters):
    """
    Creates a terminal event that fires at booster burnout.
    Terminal so we can stop, subtract booster dry mass, and restart.
    """
    def booster_sep(t, state, config, stage_index, boosters_active):
        return t - t_burn_boosters
    booster_sep.terminal = True
    booster_sep.direction = 1
    return booster_sep


def make_engine_cutoff_event(t_burn):
    """Creates a non-terminal event that records engine cutoff."""
    def engine_cutoff(t, state, config, stage_index, boosters_active):
        return t - t_burn
    engine_cutoff.terminal = False
    engine_cutoff.direction = 1
    return engine_cutoff


############################################################################
# ======================= INTEGRATION (EVENT-DRIVEN) =======================
############################################################################

def simulate_stage(stateinitial, config: FlightConfig, stage_index: int):
    """
    Simulates a single stage flight from t=0 using solve_ivp with events.

    For stage 0 with boosters: fires a terminal booster_separation event,
    subtracts booster dry mass, then continues integration.

    Returns
    -------
    dict with keys:
        't'           : time array
        'state'       : state array (N x 7)
        't_burn_end'  : time array for burn phase only
        'state_burn'  : state array for burn phase only
    """
    stage = config.stages[stage_index]
    has_boosters = (stage_index == 0 and config.boosters is not None and len(config.boosters) > 0)
    max_time = 3000.0

    # ---- Phase 1: Integration (possibly up to booster separation) ----
    events_phase1 = [make_ground_hit_event()]
    if has_boosters:
        events_phase1.append(make_booster_separation_event(config.t_burn_boosters))
    events_phase1.append(make_engine_cutoff_event(stage.t_burn))

    sol1 = solve_ivp(
        Derivatives,
        [0, max_time],
        stateinitial,
        args=(config, stage_index, has_boosters),
        events=events_phase1,
        method='RK45',
        rtol=1e-8,
        atol=1e-10,
        max_step=0.5,
        dense_output=True
    )

    # Check if booster separation event fired (event index 1 when boosters present)
    if has_boosters and len(sol1.t_events) > 1 and sol1.t_events[1].size > 0:
        # Booster separation occurred
        t_sep = sol1.t_events[1][0]
        state_sep = sol1.y_events[1][0].copy()

        # Subtract total booster dry mass
        booster_dry_mass = sum(b.m_construction for b in config.boosters)
        state_sep[6] -= booster_dry_mass

        # ---- Phase 2: Continue without boosters ----
        events_phase2 = [make_ground_hit_event(), make_engine_cutoff_event(stage.t_burn)]

        sol2 = solve_ivp(
            Derivatives,
            [t_sep, max_time],
            state_sep,
            args=(config, stage_index, False),  # boosters_active = False
            events=events_phase2,
            method='RK45',
            rtol=1e-8,
            atol=1e-10,
            max_step=0.5,
            dense_output=True
        )

        # Merge the two solutions
        t_full = np.concatenate([sol1.t, sol2.t])
        state_full = np.concatenate([sol1.y.T, sol2.y.T], axis=0)

        # Extract burn phase: everything up to t_burn
        burn_mask = t_full <= stage.t_burn
        t_burn_arr = t_full[burn_mask]
        state_burn_arr = state_full[burn_mask]

    else:
        # No booster separation (stage > 0, or no boosters)
        t_full = sol1.t
        state_full = sol1.y.T

        burn_mask = t_full <= stage.t_burn
        t_burn_arr = t_full[burn_mask]
        state_burn_arr = state_full[burn_mask]

    return {
        't': t_full,
        'state': state_full,
        't_burn': t_burn_arr,
        'state_burn': state_burn_arr,
    }


def simulate_booster_post_separation(state_at_separation, config: FlightConfig, booster_index: int):
    """
    Simulates a single booster's ballistic flight after separation.

    Parameters
    ----------
    state_at_separation : np.ndarray
        State at separation [x, y, z, vx, vy, vz, mass_of_single_booster]
    config : FlightConfig
        Flight configuration
    booster_index : int
        Index of the booster

    Returns
    -------
    dict with keys 't', 'state'
    """
    booster = config.boosters[booster_index]
    max_time = 3000.0

    # Create a minimal config for the booster's ballistic flight
    booster_stage = RocketStage(
        T_mag=0.0,          # no thrust after separation
        m_flow=0.0,
        t_burn=0.0,         # already burned out
        m_construction=booster.m_construction,
        Area_pf=booster.Area_pf,
        Area_bf=booster.Area_bf,
    )

    booster_config = FlightConfig(
        stages=[booster_stage],
        boosters=None,
        t_vertical=0.0,     # irrelevant — no thrust
        kick_angle_rad=0.0,
        Az_rad=config.Az_rad,
        Cd_crossflow=config.Cd_crossflow,
        t_burn_boosters=0.0,
    )

    events = [make_ground_hit_event()]

    sol = solve_ivp(
        Derivatives,
        [0, max_time],
        state_at_separation,
        args=(booster_config, 0, False),
        events=events,
        method='RK45',
        rtol=1e-8,
        atol=1e-10,
        max_step=0.5,
    )

    return {
        't': sol.t,
        'state': sol.y.T,
    }


def simulate_two_burn_stage(state, config: FlightConfig, stage_index: int,
                             fuel_reserve_fraction=0.20):
    """
    Two-burn approach for the last stage:
    1. First burn: use (1 - fuel_reserve_fraction) of the fuel -> transfer orbit
    2. Coast to apogee (ballistic, event: dot(r,v) = 0, direction -1)
    3. Second burn at apogee: use remaining fuel to circularize (thrust along velocity)
    
    Parameters
    ----------
    state : np.ndarray
        Initial state [x, y, z, vx, vy, vz, mass]
    config : FlightConfig
        Flight configuration  
    stage_index : int
        Index of this stage
    fuel_reserve_fraction : float
        Fraction of fuel to reserve for the second burn (0.20 = 20%)
    
    Returns
    -------
    dict with keys:
        't' : full time array (all 3 phases)
        'state' : full state array (all 3 phases)
        'burnout_state_1' : state after first burn
        'apogee_state' : state at apogee (start of second burn)
        'burnout_state_2' : state after second burn (final orbit)
        't_coast' : coast duration [s]
        'phase_times' : dict with 't_burn1', 't_coast', 't_burn2'
    """
    stage = config.stages[stage_index]
    
    # --- Phase 1: First burn (partial fuel) ---
    m_prop_total = stage.m_flow * stage.t_burn  # total propellant
    m_prop_burn1 = m_prop_total * (1.0 - fuel_reserve_fraction)
    t_burn1 = m_prop_burn1 / stage.m_flow
    
    # Create a modified stage with shorter burn time for phase 1
    stage_burn1 = RocketStage(
        T_mag=stage.T_mag,
        m_flow=stage.m_flow,
        t_burn=t_burn1,
        m_construction=stage.m_construction,
        Area_pf=stage.Area_pf,
        Area_bf=stage.Area_bf,
    )
    
    # Replace the stage in config temporarily
    original_stages = config.stages
    config_burn1_stages = list(config.stages)
    config_burn1_stages[stage_index] = stage_burn1
    config.stages = config_burn1_stages
    
    # Simulate first burn
    result1 = simulate_stage(state, config, stage_index)
    
    # Restore original config
    config.stages = original_stages
    
    # Get state after first burn
    burnout_state_1 = get_burnout_state(result1)
    
    # --- Phase 2: Coast to apogee (no thrust, wait for rdot = 0) ---
    # State for coast: same position/velocity, but mass stays the same (no fuel burn)
    coast_state = burnout_state_1.copy()
    
    def apogee_event(t, state):
        """dot(r, v) = 0 at apogee (radial velocity = 0)"""
        r_vec = state[0:3]
        v_vec = state[3:6]
        return np.dot(r_vec, v_vec)
    apogee_event.terminal = True
    apogee_event.direction = -1  # fire when rdot goes from + to - (passing through apogee)
    
    def ground_hit_coast(t, state):
        r = np.sqrt(state[0]**2 + state[1]**2 + state[2]**2)
        return r - Rplanet
    ground_hit_coast.terminal = True
    ground_hit_coast.direction = -1
    
    max_coast_time = 50000.0  # max coast time
    
    sol_coast = solve_ivp(
        _gravity_only_derivatives,
        [0, max_coast_time],
        coast_state[:6],  # _gravity_only_derivatives uses 6-element state
        events=[apogee_event, ground_hit_coast],
        method='RK45',
        rtol=1e-10,
        atol=1e-12,
        max_step=1.0,
    )
    
    t_coast = sol_coast.t[-1]
    apogee_state_6 = sol_coast.y[:, -1].copy()
    
    # Reconstruct 7-element state with mass
    apogee_state = np.append(apogee_state_6, coast_state[6])
    
    # --- Phase 3: Second burn at apogee (circularization) ---
    m_prop_burn2 = m_prop_total * fuel_reserve_fraction
    t_burn2 = m_prop_burn2 / stage.m_flow
    
    # Create a stage for the second burn
    stage_burn2 = RocketStage(
        T_mag=stage.T_mag,
        m_flow=stage.m_flow,
        t_burn=t_burn2,
        m_construction=stage.m_construction,
        Area_pf=stage.Area_pf,
        Area_bf=stage.Area_bf,
    )
    
    # Config for burn 2: stage_index doesn't matter much, but we use a high index
    # so that the thrust direction goes straight to Phase 3 (velocity-following)
    config_burn2_stages = list(config.stages)
    config_burn2_stages[stage_index] = stage_burn2
    config.stages = config_burn2_stages
    
    # Set stage_index > 0 so that Derivatives skips vertical/kick phases
    # and goes straight to velocity-following (Phase 3)
    result3 = simulate_stage(apogee_state, config, stage_index=max(stage_index, 1))
    
    # Restore original config
    config.stages = original_stages
    
    # Get state after second burn
    burnout_state_2 = get_burnout_state(result3)
    
    # --- Merge all phases ---
    # Phase 1 times
    t_phase1 = result1['t']
    state_phase1 = result1['state']
    
    # Phase 2: coast times (offset by end of phase 1)
    t_offset_coast = t_phase1[-1]
    t_phase2 = sol_coast.t + t_offset_coast
    # Add mass column to coast states
    coast_states_with_mass = np.column_stack([
        sol_coast.y.T,
        np.full(len(sol_coast.t), coast_state[6])
    ])
    
    # Phase 3: burn 2 times (offset by end of phase 2)
    t_offset_burn2 = t_phase2[-1]
    t_phase3 = result3['t'] + t_offset_burn2
    state_phase3 = result3['state']
    
    t_full = np.concatenate([t_phase1, t_phase2, t_phase3])
    state_full = np.concatenate([state_phase1, coast_states_with_mass, state_phase3], axis=0)
    
    return {
        't': t_full,
        'state': state_full,
        'burnout_state_1': burnout_state_1,
        'apogee_state': apogee_state,
        'burnout_state_2': burnout_state_2,
        't_coast': t_coast,
        'phase_times': {
            't_burn1': t_burn1,
            't_coast': t_coast,
            't_burn2': t_burn2,
        },
    }


############################################################################
# ======================= HELPER FUNCTIONS =======================
############################################################################

def extract_results(result_dict):
    """Extracts position, velocity, and mass arrays from a simulation result."""
    state = result_dict['state']
    return {
        'x': state[:, 0],
        'y': state[:, 1],
        'z': state[:, 2],
        'vx': state[:, 3],
        'vy': state[:, 4],
        'vz': state[:, 5],
        'mass': state[:, 6],
        'velmag': np.sqrt(state[:, 3]**2 + state[:, 4]**2 + state[:, 5]**2),
        'altitude': np.sqrt(state[:, 0]**2 + state[:, 1]**2 + state[:, 2]**2) - Rplanet,
    }


def get_burnout_state(result_dict):
    """Returns the state vector at the end of the burn phase."""
    state_burn = result_dict['state_burn']
    if len(state_burn) > 0:
        return state_burn[-1].copy()
    else:
        return result_dict['state'][0].copy()


############################################################################
# ======================= TIME & COORDINATE CONVERSIONS =======================
############################################################################

def compute_corrected_azimuth(lat_rad, inc_rad, h_target, mu=G * Mplanet):
    """
    Compute the corrected launch azimuth accounting for Earth's rotation.
    
    Formulas:
        Vorb  = sqrt(mu / (RE + h))
        Vpad  = (2*pi*RE / t_day) * cos(La)
        Az    = arcsin(cos(i) / cos(La))       # geometric azimuth
        VL_N  = Vorb * cos(Az)
        VL_E  = Vorb * sin(Az) - Vpad
        VL    = sqrt(VL_N^2 + VL_E^2)
        AzL   = arctan2(VL_E, VL_N)            # corrected azimuth
    
    Parameters
    ----------
    lat_rad : float
        Launch latitude [rad]
    inc_rad : float
        Target inclination [rad]
    h_target : float
        Target orbit altitude above surface [m]
    mu : float
        Gravitational parameter [m^3/s^2]
    
    Returns
    -------
    dict with keys:
        'Az_geometric' : geometric azimuth [rad]
        'Az_corrected' : corrected azimuth [rad]
        'Vorb'  : orbital velocity at target altitude [m/s]
        'Vpad'  : launch pad velocity from Earth rotation [m/s]
        'VL'    : required launch velocity magnitude [m/s]
        'dV_rot' : velocity gain from Earth rotation [m/s]
        'dV_rot_corr' : corrected velocity gain [m/s]
    """
    t_day = 86400.0  # sidereal day approximation [s]
    
    r_target = Rplanet + h_target
    Vorb = np.sqrt(mu / r_target)
    
    # Pad velocity from Earth's rotation
    Vpad = (2.0 * np.pi * Rplanet / t_day) * np.cos(lat_rad)
    
    # Geometric azimuth
    cos_i_over_cos_lat = np.cos(inc_rad) / np.cos(lat_rad)
    # Clamp to [-1, 1] for safety
    cos_i_over_cos_lat = np.clip(cos_i_over_cos_lat, -1.0, 1.0)
    Az_geo = np.arcsin(cos_i_over_cos_lat)
    
    # Velocity components in launch frame
    VL_North = Vorb * np.cos(Az_geo)
    VL_East = Vorb * np.sin(Az_geo) - Vpad
    VL = np.sqrt(VL_North**2 + VL_East**2)
    
    # Corrected azimuth
    Az_corr = np.arctan2(VL_East, VL_North)
    
    # Velocity gain from Earth rotation
    dV_rot = Vpad * np.sin(Az_geo)
    dV_rot_corr = Vorb - VL
    
    return {
        'Az_geometric': Az_geo,
        'Az_corrected': Az_corr,
        'Vorb': Vorb,
        'Vpad': Vpad,
        'VL': VL,
        'dV_rot': dV_rot,
        'dV_rot_corr': dV_rot_corr,
    }


def utc_to_julian_date(year, month, day, hour, minute, second):
    """
    Convert UTC time to Julian Date.
    
    Parameters
    ----------
    year, month, day : int
    hour, minute, second : int/float
    
    Returns
    -------
    float : Julian Date
    """
    # Fractional day
    Df = day + hour / 24.0 + minute / 1440.0 + second / 86400.0
    
    # Leap year check
    if month <= 2:
        Y_prime = year - 1
        M_prime = month + 12
    else:
        Y_prime = year
        M_prime = month
    
    # Gregorian correction
    A = int(Y_prime / 100)
    B = 2 - A + int(A / 4)
    
    # Julian Date
    JD = int(365.25 * (Y_prime + 4716)) + int(30.6001 * (M_prime + 1)) + Df + B - 1524.5
    
    return JD


def compute_gmst(jd_utc):
    """
    Compute Greenwich Mean Sidereal Time (GMST) from Julian Date.
    
    Parameters
    ----------
    jd_utc : float
        Julian Date in UTC
    
    Returns
    -------
    float : GMST in radians
    """
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
    """
    Convert state vector from ECI Greenwich-based frame to J2000 inertial frame.
    
    The ECI Greenwich frame has its x-axis along the Greenwich meridian at launch.
    J2000 has its x-axis pointing to the vernal equinox at J2000 epoch.
    The rotation angle is GMST (angle from vernal equinox to Greenwich, eastward).
    
    Parameters
    ----------
    state : np.ndarray
        State vector [x, y, z, vx, vy, vz] or [x, y, z, vx, vy, vz, mass]
    gmst_rad : float
        GMST at the time of the state, in radians
    
    Returns
    -------
    np.ndarray : state in J2000 frame
    """
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


def keplerian_to_cartesian(a, e, i_rad, RAAN_rad, omega_rad, nu_rad, mu=G * Mplanet):
    """
    Convert Keplerian orbital elements to Cartesian state vector.
    
    Parameters
    ----------
    a : float - semi-major axis [m]
    e : float - eccentricity
    i_rad : float - inclination [rad]
    RAAN_rad : float - right ascension of ascending node [rad]
    omega_rad : float - argument of periapsis [rad]
    nu_rad : float - true anomaly [rad]
    mu : float - gravitational parameter [m^3/s^2]
    
    Returns
    -------
    np.ndarray : [x, y, z, vx, vy, vz] in the reference frame
    """
    # Semi-latus rectum
    p = a * (1 - e**2)
    r_mag = p / (1 + e * np.cos(nu_rad))
    
    # Position and velocity in perifocal frame (PQW)
    r_pqw = np.array([
        r_mag * np.cos(nu_rad),
        r_mag * np.sin(nu_rad),
        0.0
    ])
    
    v_pqw = np.sqrt(mu / p) * np.array([
        -np.sin(nu_rad),
        e + np.cos(nu_rad),
        0.0
    ])
    
    # Rotation matrix from perifocal to inertial (313 rotation: -RAAN, -i, -omega)
    cos_O = np.cos(RAAN_rad)
    sin_O = np.sin(RAAN_rad)
    cos_i = np.cos(i_rad)
    sin_i = np.sin(i_rad)
    cos_w = np.cos(omega_rad)
    sin_w = np.sin(omega_rad)
    
    R = np.array([
        [cos_O * cos_w - sin_O * sin_w * cos_i,
         -cos_O * sin_w - sin_O * cos_w * cos_i,
         sin_O * sin_i],
        [sin_O * cos_w + cos_O * sin_w * cos_i,
         -sin_O * sin_w + cos_O * cos_w * cos_i,
         -cos_O * sin_i],
        [sin_w * sin_i,
         cos_w * sin_i,
         cos_i]
    ])
    
    r_inertial = R @ r_pqw
    v_inertial = R @ v_pqw
    
    return np.concatenate([r_inertial, v_inertial])


def cartesian_to_keplerian(state, mu=G * Mplanet):
    """
    Convert Cartesian state vector to Keplerian orbital elements.
    
    Parameters
    ----------
    state : np.ndarray
        [x, y, z, vx, vy, vz] (mass component ignored if present)
    mu : float
        Gravitational parameter [m^3/s^2]
    
    Returns
    -------
    dict with keys:
        'a' : semi-major axis [m]
        'e' : eccentricity
        'i' : inclination [rad]
        'RAAN' : right ascension of ascending node [rad]
        'omega' : argument of periapsis [rad]
        'nu' : true anomaly [rad]
    """
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
    e_vec = (np.cross(v_vec, h_vec) / mu) - (r_vec / r)
    e = np.linalg.norm(e_vec)
    
    # Specific orbital energy
    epsilon = v**2 / 2.0 - mu / r
    
    # Semi-major axis
    if abs(1 - e) > 1e-10:  # not parabolic
        a = -mu / (2 * epsilon)
    else:
        a = np.inf
    
    # Inclination
    i = np.arccos(np.clip(h_vec[2] / h, -1, 1))
    
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
# ======================= BACKWARD PROPAGATION =======================
############################################################################

def _gravity_only_derivatives(t, state):
    """
    Equations of motion with gravity only (no thrust, no drag).
    Used for backward propagation from target orbit.
    State: [x, y, z, vx, vy, vz]
    """
    x, y, z = state[0], state[1], state[2]
    velx, vely, velz = state[3], state[4], state[5]
    
    g = gravity(x, y, z)
    
    return np.array([velx, vely, velz, g[0], g[1], g[2]])


def backward_propagation(a_target, i_target_rad, RAAN_rad, min_altitude,
                          perigee_altitude=200000.0, mu=G * Mplanet):
    """
    Propagate backward from the apogee of a transfer orbit to find the
    burnout state at a reference (perigee) altitude.
    
    The transfer orbit has:
    - Apogee at the target orbit altitude (r_a = a_target for circular target)
    - Perigee at perigee_altitude (estimated burnout altitude)
    
    At apogee the velocity is purely tangential and less than circular velocity,
    so backward integration under gravity naturally descends to the perigee.
    
    Parameters
    ----------
    a_target : float
        Target semi-major axis [m] (= orbit radius for circular orbit)
    i_target_rad : float
        Target inclination [rad]
    RAAN_rad : float
        Right ascension of ascending node [rad]
    min_altitude : float
        Stop when altitude reaches this value [m]
    perigee_altitude : float
        Estimated burnout (perigee) altitude [m] — used to define the transfer orbit
    mu : float
        Gravitational parameter [m^3/s^2]
    
    Returns
    -------
    dict with keys:
        'target_orbit_state' : state on the target circular orbit (for reference)
        'transfer_apogee_state' : state at apogee of transfer orbit (start of backward prop)
        'burnout_target_state' : state at the reference altitude (backward endpoint)
        't_coast' : coast time from burnout to orbit (positive)
        'trajectory_t' : full backward trajectory time array
        'trajectory_state' : full backward trajectory state array
        'a_transfer' : transfer orbit SMA [m]
    """
    r_apogee = a_target  # for circular target orbit, apogee = SMA
    r_perigee = Rplanet + perigee_altitude
    
    # Transfer orbit semi-major axis
    a_transfer = (r_apogee + r_perigee) / 2.0
    
    # Transfer orbit eccentricity
    e_transfer = (r_apogee - r_perigee) / (r_apogee + r_perigee)
    
    # State on the target circular orbit (for reference only)
    target_circ_state = keplerian_to_cartesian(
        a_target, 0.0, i_target_rad, RAAN_rad, 0.0, 0.0, mu
    )
    
    # State at apogee of transfer orbit (nu = pi at apogee)
    # Using omega=0, nu=pi means spacecraft is at apogee
    transfer_apogee_state = keplerian_to_cartesian(
        a_transfer, e_transfer, i_target_rad, RAAN_rad,
        0.0, np.pi, mu  # omega=0, nu=pi (apogee)
    )
    
    # Propagate backward (negative time) from apogee
    # Half the transfer orbit period is the coast time
    T_transfer = 2 * np.pi * np.sqrt(a_transfer**3 / mu)
    max_backward_time = T_transfer  # more than enough
    
    def altitude_event(t, state):
        r = np.sqrt(state[0]**2 + state[1]**2 + state[2]**2)
        return (r - Rplanet) - min_altitude
    altitude_event.terminal = True
    altitude_event.direction = 0  # fire in either direction (needed for negative-time integration)
    
    sol = solve_ivp(
        _gravity_only_derivatives,
        [0, -max_backward_time],  # negative time
        transfer_apogee_state,
        events=altitude_event,
        method='RK45',
        rtol=1e-10,
        atol=1e-12,
        max_step=1.0,
    )
    
    # The burnout target state is the last point (either event or end)
    burnout_state = sol.y[:, -1].copy()
    t_coast = abs(sol.t[-1])  # positive coast time
    
    return {
        'target_orbit_state': target_circ_state,
        'transfer_apogee_state': transfer_apogee_state,
        'burnout_target_state': burnout_state,
        't_coast': t_coast,
        'trajectory_t': sol.t,
        'trajectory_state': sol.y.T,
        'a_transfer': a_transfer,
    }


############################################################################
# ======================= KICK ANGLE OPTIMIZATION =======================
############################################################################

def find_kick_angle(run_sim_func, data, a_target_m, e_target=0.0, mu=G * Mplanet):
    """
    Find the optimal kick angle to achieve the target orbit (SMA + eccentricity).
    
    Uses scipy.optimize.minimize_scalar to minimize a combined objective:
      error = (a_achieved - a_target)^2 / a_target^2 + w_e * (e_achieved - e_target)^2
    
    The final orbit is taken from orbit_info['keplerian'] which reflects
    the post-circularization state (after two-burn approach).
    
    Parameters
    ----------
    run_sim_func : callable
        The run_simulation function (returns 3-tuple)
    data : dict
        Simulation input data
    a_target_m : float
        Target semi-major axis [m]
    e_target : float
        Target eccentricity (0 for circular)
    mu : float
        Gravitational parameter [m^3/s^2]
    
    Returns
    -------
    dict with keys:
        'kick_angle_deg' : optimal kick angle [degrees]
        'a_achieved' : achieved SMA [m]
        'e_achieved' : achieved eccentricity
        'i_achieved_deg' : achieved inclination [deg]
        'optimization_result' : full scipy result
    """
    w_e = 100.0  # weight for eccentricity in objective
    
    def objective(kick_angle_deg):
        data_trial = data.copy()
        data_trial['kick_angle'] = kick_angle_deg
        
        try:
            _, _, orbit_info = run_sim_func(data_trial)
            kep = orbit_info['keplerian']
            
            a = kep['a']
            e = kep['e']
            i_deg = np.degrees(kep['i'])
            
            # Normalized SMA error + eccentricity penalty
            sma_err = ((a - a_target_m) / a_target_m) ** 2
            ecc_err = (e - e_target) ** 2
            
            error = sma_err + w_e * ecc_err
            print(f"  kick={kick_angle_deg:.4f} deg -> SMA={a/1000:.0f} km, e={e:.4f}, i={i_deg:.1f} deg, err={error:.2e}")
            return error
        except Exception as ex:
            print(f"  Kick angle {kick_angle_deg:.3f} deg failed: {ex}")
            return 1e20
    
    print("Optimizing kick angle...")
    result = minimize_scalar(
        objective,
        bounds=(0.01, 30.0),
        method='bounded',
        options={'xatol': 1e-3, 'maxiter': 50}
    )
    
    kick_opt = result.x
    
    # Run final simulation with optimal kick angle
    data_final = data.copy()
    data_final['kick_angle'] = kick_opt
    _, _, orbit_info = run_sim_func(data_final)
    kep = orbit_info['keplerian']
    
    a_achieved = kep['a']
    e_achieved = kep['e']
    i_achieved = np.degrees(kep['i'])
    
    print(f"\nOptimal kick angle: {kick_opt:.4f} deg")
    print(f"Target SMA:   {a_target_m/1000:.1f} km  |  Achieved: {a_achieved/1000:.1f} km")
    print(f"Target e:     {e_target:.4f}      |  Achieved: {e_achieved:.4f}")
    print(f"Inclination:  {i_achieved:.2f} deg")
    
    return {
        'kick_angle_deg': kick_opt,
        'a_achieved': a_achieved,
        'e_achieved': e_achieved,
        'i_achieved_deg': i_achieved,
        'optimization_result': result,
    }


############################################################################
############################################################################
# ======================= COORDINATE CONVERSION =======================
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
# ======================= ROCKET PARAMETER CALCULATIONS =======================
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
