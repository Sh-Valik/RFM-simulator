from functions import derivatives, mass_of_stages_eps_mode
import numpy as np
from scipy.integrate import odeint


G = 6.6742 * 10**-11  # gravitational constant [N.m^2/kg^2]
g0 = 9.80665  # standard gravitational acceleration [m/s^2]
Rplanet = 6371000  # mean radius of the Earth [m]
Mplanet = 5.97219 * 10**24  # mass of the Earth [kg]


def calculate_flight(data):
    m_payload = data["payload_mass"]
    theta_angle_deg = data["theta_angle"]
    theta_angle = np.radians(theta_angle_deg)
    stages_count = data["stages_count"]
    has_boosters = data["has_boosters"]
    booster_count = data["booster_count"]
    input_mode = data["input_mode"]
    if input_mode == "Start mass & Propellant":
        stages_data = data["stages_data_mass"]
    else:
        stages_data = mass_of_stages_eps_mode(data["stages_data_eps"]) # Placeholder for future implementation
    if has_boosters:
        boosters_data = data["boosters_data"]
        if input_mode == "Start mass & Propellant":
            boosters_data = data["boosters_data_mass"]
        else:
            boosters_data = mass_of_stages_eps_mode(data["boosters_data_eps"]) # Placeholder for future implementation
    
    launch_lat = data["launch_lat"]
    launch_lon = data["launch_lon"]
    launch_alt = data["launch_alt"]
    launch_date = data["launch_date"]
    launch_time = data["launch_time"]

    t_o_a = data["orbit_a"]
    t_o_e = data["orbit_e"]
    t_o_i = data["orbit_i"]
    t_o_Ω = data["orbit_Ω"]
    t_o_ω = data["orbit_ω"]
    t_o_ν = data["orbit_ν"]


    
