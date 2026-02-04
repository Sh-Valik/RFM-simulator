from Algorithm.functions import Derivatives, parameters_of_stages, parameters_of_boosters
import numpy as np
from scipy.integrate import odeint
import json
import os
from data_manager import load_data

G = 6.6742 * 10**-11  # gravitational constant [N.m^2/kg^2]
g0 = 9.80665  # standard gravitational acceleration [m/s^2]
Rplanet = 6371000  # mean radius of the Earth [m]
Mplanet = 5.97219 * 10**24  # mass of the Earth [kg]

############################################################################
# Rocket and mission parameters
def run_simulation(data):
    m_payload_without_boosters = data["payload_mass_without_booster"]
    m_payload_with_boosters = data["payload_mass_with_booster"]

    # theta_angle_deg = data["theta_angle"]
    # theta_angle = np.radians(theta_angle_deg)

    stages_count = data["stages_count"]
    has_boosters = data["has_boosters"]
    if has_boosters:
        booster_count = data["booster_count"]
        t_burn_ratio = data["t_burn_ratio"]

    payload_mass_ratio_total = data["payload_mass_ratio_total"]
    t_vertical_flight = data["t_vertical_flight"]

    input_mode = data["input_mode"]
    rocket_type = data["rocket_type"]

    if input_mode == "Start mass & Propellant":
        stage_data_list = data["stages_data_mass"]
        Ve_stages, mass_flow_stages, m0_stages, m_prop_stages, Vf_id_stages, Lambda_stages = parameters_of_stages(input_mode, stage_data_list, m_payload_without_boosters, payload_mass_ratio_total, rocket_type, stages_count) # Need to be changed due the logic of the algorithm will change slightly
        if has_boosters:
            booster_data_list = data["boosters_data_mass"]
            # palceholder

    else:
        stage_data_list = data["stages_data_eps"]
        Ve_stages, mass_flow_stages, diameter_stages, height_stages, m0_stages, m_prop_stages, Vf_id_stages, Lambda_stages = parameters_of_stages(input_mode, stage_data_list, m_payload_without_boosters, payload_mass_ratio_total, rocket_type, stages_count) # Need to be changed due the logic of the algorithm will change slightly
        if has_boosters:
            booster_data_list = data["boosters_data_eps"]
            Ve_boosters, mass_flow_boosters, diameter_boosters, height_boosters, m0_each_boosters, m_construction_each_boosters, m_prop_each_boosters, m0 = parameters_of_boosters(input_mode, booster_data_list, m_payload_without_boosters, m_payload_with_boosters, Vf_id_stages, m0_stages, m_prop_stages, stages_count, booster_count, Ve_stages, mass_flow_stages, t_burn_ratio)


    launch_lat = data["launch_lat"]
    launch_lon = data["launch_lon"]
    launch_alt = data["launch_alt"]
    launch_date = data["launch_date"]
    launch_time = data["launch_time"]

    t_o_a = data["orbit_a"]
    t_o_e = data["orbit_e"]
    t_o_i = data["orbit_i"]
    #######################################################

    # Cross-sectional areas
    stages_area_pf = [np.pi * (diameter_stages[i]**2) / 4 for i in range(stages_count)] # pf - propelled flight
    stages_area_bf = [diameter_stages[i] * height_stages[i] for i in range(stages_count)] # bf - ballistic flight
    
    if has_boosters:
        boosters_area_pf = [np.pi * (diameter_boosters[i]**2) / 4 for i in range(booster_count)]
        boosters_area_bf = [diameter_boosters[i] * height_boosters[i] for i in range(booster_count)]

    Cd_of_crosflow_cylinder = 1.25 # for stages and boosters in ballistic flight
    #######################################################

    # Thrust magnitude for each stage
    T_mag_stages = [mass_flow_stages[i] * Ve_stages[i] for i in range(stages_count)]
    if has_boosters:
        T_mag_boosters = [mass_flow_boosters[i] * Ve_boosters[i] for i in range(booster_count)]
    #######################################################
    
    