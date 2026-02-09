from Algorithm.functions import parameters_of_stages, parameters_of_boosters, geodetic_to_cartesian_WGS84, integration
import numpy as np
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
    theta = 45 # need to be changed 
    theta = np.radians(theta) 
    #######################################################

    # Cross-sectional areas
    stages_area_pf = [np.pi * (diameter_stages[i]**2) / 4 for i in range(stages_count)] # pf - propelled flight
    stages_area_bf = [diameter_stages[i] * height_stages[i] for i in range(stages_count)] # bf - ballistic flight
    
    if has_boosters:
        boosters_area_pf = [np.pi * (diameter_boosters[i]**2) / 4 for i in range(booster_count)]
        boosters_area_bf = [diameter_boosters[i] * height_boosters[i] for i in range(booster_count)]

    Cd_of_crosflow_cylinder = 1.25 # for stages and boosters in ballistic flight
    #######################################################
    
    # Burn time
    t_burn_stages = [m_prop_stages[i] / mass_flow_stages[i] for i in range(stages_count)]
    m_construction_stages = [m0_stages[i] - m_prop_stages[i] for i in range(stages_count)]
    if has_boosters:
        t_burn_boosters = t_burn_stages[0] * t_burn_ratio
    #######################################################    

    # Thrust magnitude for each stage
    T_mag_stages = [mass_flow_stages[i] * Ve_stages[i] for i in range(stages_count)]
    if has_boosters:
        T_mag_boosters = [mass_flow_boosters[i] * Ve_boosters[i] for i in range(booster_count)]
    #######################################################
    
    x0, y0, z0 = geodetic_to_cartesian_WGS84(launch_lat, launch_lon, launch_alt)
    

    # velx0, vely0, velz0 = launch_velocity(Rplanet, launch_lat, SMA=t_o_a, inclination=t_o_i)
    Vpad = ((2 * np.pi * Rplanet) / (24*60*60)) * np.cos(np.deg2rad(launch_lat))
    velx0 = Vpad
    vely0 = 0.0
    velz0 = 0.0

    Azimuth = np.asin((np.cos(np.deg2rad(t_o_i))) / (np.cos(np.deg2rad(launch_lat))))
    

    stages_info = [t_burn_stages, T_mag_stages, mass_flow_stages, m_prop_stages, m_construction_stages]
    boosters_info = [t_burn_boosters, T_mag_boosters, mass_flow_boosters, booster_count]
    #######################################################
    # Define all necessary arrays for each stage
    stateinitial_stages = [None] * stages_count
    q_stages = [None] * stages_count

    tout_stages = [None] * stages_count
    stateout_stages = [None] * stages_count
    tout_b_stages = [None] * stages_count
    stateout_b_stages = [None] * stages_count

    xout_stages = [None] * stages_count
    yout_stages = [None] * stages_count
    zout_stages = [None] * stages_count

    xout_b_stages = [None] * stages_count
    yout_b_stages = [None] * stages_count
    zout_b_stages = [None] * stages_count

    velxout_stages = [None] * stages_count
    velyout_stages = [None] * stages_count
    velzout_stages = [None] * stages_count
    velmag_stages = [None] * stages_count

    velxout_b_stages = [None] * stages_count
    velyout_b_stages = [None] * stages_count
    velzout_b_stages = [None] * stages_count

    massout_stages = [None] * stages_count
    
    # Define all necessary arrays for each booster
    if has_boosters:
        stateinitial_boosters = [None] * booster_count
        q_boosters = [None] * booster_count

        tout_boosters = [None] * booster_count
        stateout_boosters = [None] * booster_count
        tout_b_boosters = [None] * booster_count
        stateout_b_boosters = [None] * booster_count

        xout_boosters = [None] * booster_count
        yout_boosters = [None] * booster_count
        zout_boosters = [None] * booster_count

        xout_b_boosters = [None] * booster_count
        yout_b_boosters = [None] * booster_count
        zout_b_boosters = [None] * booster_count

        velxout_boosters = [None] * booster_count
        velyout_boosters = [None] * booster_count
        velzout_boosters = [None] * booster_count
        velmag_boosters = [None] * booster_count

        velxout_b_boosters = [None] * booster_count
        velyout_b_boosters = [None] * booster_count
        velzout_b_boosters = [None] * booster_count

        massout_boosters = [None] * booster_count

#######################################################
    for i in range(stages_count):
        if has_boosters:
            if i == 0:
                stateinitial_stages[i] = np.array([x0, y0, z0, velx0, vely0, velz0, m0[i]])
                
                tout_stages[i], stateout_stages[i], tout_b_stages[i], stateout_b_stages[i] = integration(stateinitial_stages[i], stages_info, boosters_info, m_construction_stages[i], stages_area_pf[i], stages_area_bf[i], Cd_of_crosflow_cylinder, t_vertical_flight, theta, Azimuth, t_burn_stages[i], T_mag_stages[i], mass_flow_stages[i], has_boosters)

                xout_stages[i] = stateout_stages[i][:, 0]
                yout_stages[i] = stateout_stages[i][:, 1]
                zout_stages[i] = stateout_stages[i][:, 2]

                velxout_stages[i] = stateout_stages[i][:, 3]
                velyout_stages[i] = stateout_stages[i][:, 4]
                velzout_stages[i] = stateout_stages[i][:, 5]

                velmag_stages[i] = np.sqrt(velxout_stages[i]**2 + velyout_stages[i]**2 + velzout_stages[i]**2)

                xout_b_stages[i] = stateout_b_stages[i][:, 0]
                yout_b_stages[i] = stateout_b_stages[i][:, 1]
                zout_b_stages[i] = stateout_b_stages[i][:, 2]

                velxout_b_stages[i] = stateout_b_stages[i][:, 3]
                velyout_b_stages[i] = stateout_b_stages[i][:, 4]
                velzout_b_stages[i] = stateout_b_stages[i][:, 5]

                massout_stages[i] = stateout_b_stages[i][:, 6]

            else:
                stateinitial_stages[i] = np.array([xout_b_stages[i-1][-1], yout_b_stages[i-1][-1], zout_b_stages[i-1][-1], velxout_b_stages[i-1][-1], velyout_b_stages[i-1][-1], velzout_b_stages[i-1][-1], m0[i]])

                tout_stages[i], stateout_stages[i], tout_b_stages[i], stateout_b_stages[i] = integration(stateinitial_stages[i], stages_info, boosters_info, m_construction_stages[i], stages_area_pf[i], stages_area_bf[i], Cd_of_crosflow_cylinder, t_vertical_flight, theta, Azimuth, t_burn_stages[i], T_mag_stages[i], mass_flow_stages[i], has_boosters=False)

                xout_stages[i] = stateout_stages[i][:, 0]
                yout_stages[i] = stateout_stages[i][:, 1]
                zout_stages[i] = stateout_stages[i][:, 2]

                velxout_stages[i] = stateout_stages[i][:, 3]
                velyout_stages[i] = stateout_stages[i][:, 4]
                velzout_stages[i] = stateout_stages[i][:, 5]

                velmag_stages[i] = np.sqrt(velxout_stages[i]**2 + velyout_stages[i]**2 + velzout_stages[i]**2)

                xout_b_stages[i] = stateout_b_stages[i][:, 0]
                yout_b_stages[i] = stateout_b_stages[i][:, 1]
                zout_b_stages[i] = stateout_b_stages[i][:, 2]

                velxout_b_stages[i] = stateout_b_stages[i][:, 3]
                velyout_b_stages[i] = stateout_b_stages[i][:, 4]
                velzout_b_stages[i] = stateout_b_stages[i][:, 5]

                massout_stages[i] = stateout_b_stages[i][:, 6]
        else:
            if i == 0:
                stateinitial_stages[i] = np.array([x0, y0, z0, velx0, vely0, velz0, m0_stages[i]])
                
                tout_stages[i], stateout_stages[i], tout_b_stages[i], stateout_b_stages[i] = integration(stateinitial_stages[i], stages_info, boosters_info, m_construction_stages[i], stages_area_pf[i], stages_area_bf[i], Cd_of_crosflow_cylinder, t_vertical_flight, theta, Azimuth, t_burn_stages[i], T_mag_stages[i], mass_flow_stages[i], has_boosters=False)

                xout_stages[i] = stateout_stages[i][:, 0]
                yout_stages[i] = stateout_stages[i][:, 1]
                zout_stages[i] = stateout_stages[i][:, 2]

                velxout_stages[i] = stateout_stages[i][:, 3]
                velyout_stages[i] = stateout_stages[i][:, 4]
                velzout_stages[i] = stateout_stages[i][:, 5]

                velmag_stages[i] = np.sqrt(velxout_stages[i]**2 + velyout_stages[i]**2 + velzout_stages[i]**2)

                xout_b_stages[i] = stateout_b_stages[i][:, 0]
                yout_b_stages[i] = stateout_b_stages[i][:, 1]
                zout_b_stages[i] = stateout_b_stages[i][:, 2]

                velxout_b_stages[i] = stateout_b_stages[i][:, 3]
                velyout_b_stages[i] = stateout_b_stages[i][:, 4]
                velzout_b_stages[i] = stateout_b_stages[i][:, 5]

                massout_stages[i] = stateout_b_stages[i][:, 6]

            else:
                stateinitial_stages[i] = np.array([xout_b_stages[i-1][-1], yout_b_stages[i-1][-1], zout_b_stages[i-1][-1], velxout_b_stages[i-1][-1], velyout_b_stages[i-1][-1], velzout_b_stages[i-1][-1], m0_stages[i]])

                tout_stages[i], stateout_stages[i], tout_b_stages[i], stateout_b_stages[i] = integration(stateinitial_stages[i], stages_info, boosters_info, m_construction_stages[i], stages_area_pf[i], stages_area_bf[i], Cd_of_crosflow_cylinder, t_vertical_flight, theta, Azimuth, t_burn_stages[i], T_mag_stages[i], mass_flow_stages[i], has_boosters=False)

                xout_stages[i] = stateout_stages[i][:, 0]
                yout_stages[i] = stateout_stages[i][:, 1]
                zout_stages[i] = stateout_stages[i][:, 2]

                velxout_stages[i] = stateout_stages[i][:, 3]
                velyout_stages[i] = stateout_stages[i][:, 4]
                velzout_stages[i] = stateout_stages[i][:, 5]

                velmag_stages[i] = np.sqrt(velxout_stages[i]**2 + velyout_stages[i]**2 + velzout_stages[i]**2)

                xout_b_stages[i] = stateout_b_stages[i][:, 0]
                yout_b_stages[i] = stateout_b_stages[i][:, 1]
                zout_b_stages[i] = stateout_b_stages[i][:, 2]

                velxout_b_stages[i] = stateout_b_stages[i][:, 3]
                velyout_b_stages[i] = stateout_b_stages[i][:, 4]
                velzout_b_stages[i] = stateout_b_stages[i][:, 5]

                massout_stages[i] = stateout_b_stages[i][:, 6]

    if has_boosters:
        for i in range(booster_count):
            stateinitial_boosters[i] = np.array([x0, y0, z0, velx0, vely0, velz0, m0_each_boosters[i]])

            tout_boosters[i], stateout_boosters[i], tout_b_boosters[i], stateout_b_boosters[i] = integration(stateinitial_boosters[i], stages_info, boosters_info, m_construction_each_boosters[i], boosters_area_pf[i], boosters_area_bf[i], Cd_of_crosflow_cylinder, t_vertical_flight, theta, Azimuth, t_burn_boosters, T_mag_boosters[i], mass_flow_boosters[i], has_boosters=False)
            
            xout_boosters[i] = stateout_boosters[i][:, 0]
            yout_boosters[i] = stateout_boosters[i][:, 1]
            zout_boosters[i] = stateout_boosters[i][:, 2]

            velxout_boosters[i] = stateout_boosters[i][:, 3]
            velyout_boosters[i] = stateout_boosters[i][:, 4]
            velzout_boosters[i] = stateout_boosters[i][:, 5]

            velmag_boosters[i] = np.sqrt(velxout_boosters[i]**2 + velyout_boosters[i]**2 + velzout_boosters[i]**2)

            xout_b_boosters[i] = stateout_b_boosters[i][:, 0]
            yout_b_boosters[i] = stateout_b_boosters[i][:, 1]
            zout_b_boosters[i] = stateout_b_boosters[i][:, 2]

            velxout_b_boosters[i] = stateout_b_boosters[i][:, 3]
            velyout_b_boosters[i] = stateout_b_boosters[i][:, 4]
            velzout_b_boosters[i] = stateout_b_boosters[i][:, 5]

            massout_boosters[i] = stateout_b_boosters[i][:, 6]

    return tout_stages, massout_stages, stages_count