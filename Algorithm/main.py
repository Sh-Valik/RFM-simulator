from Algorithm.functions import (
    parameters_of_stages, parameters_of_boosters, geodetic_to_cartesian_WGS84,
    simulate_stage, simulate_two_burn_stage, extract_results, get_burnout_state,
    RocketStage, FlightConfig,
    utc_to_julian_date, compute_gmst, eci_greenwich_to_j2000,
    cartesian_to_keplerian, backward_propagation, find_kick_angle,
    compute_corrected_azimuth
)
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
        Ve_stages, mass_flow_stages, m0_stages, m_prop_stages, Vf_id_stages, Lambda_stages = parameters_of_stages(input_mode, stage_data_list, m_payload_without_boosters, payload_mass_ratio_total, rocket_type, stages_count)
        if has_boosters:
            booster_data_list = data["boosters_data_mass"]
            # placeholder

    else:
        stage_data_list = data["stages_data_eps"]
        Ve_stages, mass_flow_stages, diameter_stages, height_stages, m0_stages, m_prop_stages, Vf_id_stages, Lambda_stages = parameters_of_stages(input_mode, stage_data_list, m_payload_without_boosters, payload_mass_ratio_total, rocket_type, stages_count)
        if has_boosters:
            booster_data_list = data["boosters_data_eps"]
            Ve_boosters, mass_flow_boosters, diameter_boosters, height_boosters, m0_each_boosters, m_construction_each_boosters, m_prop_each_boosters, m0 = parameters_of_boosters(input_mode, booster_data_list, m_payload_without_boosters, m_payload_with_boosters, Vf_id_stages, m0_stages, m_prop_stages, stages_count, booster_count, Ve_stages, mass_flow_stages, t_burn_ratio)


    launch_lat = data["launch_lat"]
    launch_lon = data["launch_lon"]
    launch_alt = data["launch_alt"]
    launch_date_str = data["launch_date"]
    launch_time_str = data["launch_time"]

    t_o_a = data["orbit_a"]        # target SMA [km]
    t_o_e = data["orbit_e"]        # target eccentricity
    t_o_i = data["orbit_i"]        # target inclination [deg]

    # Kick angle: from data (degrees) or default
    kick_angle_deg = data.get("kick_angle", 0.5)  # default ~0.5°
    kick_angle_rad = np.radians(kick_angle_deg)
    kick_duration = data.get('kick_duration', 2.0)

    #######################################################
    # Compute GMST at launch time
    #######################################################
    date_parts = launch_date_str.split("-")
    time_parts = launch_time_str.split(":")
    launch_year = int(date_parts[0])
    launch_month = int(date_parts[1])
    launch_day = int(date_parts[2])
    launch_hour = int(time_parts[0])
    launch_minute = int(time_parts[1])
    launch_second = float(time_parts[2])

    jd_launch = utc_to_julian_date(launch_year, launch_month, launch_day,
                                     launch_hour, launch_minute, launch_second)
    gmst_launch = compute_gmst(jd_launch)
    #######################################################

    # Cross-sectional areas
    stages_area_pf = [np.pi * (diameter_stages[i]**2) / 4 for i in range(stages_count)]
    stages_area_bf = [diameter_stages[i] * height_stages[i] for i in range(stages_count)]
    
    if has_boosters:
        boosters_area_pf = [np.pi * (diameter_boosters[i]**2) / 4 for i in range(booster_count)]
        boosters_area_bf = [diameter_boosters[i] * height_boosters[i] for i in range(booster_count)]

    Cd_of_crosflow_cylinder = 1.25
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
    
    # Launch pad velocity due to Earth's rotation: v = ω × r
    omega_earth = 7.2921159e-5  # Earth rotation rate [rad/s]
    velx0 = -omega_earth * y0
    vely0 =  omega_earth * x0
    velz0 = 0.0

    # Compute corrected azimuth (accounting for Earth's rotation)
    h_target = (t_o_a * 1000.0) - Rplanet  # target altitude [m]
    az_data = compute_corrected_azimuth(
        lat_rad=np.deg2rad(launch_lat),
        inc_rad=np.deg2rad(t_o_i),
        h_target=h_target,
    )
    Azimuth = az_data['Az_corrected']  # corrected azimuth [rad]
    
    #######################################################
    # Build RocketStage objects
    #######################################################
    stage_objects = []
    for i in range(stages_count):
        stage_objects.append(RocketStage(
            T_mag=T_mag_stages[i],
            m_flow=mass_flow_stages[i],
            t_burn=t_burn_stages[i],
            m_construction=m_construction_stages[i],
            Area_pf=stages_area_pf[i],
            Area_bf=stages_area_bf[i],
        ))
    
    booster_objects = None
    if has_boosters:
        booster_objects = []
        for i in range(booster_count):
            booster_objects.append(RocketStage(
                T_mag=T_mag_boosters[i],
                m_flow=mass_flow_boosters[i],
                t_burn=t_burn_boosters,
                m_construction=m_construction_each_boosters[i],
                Area_pf=boosters_area_pf[i],
                Area_bf=boosters_area_bf[i],
            ))
    
    config = FlightConfig(
        stages=stage_objects,
        boosters=booster_objects,
        t_vertical=t_vertical_flight,
        kick_angle_rad=kick_angle_rad,
        Az_rad=Azimuth,
        Cd_crossflow=Cd_of_crosflow_cylinder,
        t_burn_boosters=t_burn_boosters if has_boosters else 0.0,
        kick_duration=kick_duration,
    )

    #######################################################
    # Run simulation for each stage
    #######################################################
    tout_stages = [None] * stages_count
    massout_stages = [None] * stages_count
    velmag_stages = [None] * stages_count
    stage_results = [None] * stages_count
    burnout_states = [None] * stages_count  # burnout state for each stage

    # Initial state for the first stage
    if has_boosters:
        initial_mass = m0[0]  # m0 from parameters_of_boosters (includes booster mass)
    else:
        initial_mass = m0_stages[0]

    state = np.array([x0, y0, z0, velx0, vely0, velz0, initial_mass])

    two_burn_result = None  # will be set for last stage
    fuel_reserve_fraction = data.get('fuel_reserve_fraction', 0.20)

    for i in range(stages_count):
        if i == stages_count - 1:
            # Last stage: use two-burn approach (transfer + circularization)
            two_burn_result = simulate_two_burn_stage(
                state, config, stage_index=i,
                fuel_reserve_fraction=fuel_reserve_fraction
            )
            stage_results[i] = two_burn_result
            tout_stages[i] = two_burn_result['t']
            extracted = extract_results(two_burn_result)
            massout_stages[i] = extracted['mass']
            velmag_stages[i] = extracted['velmag']
            # Final burnout state is after the second burn (circularization)
            burnout_states[i] = two_burn_result['burnout_state_2']
        else:
            result = simulate_stage(state, config, stage_index=i)
            stage_results[i] = result

            extracted = extract_results(result)
            tout_stages[i] = result['t']
            massout_stages[i] = extracted['mass']
            velmag_stages[i] = extracted['velmag']

            # Save burnout state for this stage
            burnout_states[i] = get_burnout_state(result)

            # Chain to next stage: use burnout state, replace mass with next stage mass
            state = get_burnout_state(result)
            if has_boosters:
                state[6] = m0[i + 1]
            else:
                state[6] = m0_stages[i + 1]

    #######################################################
    # Run simulation for each booster (FULL flight: launch → burn → ballistic)
    #######################################################
    tout_boosters = [None] * booster_count if has_boosters else []
    massout_boosters = [None] * booster_count if has_boosters else []
    velmag_boosters = [None] * booster_count if has_boosters else []

    if has_boosters:
        for i in range(booster_count):
            # Each booster is simulated as an independent single-stage rocket
            # from t=0 with its full mass (propellant + construction)
            booster_as_stage = RocketStage(
                T_mag=T_mag_boosters[i],
                m_flow=mass_flow_boosters[i],
                t_burn=t_burn_boosters,
                m_construction=m_construction_each_boosters[i],
                Area_pf=boosters_area_pf[i],
                Area_bf=boosters_area_bf[i],
            )

            booster_config = FlightConfig(
                stages=[booster_as_stage],
                boosters=None,           # no sub-boosters
                t_vertical=t_vertical_flight,
                kick_angle_rad=kick_angle_rad,
                Az_rad=Azimuth,
                Cd_crossflow=Cd_of_crosflow_cylinder,
                t_burn_boosters=0.0,
            )

            # Booster starts from same launch position/velocity with its own full mass
            booster_state = np.array([x0, y0, z0, velx0, vely0, velz0, m0_each_boosters[i]])

            booster_result = simulate_stage(booster_state, booster_config, stage_index=0)

            booster_extracted = extract_results(booster_result)
            tout_boosters[i] = booster_result['t']
            massout_boosters[i] = booster_extracted['mass']
            velmag_boosters[i] = booster_extracted['velmag']

    #######################################################
    # Convert final state (after circularization) to J2000 and Keplerian
    #######################################################
    last_burnout = burnout_states[stages_count - 1]
    
    # Compute total flight time: sum all stage burn times + coast + burn2
    t_total_flight = 0.0
    for i in range(stages_count - 1):
        t_total_flight += t_burn_stages[i]
    if two_burn_result is not None:
        pt = two_burn_result['phase_times']
        t_total_flight += pt['t_burn1'] + pt['t_coast'] + pt['t_burn2']
    else:
        t_total_flight += t_burn_stages[stages_count - 1]
    
    # GMST at end of flight = GMST at launch + Earth rotation during flight
    gmst_final = gmst_launch + omega_earth * t_total_flight
    
    # Convert final state from ECI Greenwich to J2000
    final_j2000 = eci_greenwich_to_j2000(last_burnout, gmst_final)
    
    # Convert to Keplerian elements
    keplerian = cartesian_to_keplerian(final_j2000)

    #######################################################
    # Return results
    #######################################################
    stages_return = [tout_stages, massout_stages, m0_stages, t_burn_stages, m_prop_stages, m_construction_stages, burnout_states]
    
    if has_boosters:
        boosters_return = [tout_boosters, massout_boosters, m0, t_burn_boosters, m_prop_each_boosters, m_construction_each_boosters]
    else:
        boosters_return = [[], [], [], 0, [], []]

    orbit_info = {
        'keplerian': keplerian,
        'final_j2000': final_j2000,
        'gmst_launch_rad': gmst_launch,
        'gmst_final_rad': gmst_final,
        'jd_launch': jd_launch,
        'kick_angle_deg': kick_angle_deg,
        'azimuth_geometric_rad': az_data['Az_geometric'],
        'azimuth_corrected_rad': Azimuth,
        'azimuth_data': az_data,
        'two_burn': two_burn_result['phase_times'] if two_burn_result else None,
    }

    return stages_return, boosters_return, orbit_info