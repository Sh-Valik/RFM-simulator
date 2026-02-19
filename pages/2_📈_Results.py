import os
import streamlit as st
from sidebar import sidebar
from Algorithm.functions import print_output_parameters, plot_3d_orbit, projection_test, plot_velocity_vs_time, plot_altitude_vs_time, plot_mass_vs_time, plot_temperature_profile, plot_drag_coefficient_vs_mach
from data_manager import load_data
from Algorithm.main import run_simulation
##############################################################

st.title("📈 Results")
selected_option = sidebar("results_page")
# tout_stages, velmag_stages, stages_count = run_simulation(load_data())
simulation_results = run_simulation(load_data())

trajectoiries = simulation_results["trajectories"]
stages_count = simulation_results["stages_count"]
booster_count = simulation_results["booster_count"]
velmag_stages = simulation_results["velmag_stages"]
velmag_boosters = simulation_results["velmag_boosters"]
tout_stages = simulation_results["tout_stages"]
tout_boosters = simulation_results["tout_boosters"]
altititude_stages = simulation_results["altitude_stages"]
altitude_boosters = simulation_results["altitude_boosters"]
massout_stages = simulation_results["massout_stages"]
massout_boosters = simulation_results["massout_boosters"]
orbital_elements = simulation_results["orbital_elements"]
rocket_parameters = simulation_results["rocket_parameters"]


##############################################################
if selected_option == "Parameters":
    print_output_parameters(rocket_parameters, orbital_elements)
elif selected_option == "3D Orbit":
    plot_3d_orbit(trajectoiries, stages_count, booster_count)
elif selected_option == "2D Projection":
    projection_test(trajectoiries, stages_count)
elif selected_option == "Velocity vs Time":
    st.header("Velocity vs Time")
    plot_velocity_vs_time(tout_stages, velmag_stages, stages_count, tout_boosters, velmag_boosters)
elif selected_option == "Altitude vs Time":
    plot_altitude_vs_time(tout_stages, altititude_stages, stages_count, tout_boosters, altitude_boosters)
elif selected_option == "Mass vs Time":
    plot_mass_vs_time(tout_stages, massout_stages, stages_count)
elif selected_option == "Temperature profile":
    plot_temperature_profile()
elif selected_option == "Drag coefficient vs Mach":
    plot_drag_coefficient_vs_mach()
