import os
import streamlit as st
from data_manager import load_data
from sidebar import sidebar


from Algorithm.functions import plot_3d_orbit, plot_velocity_vs_time, plot_altitude_vs_time, plot_mass_vs_time, plot_density_vs_altitude, plot_temperature_profile, plot_drag_coefficient_vs_mach

from Algorithm import main
##############################################################

st.title("📈 Results")


selected_option = sidebar("results_page")

##############################################################
if selected_option == "Parameters":
    pass
elif selected_option == "3D Orbit":
    # plot_3d_orbit(results_data)
    pass
elif selected_option == "Velocity vs Time":
    tab_velocity_m_s, tab_velocity_mach = st.tabs(["Velocity (m/s)", "Velocity (Mach)"])
    with tab_velocity_m_s:
        st.header("Velocity vs Time (m/s)")
        # Code for plotting velocity vs time in m/s
    with tab_velocity_mach:
        st.header("Velocity vs Time (Mach)")
        # Code for plotting velocity vs time in Mach numbers
    # plot_velocity_vs_time(results_data)
elif selected_option == "Altitude vs Time":
    pass
    # plot_altitude_vs_time(results_data)
elif selected_option == "Mass vs Time":
    # plot_mass_vs_time(results_data)
    st.write("Plotting Mass vs Time...")
elif selected_option == "Density vs Altitude":
    # plot_density_vs_altitude(results_data)
    st.write("Plotting Density vs Altitude...")
elif selected_option == "Temperature profile":
    plot_temperature_profile()
elif selected_option == "Drag coefficient vs Mach":
    plot_drag_coefficient_vs_mach()
