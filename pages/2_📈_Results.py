import os
import streamlit as st
from data_manager import load_data
from sidebar import sidebar
import numpy as np
from scipy.interpolate import interp1d
import pandas as pd
import plotly.express as px
##############################################################


current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
cd_file_path = os.path.join(project_root, 'resources', 'CD-Mach_relation.txt')
temperature_file_path = os.path.join(project_root, 'resources', 'temperature_profile_0_to_600km.txt')

##############################################################
if not os.path.exists(cd_file_path):
    st.error("Error: 'CD-Mach_relation.txt' file not found in resources folder.")
    st.stop()
else:
    cd_file = np.loadtxt(cd_file_path)
    Mach_ref = cd_file[:, 0]
    Cd_ref = cd_file[:, 1]

    drag_interp = interp1d(Mach_ref, Cd_ref, kind='linear')
    new_Mach = np.linspace(min(Mach_ref), max(Mach_ref), 200)
    new_Cd = drag_interp(new_Mach)
    cd_data = pd.DataFrame({'Mach': new_Mach, 'Cd': new_Cd})
##############################################################

##############################################################
if not os.path.exists(temperature_file_path):
    st.error("Error: 'temperature_profile_0_to_600km.txt' file not found in resources folder.")
    st.stop()
else:
    temp_file = np.loadtxt(temperature_file_path)
    alt_ref = temp_file[:, 0]
    temp_ref = temp_file[:, 1]
    temp_interp = interp1d(alt_ref, temp_ref, kind = 'linear', fill_value='extrapolate')

    new_alt = np.linspace(min(alt_ref), max(alt_ref), 200)
    new_temp = temp_interp(new_alt)
    temp_data = pd.DataFrame({'Temperature (K)': new_temp, 'Altitude (km)': new_alt/1000})
##############################################################



st.title("📈 Results")
selected_plot = sidebar("results_page")


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


##############################################################
if selected_plot == "3D Orbit":
    # plot_3d_orbit(results_data)
    pass
elif selected_plot == "Velocity vs Time":
    tab_velocity_m_s, tab_velocity_mach = st.tabs(["Velocity (m/s)", "Velocity (Mach)"])
    with tab_velocity_m_s:
        st.header("Velocity vs Time (m/s)")
        # Code for plotting velocity vs time in m/s
    with tab_velocity_mach:
        st.header("Velocity vs Time (Mach)")
        # Code for plotting velocity vs time in Mach numbers
    # plot_velocity_vs_time(results_data)
elif selected_plot == "Altitude vs Time":
    pass
    # plot_altitude_vs_time(results_data)
elif selected_plot == "Mass vs Time":
    # plot_mass_vs_time(results_data)
    st.write("Plotting Mass vs Time...")
elif selected_plot == "Density vs Altitude":
    # plot_density_vs_altitude(results_data)
    st.write("Plotting Density vs Altitude...")
elif selected_plot == "Temperature profile":
    plot_temperature_profile()
elif selected_plot == "Drag coefficient vs Mach":
    plot_drag_coefficient_vs_mach()
