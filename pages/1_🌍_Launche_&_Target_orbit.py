import streamlit as st
import numpy as np
import plotly.graph_objects as go
from sidebar import sidebar
from data_manager import load_data, update_field
import pandas as pd


sidebar("input_page")
data = load_data()



EARTH_RADIUS = 6371.0  # км

def get_orbit_coords(a, e, i, Omega, omega, n_points=200):
    """Calculate orbit coordinates in ECI frame"""
    # 1. Create an array of true anomalies
    nu = np.linspace(0, 2 * np.pi, n_points)
    
    # 2. Calculate radius for each point (ellipse equation in polar coordinates)
    r = a * (1 - e**2) / (1 + e * np.cos(nu))
    
    # 3. Coordinates in the perifocal system (orbit plane)
    x_p = r * np.cos(nu)
    y_p = r * np.sin(nu)
    z_p = np.zeros_like(nu)
    
    # 4. Rotation matrices
    i_rad = np.radians(i)
    Omega_rad = np.radians(Omega)
    omega_rad = np.radians(omega)
    
    # Matrix transition from perifocal system to equatorial (ECI)
    # Rotation on omega around Z, then on i around X, then on Omega around Z
    cos_O, sin_O = np.cos(Omega_rad), np.sin(Omega_rad)
    cos_i, sin_i = np.cos(i_rad), np.sin(i_rad)
    cos_w, sin_w = np.cos(omega_rad), np.sin(omega_rad)

    # Elements of the transformation matrix
    X = (cos_O*cos_w - sin_O*sin_w*cos_i)*x_p + (-cos_O*sin_w - sin_O*cos_w*cos_i)*y_p
    Y = (sin_O*cos_w + cos_O*sin_w*cos_i)*x_p + (-sin_O*sin_w + cos_O*cos_w*cos_i)*y_p
    Z = (sin_w*sin_i)*x_p + (cos_w*sin_i)*y_p
    
    return X, Y, Z

def get_single_point(a, e, i, Omega, omega, nu_deg):
    """Calculate position of a specific point (rocket)"""
    nu_rad = np.radians(nu_deg)
    r = a * (1 - e**2) / (1 + e * np.cos(nu_rad))
    x_p = r * np.cos(nu_rad)
    y_p = r * np.sin(nu_rad)
    
    i_rad, Omega_rad, omega_rad = np.radians(i), np.radians(Omega), np.radians(omega)
    cos_O, sin_O = np.cos(Omega_rad), np.sin(Omega_rad)
    cos_i, sin_i = np.cos(i_rad), np.sin(i_rad)
    cos_w, sin_w = np.cos(omega_rad), np.sin(omega_rad)
    
    x = (cos_O*cos_w - sin_O*sin_w*cos_i)*x_p + (-cos_O*sin_w - sin_O*cos_w*cos_i)*y_p
    y = (sin_O*cos_w + cos_O*sin_w*cos_i)*x_p + (-sin_O*sin_w + cos_O*cos_w*cos_i)*y_p
    z = (sin_w*sin_i)*x_p + (cos_w*sin_i)*y_p
    return x, y, z

st.set_page_config(layout="wide")

st.title("🌍 Launch & Target orbit")

tab_launch, tab_orbit = st.tabs(["🚀 Launch Parameters", "🛰️ Target Orbit Parameters"])
with tab_launch:
    st.header("Launch Parameters")
    col1, col2 = st.columns(2)
    with col1:
        lat = st.number_input("Latitude (deg)", -90.0, 90.0, value=data["launch_lat"], on_change=lambda: update_field("launch_lat", lat), format="%.2f", key="launch_lat")
        lon = st.number_input("Longitude (deg)", -180.0, 180.0, value=data["launch_lon"], on_change=lambda: update_field("launch_lon", lon), format="%.2f", key="launch_lon")
        alt = st.number_input("Altitude (m)", 0.0, 10000.0, value=data["launch_alt"], step=10.0, on_change=lambda: update_field("launch_alt", alt), key="launch_alt")
    with col2:
        ldate = st.date_input("Launch Date", value=pd.to_datetime(data["launch_date"]).date(), on_change=lambda: update_field("launch_date", ldate.strftime("%Y-%m-%d")), key="launch_date")
        ltime = st.time_input("Launch Time", value=pd.to_datetime(data["launch_time"]).time(), on_change=lambda: update_field("launch_time", ltime.strftime("%H:%M:%S")), key="launch_time")

with tab_orbit:
    st.header("Target Orbit Parameters (Keplerian Elements)")
    col1, col2, col3 = st.columns(3)
    with col1:
        a = st.number_input("Semi-major axis (a) (km)", 0.0, value=data["orbit_a"], on_change=lambda: update_field("orbit_a", a), step=100.0, key="orbit_a")
        e = st.number_input("Eccentricity (e)", 0.0, 1.0, value=data["orbit_e"], on_change=lambda: update_field("orbit_e", e), step=0.01, key="orbit_e")
    with col2:
        i = st.number_input("Inclination (i) (deg)", 0.0, 360.0, value=data["orbit_i"], on_change=lambda: update_field("orbit_i", i), key="orbit_i")
        Ω = st.number_input("Right Ascension of Ascending Node (Ω) (deg)", 0.0, 360.0, value=data["orbit_Ω"], on_change=lambda: update_field("orbit_Ω", Ω), key="orbit_Ω")
    with col3:
        ω = st.number_input("Argument of Periapsis (ω) (deg)", 0.0, 360.0, value=data["orbit_ω"], on_change=lambda: update_field("orbit_ω", ω), key="orbit_ω")
        ν = st.number_input("True Anomaly (ν) (deg)", 0.0, 360.0, value=data["orbit_ν"], on_change=lambda: update_field("orbit_ν", ν), key="orbit_ν")

    st.divider()
    st.header("Orbit Visualization")
    orbit_x, orbit_y, orbit_z = get_orbit_coords(a, e, i, Ω, ω)
    pos_x, pos_y, pos_z = get_single_point(a, e, i, Ω, ω, ν)

    # 3. Create 3D plot with Plotly
    fig = go.Figure()

    # Draw Earth
    u, v = np.mgrid[0:2*np.pi:30j, 0:np.pi:20j]
    earth_x = EARTH_RADIUS * np.cos(u) * np.sin(v)
    earth_y = EARTH_RADIUS * np.sin(u) * np.sin(v)
    earth_z = EARTH_RADIUS * np.cos(v)
    fig.add_trace(go.Surface(x=earth_x, y=earth_y, z=earth_z, colorscale='Blues', showscale=False, opacity=0.8, name='Earth'))

    # Draw orbit
    fig.add_trace(go.Scatter3d(x=orbit_x, y=orbit_y, z=orbit_z, mode='lines', line=dict(color='yellow', width=5), name='Orbit Path'))

    # Draw rocket (current position)
    # fig.add_trace(go.Scatter3d(x=[pos_x], y=[pos_y], z=[pos_z], mode='markers', marker=dict(size=8, color='red', symbol='diamond'), name='Current Position'))

    # Axis settings
    max_val = max(a * (1+e), EARTH_RADIUS * 1.5) # Scale according to orbit size
    fig.update_layout(
        scene=dict(
            xaxis=dict(title='X (km)', range=[-max_val, max_val]),
            yaxis=dict(title='Y (km)', range=[-max_val, max_val]),
            zaxis=dict(title='Z (km)', range=[-max_val, max_val]),
            aspectmode='manual',
            aspectratio=dict(x=1, y=1, z=1)
        ),
        margin=dict(l=0, r=0, b=0, t=0),
        height=600
    )

    st.plotly_chart(fig, use_container_width=True)
