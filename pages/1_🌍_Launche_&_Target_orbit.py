import streamlit as st
import numpy as np
import plotly.graph_objects as go
from sidebar import sidebar
from data_manager import load_data, update_field
import pandas as pd

sidebar("input_page")
data = load_data()



EARTH_RADIUS = 6371.0  # км




st.set_page_config(layout="wide")

st.title("🌍 Launch & Target orbit")

tab_launch, tab_orbit = st.tabs(["🚀 Launch Parameters", "🛰️ Target Orbit Parameters"])
with tab_launch:
    st.header("Launch Parameters")
    col1, col2 = st.columns(2)
    with col1:
        lat = st.number_input("Latitude (deg)", -90.0, 90.0, value=data["launch_lat"], on_change=lambda: update_field("launch_lat", lat), format="%.6f", key="launch_lat")
        lon = st.number_input("Longitude (deg)", -180.0, 180.0, value=data["launch_lon"], on_change=lambda: update_field("launch_lon", lon), format="%.6f", key="launch_lon")
        alt = st.number_input("Altitude (m)", 0.0, 1000.0, value=data["launch_alt"], step=0.1, on_change=lambda: update_field("launch_alt", alt), key="launch_alt")
    with col2:
        ldate = st.date_input("Launch Date", value=pd.to_datetime(data["launch_date"]).date(), on_change=lambda: update_field("launch_date", ldate.strftime("%Y-%m-%d")), key="launch_date")
        ltime = st.time_input("Launch Time", value=pd.to_datetime(data["launch_time"]).time(), on_change=lambda: update_field("launch_time", ltime.strftime("%H:%M:%S")), key="launch_time")

with tab_orbit:
    st.header("Target Orbit Parameters (Keplerian Elements)")
    col1, col2, col3 = st.columns(3)
    with col1:
        a = st.number_input("Semi-major axis (a) (km)", 0.0, value=data["orbit_a"], on_change=lambda: update_field("orbit_a", a), step=100.0, key="orbit_a")
    with col2:
        e = st.number_input("Eccentricity (e)", 0.0, 1.0, value=data["orbit_e"], on_change=lambda: update_field("orbit_e", e), step=0.01, key="orbit_e")
    with col3:
        i = st.number_input("Inclination (i) (deg)", 0.0, 360.0, value=data["orbit_i"], on_change=lambda: update_field("orbit_i", i), key="orbit_i")
