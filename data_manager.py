import json
import os
import time
import streamlit as st
from streamlit.runtime.scriptrunner import get_script_run_ctx

# Folder to store session files
SESSION_DIR = "temp_sessions"
MAX_FILE_AGE_SECONDS = 3600  # 1 час

if not os.path.exists(SESSION_DIR):
    os.makedirs(SESSION_DIR)

# Default data structure
DEFAULT_DATA = {
    "payload_mass_without_booster": 2880.0, # Stoil
    "payload_mass_with_booster": 3880.0, # Stoil
    "theta_angle": 80.0,
    "stages_count": 2, # Stoil
    "has_boosters": True,
    "booster_count": 2,
    "t_burn_ratio": 0.71,
    "payload_mass_ratio_total": 0.03,
    "t_vertical_flight": 15.0,
    "input_mode": "EPS & lambda",
    "rocket_type": "Optimal",
    "stages_data_mass": [
        {"Start Mass (kg)": 1000.0, "Propellant (kg)": 1000.0, "Mass flow (kg/s)": 300.0, "Ve (m/s)": 3400.0},
        {"Start Mass (kg)": 1000.0, "Propellant (kg)": 1000.0, "Mass flow (kg/s)": 300.0, "Ve (m/s)": 3900.0}
    ],
    "boosters_data_mass": [
        {"Start Mass (kg)": 1555.56, "Propellant (kg)": 1444.44, "Mass flow (kg/s)": 322.22, "Ve (m/s)": 2900.00},
        {"Start Mass (kg)": 1555.56, "Propellant (kg)": 1444.44, "Mass flow (kg/s)": 322.22, "Ve (m/s)": 2900.00}
    ],
    "stages_data_eps": [
        {"EPS": 0.11, "Mass_flow (kg/s)": 170.0, "Ve (m/s)": 3400.0}, # EPS and Ve - Stoil
        {"EPS": 0.08, "Mass_flow (kg/s)": 170.0, "Ve (m/s)": 3900.0} # EPS and Ve - Stoil
    ],
    "boosters_data_eps": [
        {"EPS": 0.1, "Mass_flow (kg/s)": 322.22, "Ve (m/s)": 2900.00}, # EPS and Ve - Stoil
        {"EPS": 0.1, "Mass_flow (kg/s)": 322.22, "Ve (m/s)": 2900.00} # EPS and Ve - Stoil
    ],
    "launch_lat": 25.991389, "launch_lon": -97.183611, "launch_alt": 0.91,
    "launch_date": "2026-02-10", "launch_time": "07:15:00", # Stoil. Time UTC
    "orbit_a": 26571.0, "orbit_e": 0.0, "orbit_i": 55.0
}

def cleanup_old_sessions():
    """Delete session files older than MAX_FILE_AGE_SECONDS."""
    now = time.time()
    try:
        for filename in os.listdir(SESSION_DIR):
            file_path = os.path.join(SESSION_DIR, filename)
            # Check if file is older than MAX_FILE_AGE_SECONDS
            if os.path.getmtime(file_path) < now - MAX_FILE_AGE_SECONDS:
                os.remove(file_path)
    except Exception as e:
        print(f"Error during cleanup: {e}")

def get_session_id():
    """Get the current Streamlit session ID."""
    ctx = get_script_run_ctx()
    return ctx.session_id if ctx else "default"

def get_file_path():
    """Get the file path for the current session's data."""
    cleanup_old_sessions()
    return os.path.join(SESSION_DIR, f"{get_session_id()}.json")

def load_data():
    """Load data for the current session, or create default if not exists."""
    path = get_file_path()
    if not os.path.exists(path):
        save_data(DEFAULT_DATA) # Create default data file
        return DEFAULT_DATA
    with open(path, 'r') as f:
        return json.load(f)

def save_data(data):
    """Save data for the current session."""
    path = get_file_path()
    with open(path, 'w') as f:
        json.dump(data, f, indent=4)

def update_field(key, value):
    """Update a specific field in the session data."""
    data = load_data()
    data[key] = value
    save_data(data)