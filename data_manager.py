import json
import os
import time
import streamlit as st
from streamlit.runtime.scriptrunner import get_script_run_ctx

# Папка для хранения сессий
SESSION_DIR = "temp_sessions"
MAX_FILE_AGE_SECONDS = 3600  # 1 час

if not os.path.exists(SESSION_DIR):
    os.makedirs(SESSION_DIR)

# Дефолтные данные
DEFAULT_DATA = {
    "payload_mass": 1000.0,
    "stages_count": 2,
    "has_boosters": False,
    "booster_count": 1,
    "stages_data": [
        {"Propellant (kg)": 1000.0, "Mass flow (kg/s)": 300.0, "Ve (m/s)": 300.0, "Structural Mass (kg)": 100.0},
        {"Propellant (kg)": 1000.0, "Mass flow (kg/s)": 300.0, "Ve (m/s)": 300.0, "Structural Mass (kg)": 100.0}
    ],
    "boosters_data": [
        {"Propellant (kg)": 1000.0, "Mass flow (kg/s)": 300.0, "Ve (m/s)": 300.0, "Structural Mass (kg)": 100.0}
    ],
    "launch_lat": 0.0, "launch_lon": 0.0, "launch_alt": 0.0,
    "launch_date": "2025-02-10", "launch_time": "00:00:00",
    "orbit_a": 7000.0, "orbit_e": 0.0, "orbit_i": 0.0, "orbit_Ω": 0.0, "orbit_ω": 0.0, "orbit_ν": 0.0
}

def cleanup_old_sessions():
    """Удаляет файлы сессий, которые старше чем MAX_FILE_AGE_SECONDS"""
    now = time.time()
    try:
        for filename in os.listdir(SESSION_DIR):
            file_path = os.path.join(SESSION_DIR, filename)
            # Проверяем время последнего изменения файла
            if os.path.getmtime(file_path) < now - MAX_FILE_AGE_SECONDS:
                os.remove(file_path)
    except Exception as e:
        print(f"Error during cleanup: {e}")

def get_session_id():
    ctx = get_script_run_ctx()
    return ctx.session_id if ctx else "default"

def get_file_path():
    cleanup_old_sessions()
    return os.path.join(SESSION_DIR, f"{get_session_id()}.json")

def load_data():
    path = get_file_path()
    if not os.path.exists(path):
        save_data(DEFAULT_DATA) # Создаем файл с дефолтами при первом обращении
        return DEFAULT_DATA
    with open(path, 'r') as f:
        return json.load(f)

def save_data(data):
    path = get_file_path()
    with open(path, 'w') as f:
        json.dump(data, f, indent=4)

def update_field(key, value):
    data = load_data()
    data[key] = value
    save_data(data)