import streamlit as st
import numpy as np
import plotly.graph_objects as go
from sidebar import sidebar

sidebar("input_page")


EARTH_RADIUS = 6371.0  # км

def get_orbit_coords(a, e, i, Omega, omega, n_points=200):
    """Вычисляет координаты точек орбиты"""
    # 1. Создаем массив углов (истинная аномалия) от 0 до 2pi
    nu = np.linspace(0, 2 * np.pi, n_points)
    
    # 2. Вычисляем радиус для каждой точки (уравнение эллипса в полярных координатах)
    r = a * (1 - e**2) / (1 + e * np.cos(nu))
    
    # 3. Координаты в перифокальной системе (плоскость орбиты)
    x_p = r * np.cos(nu)
    y_p = r * np.sin(nu)
    z_p = np.zeros_like(nu)
    
    # 4. Матрицы поворота
    i_rad = np.radians(i)
    Omega_rad = np.radians(Omega)
    omega_rad = np.radians(omega)
    
    # Матрица перехода из перифокальной системы в экваториальную (ECI)
    # Поворот на omega вокруг Z, затем на i вокруг X, затем на Omega вокруг Z
    cos_O, sin_O = np.cos(Omega_rad), np.sin(Omega_rad)
    cos_i, sin_i = np.cos(i_rad), np.sin(i_rad)
    cos_w, sin_w = np.cos(omega_rad), np.sin(omega_rad)
    
    # Элементы матрицы трансформации
    X = (cos_O*cos_w - sin_O*sin_w*cos_i)*x_p + (-cos_O*sin_w - sin_O*cos_w*cos_i)*y_p
    Y = (sin_O*cos_w + cos_O*sin_w*cos_i)*x_p + (-sin_O*sin_w + cos_O*cos_w*cos_i)*y_p
    Z = (sin_w*sin_i)*x_p + (cos_w*sin_i)*y_p
    
    return X, Y, Z

def get_single_point(a, e, i, Omega, omega, nu_deg):
    """Вычисляет положение конкретной точки (ракеты)"""
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
        lat = st.number_input("Latitude (deg)", -90.0, 90.0, 0.0, format="%.2f")
        lon = st.number_input("Longitude (deg)", -180.0, 180.0, 0.0, format="%.2f")
        alt = st.number_input("Altitude (m)", 0.0, 10000.0, 0.0, step=10.0)
    with col2:
        launch_date = st.date_input("Launch Date")
        launch_time = st.time_input("Launch Time")

with tab_orbit:
    st.header("Target Orbit Parameters (Keplerian Elements)")
    col1, col2, col3 = st.columns(3)
    with col1:
        a = st.number_input("Semi-major axis (a) (km)", 0.0, value=7000.0, step=100.0)
        e = st.number_input("Eccentricity (e)", 0.0, 1.0, value=0.0, step=0.01)
    with col2:
        i = st.number_input("Inclination (i) (deg)", 0.0, 360.0, value=0.0)
        Ω = st.number_input("Right Ascension of Ascending Node (Ω) (deg)", 0.0, 360.0, value=0.0)
    with col3:
        ω = st.number_input("Argument of Periapsis (ω) (deg)", 0.0, 360.0, value=0.0)
        ν = st.number_input("True Anomaly (ν) (deg)", 0.0, 360.0, value=0.0)



    st.divider()
    st.header("Orbit Visualization")
    orbit_x, orbit_y, orbit_z = get_orbit_coords(a, e, i, Ω, ω)
    pos_x, pos_y, pos_z = get_single_point(a, e, i, Ω, ω, ν)

    # 3. Создание 3D сцены
    fig = go.Figure()

    # Рисуем Землю (сфера)
    u, v = np.mgrid[0:2*np.pi:30j, 0:np.pi:20j]
    earth_x = EARTH_RADIUS * np.cos(u) * np.sin(v)
    earth_y = EARTH_RADIUS * np.sin(u) * np.sin(v)
    earth_z = EARTH_RADIUS * np.cos(v)
    fig.add_trace(go.Surface(x=earth_x, y=earth_y, z=earth_z, colorscale='Blues', showscale=False, opacity=0.8, name='Earth'))

    # Рисуем орбиту
    fig.add_trace(go.Scatter3d(x=orbit_x, y=orbit_y, z=orbit_z, mode='lines', line=dict(color='yellow', width=5), name='Orbit Path'))

    # Рисуем ракету (текущее положение)
    # fig.add_trace(go.Scatter3d(x=[pos_x], y=[pos_y], z=[pos_z], mode='markers', marker=dict(size=8, color='red', symbol='diamond'), name='Current Position'))

    # Настройка осей и масштаба
    max_val = max(a * (1+e), EARTH_RADIUS * 1.5) # Масштаб по апоцентру
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
