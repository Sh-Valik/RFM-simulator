import matplotlib.pyplot as plt
import numpy as np
from Algorithm.main import run_simulation

DEFAULT_DATA = {
    "payload_mass_without_booster": 2880.0, # Stoil
    "payload_mass_with_booster": 3880.0, # Stoil
    #"theta_angle": 80.0,
    "stages_count": 2, # Stoil
    "has_boosters": True,
    "booster_count": 2,
    "t_burn_ratio": 0.71,
    "payload_mass_ratio_total": 0.03,
    "t_vertical_flight": 5.0,
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
        {"EPS": 0.11, "Mass_flow (kg/s)": 250.0, "Ve (m/s)": 3400.0, "Diameter (m)": 3.7, "Height (m)": 42.6}, # EPS and Ve - Stoil
        {"EPS": 0.08, "Mass_flow (kg/s)": 100.0, "Ve (m/s)": 3900.0, "Diameter (m)": 3.7, "Height (m)": 12.6} # EPS and Ve - Stoil
    ],
    "boosters_data_eps": [
        {"EPS": 0.1, "Mass_flow (kg/s)": 250, "Ve (m/s)": 2900.00, "Diameter (m)": 3.7, "Height (m)": 42.6}, # EPS and Ve - Stoil
        {"EPS": 0.1, "Mass_flow (kg/s)": 250, "Ve (m/s)": 2900.00, "Diameter (m)": 3.7, "Height (m)": 42.6} # EPS and Ve - Stoil
    ],
    "launch_lat": 25.991389, "launch_lon": -97.183611, "launch_alt": 0.91,
    "launch_date": "2026-02-10", "launch_time": "07:15:00", # Stoil. Time UTC
    "orbit_a": 26571.0, "orbit_e": 0.0, "orbit_i": 55.0
}

stages_return, boosters_return, orbital_elements = run_simulation(DEFAULT_DATA)
print("Orbital elements:", orbital_elements)

tout_stages, massout_stages, xout_stages, yout_stages, zout_stages = stages_return
tout_boosters, massout_boosters, xout_boosters, yout_boosters, zout_boosters = boosters_return



# line_type = ['b-', 'c-', 'g-', 'r-', 'm-', 'y-']
# plt.figure()

# for i in range(len(tout_stages)):
#     plt.plot(tout_stages[i], massout_stages[i], label=f"Mass of stage{i+1}")
#     # plt.plot(tout_boosters[i], massout_boosters[i], line_type[i+3], label=f"Mass of boostes{i+1}")

# plt.grid()
# plt.legend()
# plt.show()

plt.figure()

for i in range(len(tout_boosters)):
    plt.plot(tout_boosters[i], massout_boosters[i], label=f"Mass of booster{i+1}")

plt.grid()
plt.legend()
plt.show()


G = 6.6742 * 10**-11  # gravitational constant [N.m^2/kg^2]
g0 = 9.80665  # standard gravitational acceleration [m/s^2]
Rplanet = 6371000  # mean radius of the Earth [m]
Mplanet = 5.97219 * 10**24  # mass of the Earth [kg]


# plt.figure()
fig = plt.figure('3D trajectory')
ax = fig.add_subplot(111, projection = '3d')
u, v_ = np.mgrid[0:2 * np.pi:50j, 0:np.pi:25j]
x_sphere = Rplanet * np.cos(u) * np.sin(v_)
y_sphere = Rplanet * np.sin(u) * np.sin(v_)
z_sphere = Rplanet * np.cos(v_)
ax.plot_surface(x_sphere, y_sphere, z_sphere, color = 'lightblue', alpha = 0.3)
ax.set_box_aspect([1, 1, 1])

for i in range(len(tout_stages)):
    ax.plot(xout_stages[i], yout_stages[i], zout_stages[i], label=f"Position of stage{i+1}")
for i in range(len(tout_boosters)):
    ax.plot(xout_boosters[i], yout_boosters[i], zout_boosters[i], label=f"Position of booster{i+1}")
ax.axis('equal')
ax.set_box_aspect([1, 1, 1])
ax.legend()
plt.show()







