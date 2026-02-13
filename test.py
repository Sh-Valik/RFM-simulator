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
        {"EPS": 0.11, "Mass_flow (kg/s)": 2000.0, "Ve (m/s)": 3400.0, "Diameter (m)": 3.7, "Height (m)": 42.6}, # EPS and Ve - Stoil
        {"EPS": 0.08, "Mass_flow (kg/s)": 270.0, "Ve (m/s)": 3900.0, "Diameter (m)": 3.7, "Height (m)": 12.6} # EPS and Ve - Stoil
    ],
    "boosters_data_eps": [
        {"EPS": 0.1, "Mass_flow (kg/s)": 2000, "Ve (m/s)": 2900.00, "Diameter (m)": 3.7, "Height (m)": 42.6}, # EPS and Ve - Stoil
        {"EPS": 0.1, "Mass_flow (kg/s)": 2000, "Ve (m/s)": 2900.00, "Diameter (m)": 3.7, "Height (m)": 42.6} # EPS and Ve - Stoil
    ],
    "launch_lat": 25.991389, "launch_lon": -97.183611, "launch_alt": 0.91,
    "launch_date": "2026-02-10", "launch_time": "07:15:00", # Stoil. Time UTC
    "orbit_a": 26571.0, "orbit_e": 0.0, "orbit_i": 55.0
}

stages_return, boosters_return = run_simulation(DEFAULT_DATA)

tout_stages, massout_stages, m0_stages, t_burn_stages, m_prop_stages, m_construction_stages = stages_return
tout_boosters, massout_boosters, m0, t_burn_boosters, m_prop_each_boosters, m_construction_each_boosters = boosters_return

print(f'mo_stages: {m0_stages}')
print(f't_burn_stages: {t_burn_stages}')
print(f'm_prop_stages: {m_prop_stages}')
print(f'm_construction_stages: {m_construction_stages}')


print(f'm0: {m0}')
print(f't_burn_boosters: {t_burn_boosters}')
print(f'm_prop_each_boosters: {m_prop_each_boosters}')
print(f'm_construction_each_boosters: {m_construction_each_boosters}')









line_type = ['b-', 'c-', 'g-', 'r-', 'm-', 'y-']
plt.figure()

for i in range(len(tout_stages)):
    plt.plot(tout_stages[i], massout_stages[i], line_type[i], label=f"Mass of stage{i+1}")
    # plt.plot(tout_boosters[i], massout_boosters[i], line_type[i+3], label=f"Mass of boostes{i+1}")

plt.grid()
plt.legend()
plt.show()