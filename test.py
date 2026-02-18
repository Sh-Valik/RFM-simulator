import matplotlib.pyplot as plt
import numpy as np
from Algorithm.main import run_simulation
from scipy.optimize import minimize

DEFAULT_DATA = {
    # Для автоподбора:
    # t_vertical_flight, theta_angle_deg, fuel_reserve_fraction
    # Можно варьировать вручную или через цикл
    "payload_mass_without_booster": 2880.0, # Stoil
    "payload_mass_with_booster": 3880.0, # Stoil
    "stages_count": 2, # Stoil
    "has_boosters": True,
    "booster_count": 2,
    "t_burn_ratio": 0.71,
    "payload_mass_ratio_total": 0.005,
    "t_vertical_flight": 64.994,
    "theta_angle": 67.7,
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
        {"EPS": 0.11, "Mass_flow (kg/s)": 2700.0, "Ve (m/s)": 3400.0, "Diameter (m)": 3.7, "Height (m)": 42.6}, # EPS and Ve - Stoil
        {"EPS": 0.08, "Mass_flow (kg/s)": 500.0, "Ve (m/s)": 3900.0, "Diameter (m)": 3.7, "Height (m)": 12.6} # EPS and Ve - Stoil
    ],
    "boosters_data_eps": [
        {"EPS": 0.1, "Mass_flow (kg/s)": 2000, "Ve (m/s)": 2900.00, "Diameter (m)": 3.7, "Height (m)": 42.6}, # EPS and Ve - Stoil
        {"EPS": 0.1, "Mass_flow (kg/s)": 2000, "Ve (m/s)": 2900.00, "Diameter (m)": 3.7, "Height (m)": 42.6} # EPS and Ve - Stoil
    ],
    "launch_lat": 25.991389, "launch_lon": -97.183611, "launch_alt": 0.91,
    "launch_date": "2026-02-10", "launch_time": "07:15:00", # Stoil. Time UTC
    "orbit_a": 26571.0, "orbit_e": 0.0, "orbit_i": 55.0
}

return_data = run_simulation(DEFAULT_DATA)
orbital_elements = return_data['orbital_elements']
trajectories = return_data['trajectories']
stages_trajectories = trajectories[0]
boosters_trajectories = trajectories[1]
rocket_parameters = return_data['rocket_parameters']
boosters_parameters = rocket_parameters['boosters_parameters']
t_burn_boosters = boosters_parameters['t_burn_boosters']
tout_booster = np.linspace(0, 3000, 12000)
altitude_boosters = return_data['altitude_boosters']
tout_boosters = return_data['tout_boosters']
velmag_boosters = return_data['velmag_boosters']
Rplanet = 6371000


tout_stages = return_data['tout_stages'][-1]
idx1 = np.abs(tout_stages - 262).argmin() # bout 1
idx2 = np.abs(tout_stages - 10262).argmin() # apogee
idx3 = np.abs(tout_stages - (10262+6.21)).argmin() # bout 2
print('idx1', idx1)
print('idx2', idx2)
print('idx3', idx3)



# plt.figure()
# plt.plot(tout_boosters, velmag_boosters)
# plt.grid()
# plt.show()



# fig = plt.figure('3D trajectory')
# ax = fig.add_subplot(111, projection = '3d')
# u, v_ = np.mgrid[0:2 * np.pi:50j, 0:np.pi:25j]
# x_sphere = Rplanet * np.cos(u) * np.sin(v_)
# y_sphere = Rplanet * np.sin(u) * np.sin(v_)
# z_sphere = Rplanet * np.cos(v_)
# ax.plot_surface(x_sphere, y_sphere, z_sphere, color = 'lightblue', alpha = 0.3)
# ax.set_box_aspect([1, 1, 1])

# ax.plot(stages_trajectories[3][0], stages_trajectories[4][0], stages_trajectories[5][0], label="Stage 1")
# ax.plot(boosters_trajectories[3], boosters_trajectories[4], boosters_trajectories[5], label="Boosters")


# ax.axis('equal')
# ax.set_box_aspect([1, 1, 1])
# ax.legend()
# plt.show()





# tout_stages, massout_stages, xout_stages, yout_stages, zout_stages = stages_return

# altitude_stages = []
# for i in range(len(tout_stages)):
#     current_alt = np.sqrt(xout_stages[i]**2 + yout_stages[i]**2 + zout_stages[i]**2) - 6371000
#     altitude_stages.append(current_alt)

# plt.figure()
# for i in range(len(tout_stages)):
#     plt.plot(tout_stages[i], altitude_stages[i], label=f"Altitude of stage{i+1}")
# plt.grid()
# plt.legend()
# plt.show()


# line_type = ['b-', 'c-', 'g-', 'r-', 'm-', 'y-']
# plt.figure()

# for i in range(len(tout_stages)):
#     plt.plot(tout_stages[i], massout_stages[i], label=f"Mass of stage{i+1}")
#     # plt.plot(tout_boosters[i], massout_boosters[i], line_type[i+3], label=f"Mass of boostes{i+1}")

# plt.grid()
# plt.legend()
# plt.show()

# plt.figure()

# for i in range(len(tout_boosters)):
#     plt.plot(tout_boosters[i], massout_boosters[i], label=f"Mass of booster{i+1}")

# plt.grid()
# plt.legend()
# plt.show()


# G = 6.6742 * 10**-11  # gravitational constant [N.m^2/kg^2]
# g0 = 9.80665  # standard gravitational acceleration [m/s^2]
# Rplanet = 6371000  # mean radius of the Earth [m]
# Mplanet = 5.97219 * 10**24  # mass of the Earth [kg]


# fig = plt.figure('3D trajectory')
# ax = fig.add_subplot(111, projection = '3d')
# u, v_ = np.mgrid[0:2 * np.pi:50j, 0:np.pi:25j]
# x_sphere = Rplanet * np.cos(u) * np.sin(v_)
# y_sphere = Rplanet * np.sin(u) * np.sin(v_)
# z_sphere = Rplanet * np.cos(v_)
# ax.plot_surface(x_sphere, y_sphere, z_sphere, color = 'lightblue', alpha = 0.3)
# ax.set_box_aspect([1, 1, 1])

# for i in range(len(tout_stages)):
#     ax.plot(xout_stages[i], yout_stages[i], zout_stages[i], label=f"Position of stage{i+1}")
# ax.axis('equal')
# ax.set_box_aspect([1, 1, 1])
# ax.legend()
# plt.show()



def objective_function(params, data, target_orbit):
    """
    Функция стоимости (Cost function). 
    Оптимизатор вызывает её, меняя params, чтобы вернуть как можно меньшее число.
    """
    # 1. Распаковываем параметры, которые подбирает оптимизатор
    t_vertical_guess, theta_angle_guess = params

    # 2. Обновляем входные данные для симуляции
    # Мы создаем копию, чтобы не ломать исходный словарь
    sim_data = data.copy()
    sim_data["t_vertical_flight"] = t_vertical_guess
    sim_data["theta_angle"] = theta_angle_guess

    # 3. Запускаем симуляцию
    # Важно: run_simulation должна возвращать orbital_elements третьим аргументом!
    try:
        return_data = run_simulation(sim_data)
        orbital_elements = return_data['orbital_elements']
        
        # Полученные значения
        calc_a = orbital_elements['a']  # км
        calc_e = orbital_elements['e']  # безразмерный
        calc_i = orbital_elements['i']  # градусы

        # Целевые значения
        target_a = target_orbit['a']
        target_e = target_orbit['e']
        target_i = target_orbit['i']

        # 4. Считаем ошибку (штраф)
        # Нормализуем веса, так как 'a' измеряется тысячами км, а 'e' — долями единицы.
        
        # Штраф за большую полуось (в км)
        # Если a = inf (парабола/гипербола), даем огромный штраф
        if np.isinf(calc_a) or calc_a < 0:
            return 1e9

        error_a = ((calc_a - target_a) / target_a) ** 2  # Относительная квадратичная ошибка
        
        # Штраф за эксцентриситет (умножаем на вес, т.к. значение маленькое)
        error_e = (calc_e - target_e) ** 2 * 1000 
        
        # Штраф за наклонение (обычно зависит от азимута, но кик тоже влияет)
        error_i = (calc_i - target_i) ** 2 * 10 

        total_cost = error_a + error_e + error_i
        
        return total_cost

    except Exception as e:
        # Если симуляция упала (например, ракета врезалась в землю), возвращаем огромный штраф
        return 1e9


def find_optimal_parameters(initial_data):
    """
    Основная функция для поиска параметров.
    """
    print("Начинаем подбор параметров орбиты...")

    # Целевые параметры орбиты берем из исходного json/словаря
    target_orbit = {
        'a': initial_data["orbit_a"], # Ожидается в км
        'e': initial_data["orbit_e"],
        'i': initial_data["orbit_i"]
    }

    # Начальное предположение [t_vertical, theta_angle]
    # Берем то, что было в конфиге изначально
    x0 = [initial_data["t_vertical_flight"], 72.0] 

    # Границы поиска (Bounds):
    # t_vertical: от 1 сек до 30 сек (пример)
    # theta_angle: от 0 град (вертикально) до 89.9 град (горизонтально)
    # Примечание: theta_angle обычно отсчитывается от вертикали. 
    # Если у вас 0 - это горизонт, поменяйте границы.
    bounds = [(1.0, 50.0), (45.0, 89.0)] 

    # ЗАПУСК ОПТИМИЗАТОРА
    # Используем метод Nelder-Mead, так как он хорошо работает с негладкими функциями (симуляциями)
    # или 'L-BFGS-B' если нужны строгие границы.
    result = minimize(
        objective_function, 
        x0, 
        args=(initial_data, target_orbit),
        method='Nelder-Mead', 
        bounds=bounds if 'L-BFGS-B' in ['L-BFGS-B', 'TNC', 'SLSQP'] else None, # Nelder-Mead не всегда поддерживает bounds напрямую через этот интерфейс в старых версиях, но в новых да.
        tol=1e-4,
        options={'maxiter': 100, 'disp': True}
    )

    print("\n--- Результаты оптимизации ---")
    if result.success:
        print("Оптимизация успешна!")
    else:
        print("Оптимизатор завершил работу (возможно, локальный минимум).")

    best_t_vertical, best_theta_angle = result.x
    
    print(f"Оптимальное время вертикального полета: {best_t_vertical:.4f} с")
    print(f"Оптимальный угол (Kick Angle): {best_theta_angle:.4f} град")
    print(f"Финальная ошибка (cost): {result.fun:.6f}")

    # Запускаем финальную симуляцию с лучшими параметрами, чтобы получить графики
    print("\nЗапуск контрольной симуляции...")
    initial_data["t_vertical_flight"] = best_t_vertical
    initial_data["theta_angle"] = best_theta_angle
    
    return_data = run_simulation(initial_data)
    elements = return_data['orbital_elements']
    
    print(f"Полученная орбита:\n SMA (a): {elements['a']:.2f} km\n ECC (e): {elements['e']:.4f}\n INC (i): {elements['i']:.2f} deg")

    return best_t_vertical, best_theta_angle, elements

# best_t, best_angle, final_orbit = find_optimal_parameters(DEFAULT_DATA)