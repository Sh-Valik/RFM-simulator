import streamlit as st
import pandas as pd
from sidebar import sidebar
from data_manager import load_data, update_field

st.set_page_config(layout="wide")
sidebar("input_page")

data = load_data()

def counter_sync(count_key, category):
    """
    count_key: ключ session_state с числом (например, "stages_count")
    category: префикс данных ("stages" или "boosters")
    """
    new_count = st.session_state[count_key]
    current_data = load_data()

    # Определяем текущий список данных
    mode = current_data.get("input_mode", "EPS & lambda")

    if mode == "Start mass & Propellant":
        data_key = f"{category}_data_mass"
        new_row = {"Start Mass (kg)": 1000.0, "Propellant (kg)": 1000.0, "Mass flow (kg/s)": 300.0, "Ve (m/s)": 300.0}
    else:
        data_key = f"{category}_data_eps"
        new_row = {"EPS": 0.08, "lambda": 0.03, "Mass_flow (kg/s)": 170.0, "Ve (m/s)": 300.0}
    

    current_list = current_data.get(data_key, [])

    if len(current_list) < new_count:
        # Добавляем новые строки
        for _ in range(new_count - len(current_list)):
            current_list.append(new_row.copy())
    elif len(current_list) > new_count:
        # Удаляем лишние строки
        current_list = current_list[:new_count]

    # Сохраняем обновленные данные
    update_field(data_key, current_list)
    update_field(count_key, new_count)
    

######################
st.title("🚀 Rocket Parameters")
st.header("General Parameters")
col1, col2, col3 = st.columns(3)

with col1:
    st.number_input("Payload Mass (kg)", min_value=0.0, value=data["payload_mass"], step=100.0, key="payload_mass", on_change=lambda: update_field("payload_mass", st.session_state.payload_mass), help="Mass of the payload to be delivered to orbit.")
    st.checkbox("Has Boosters", value=data["has_boosters"], key="has_boosters", on_change=lambda: update_field("has_boosters", st.session_state.has_boosters), help="Check if the rocket has additional boosters for extra thrust.")



    input_mode = ["EPS & lambda", "Start mass & Propellant"]
    current_mode_value = data.get("input_mode", "EPS & lambda")
    try:
        mode_index = input_mode.index(current_mode_value)
    except ValueError:
        mode_index = 0

    st.radio("Rocket parameters input mode", options=input_mode, index=mode_index, key="input_mode", on_change=lambda: update_field("input_mode", st.session_state.input_mode), help="Select the method for inputting rocket parameters. EPS is construction mass ratio, lambda is payload mass ratio.")


with col2:
    st.number_input("Theta angle (degrees)", min_value=0.0, max_value=90.0, value=data["theta_angle"], step=1.0, key="theta_angle", on_change=lambda: update_field("theta_angle", st.session_state.theta_angle), help="Launch angle of the rocket relative to the horizontal plane.")
    if data["has_boosters"]:
        st.number_input("Number of rocket boosters", min_value=1, max_value=10, value=data["booster_count"], step=1, key="booster_count", on_change=lambda: counter_sync("booster_count", "boosters"), help="Number of additional boosters attached to the rocket for extra thrust.")

    if st.session_state.input_mode == "EPS & lambda":
        r_types = ["Optimal", "Non-optimal"]
        r_idx = r_types.index(data.get("rocket_type", "Optimal")) if data.get("rocket_type", "Optimal") in r_types else 0
        st.radio("Type of a Rocket", options=r_types, index=r_idx, key="rocket_type", on_change=lambda: update_field("rocket_type", st.session_state.rocket_type), help="Select whether the rocket is designed for optimal or non-optimal performance.")

        
with col3:
    st.number_input("Number of rocket stages", min_value=1, max_value=10, value=data["stages_count"], step=1, key="stages_count", on_change=lambda: counter_sync("stages_count", "stages"), help="Number of stages in the rocket. Each stage is jettisoned when its fuel is expended.")
    

    if st.session_state.input_mode == "EPS & lambda":
        st.number_input("Total payload mass ratio", min_value=0.0, max_value=1.0, value=data["payload_mass_ratio_total"], step=0.01, key="payload_mass_ratio", on_change=lambda: update_field("payload_mass_ratio", st.session_state.payload_mass_ratio), help="Overall payload mass ratio of the rocket. Optimal rocket have optimal payload mass ratio distribution. Non-optimal rocket have equal payload mass ratio for each stage.")


st.divider()
st.header("Stages Details")


current_mode = data.get("input_mode", "EPS & lambda")
key_suffix = "mass" if current_mode == "Start mass & Propellant" else "eps"

stage_real_key = f"stages_data_{key_suffix}"
df_stages = pd.DataFrame(data.get(stage_real_key, []))
df_stages.index += 1  # Начинаем нумерацию с 1
df_stages.index.name = "Stage"

edited_stages = st.data_editor(df_stages, use_container_width=True, key="ed_stages")

# Сохраняем таблицу, если она изменилась
if st.session_state.ed_stages:
    update_field(stage_real_key, edited_stages.to_dict('records'))

if data["has_boosters"]:
    st.divider()
    st.header("Boosters Details")

    boosters_real_key = f"boosters_data_{key_suffix}"

    df_boosters = pd.DataFrame(data.get(boosters_real_key, []))
    df_boosters.index += 1  # Начинаем нумерацию с 1
    df_boosters.index.name = "Booster"
    
    edited_boosters = st.data_editor(df_boosters, use_container_width=True, key="ed_boosters")

    # Сохраняем таблицу, если она изменилась
    if st.session_state.ed_boosters:
        update_field(boosters_real_key, edited_boosters.to_dict('records'))