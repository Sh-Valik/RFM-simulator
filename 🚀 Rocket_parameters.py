import streamlit as st
import pandas as pd
from sidebar import sidebar
from data_manager import load_data, update_field

st.set_page_config(layout="wide")
sidebar("input_page")

data = load_data()

def counter_sync(session_state_of_count, count_data):
    """Синхронизирует количество строк в таблице с числом степеней или бустеров"""
    new_count = st.session_state[session_state_of_count]
    current_data = load_data()
    current_list = current_data.get(count_data, [])

    if len(current_list) < new_count:
        # Добавляем пустые строки
        for _ in range(new_count - len(current_list)):
            current_list.append({"Propellant (kg)": 1000.0, "Mass flow (kg/s)": 300.0, "Ve (m/s)": 300.0, "Structural Mass (kg)": 100.0})
    elif len(current_list) > new_count:
        # Удаляем лишние строки
        current_list = current_list[:new_count]

    current_data[session_state_of_count] = new_count
    current_data[count_data] = current_list
    update_field(session_state_of_count, new_count)
    update_field(count_data, current_list)
    

######################
st.title("🚀 Rocket Parameters")
st.header("General Parameters")
col1, col2, col3 = st.columns(3)

with col1:
    st.number_input("Payload Mass (kg)", min_value=0.0, value=data["payload_mass"], step=100.0, key="payload_mass", on_change=lambda: update_field("payload_mass", st.session_state.payload_mass))
    st.checkbox("Has Boosters", value=data["has_boosters"], key="has_boosters", on_change=lambda: update_field("has_boosters", st.session_state.has_boosters))

with col2:
    st.number_input("Stages Count", min_value=1, max_value=10, value=data["stages_count"], step=1, key="stages_count", on_change=lambda: counter_sync("stages_count", "stages_data"))
with col3:
    if data["has_boosters"]:
        st.number_input("Booster Count", min_value=1, max_value=10, value=data["booster_count"], step=1, key="booster_count", on_change=lambda: counter_sync("booster_count", "boosters_data"))

st.divider()
st.header("Stages Details")
df_stages = pd.DataFrame(data["stages_data"])
df_stages.index += 1  # Начинаем нумерацию с 1
df_stages.index.name = "Stage"
edited_stages = st.data_editor(df_stages, use_container_width=True, key="ed_stages")

# Сохраняем таблицу, если она изменилась
if st.session_state.ed_stages:
    update_field("stages_data", edited_stages.to_dict('records'))

if data["has_boosters"]:
    st.divider()
    st.header("Boosters Details")
    df_boosters = pd.DataFrame(data["boosters_data"])
    df_boosters.index += 1  # Начинаем нумерацию с 1
    df_boosters.index.name = "Booster"
    edited_boosters = st.data_editor(df_boosters, use_container_width=True, key="ed_boosters")

    # Сохраняем таблицу, если она изменилась
    if st.session_state.ed_boosters:
        update_field("boosters_data", edited_boosters.to_dict('records'))