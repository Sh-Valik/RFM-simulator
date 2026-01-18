import streamlit as st
import pandas as pd
from sidebar import sidebar

st.set_page_config(layout="wide")

st.title("🚀 Rocket Parameters")

sidebar("input_page")
st.header("General Parameters")

col1, col2, col3 = st.columns(3)
with col1:
    payload_mass = st.number_input("Payload Mass (kg)", min_value=0.0, value=1000.0, step=100.0)
    has_boosters = st.checkbox("Has Boosters")
    booster_count = 0   
with col2:
    stages_count = st.number_input("Stages Count", min_value=1, max_value=10, value=2, step=1)
with col3:
    if has_boosters:
        booster_count = st.number_input("Booster Count", min_value=0, max_value=10, value=0, step=1)



st.divider()
st.header("Stages Details")
default_data = {
    "Propellant (kg)": [1000.0 for _ in range(stages_count)],
    "Mass flow (kg/s)": [300.0 for _ in range(stages_count)],
    "Ve (m/s)": [300.0 for _ in range(stages_count)],
    "Structural Mass (kg)": [100.0 for _ in range(stages_count)],
}
df_stages = pd.DataFrame(default_data)
df_stages.index.name = "Stage"
df_stages.index += 1  # Start index from 1 for stages
edited_df = st.data_editor(
    df_stages,
    use_container_width=True,
    num_rows="fixed",
    hide_index=False
)

if has_boosters:
    st.divider()
    st.header("Boosters Details")
    default_booster_data = {
        "Propellant (kg)": [500.0 for _ in range(booster_count)],
        "Mass flow (kg/s)": [200.0 for _ in range(booster_count)],
        "Ve (m/s)": [250.0 for _ in range(booster_count)],
        "Structural Mass (kg)": [50.0 for _ in range(booster_count)],
    }
    df_boosters = pd.DataFrame(default_booster_data)
    df_boosters.index.name = "Booster"
    df_boosters.index += 1  # Start index from 1 for boosters
    edited_booster_df = st.data_editor(
        df_boosters,
        use_container_width=True,
        num_rows="fixed",
        hide_index=False
    )