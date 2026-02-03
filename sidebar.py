import streamlit as st
from Algorithm.main import run_simulation
from data_manager import load_data

def sidebar(page_type):
    """Return selected options from the sidebar based on page type."""

    if page_type == "input_page":
        st.sidebar.subheader("Calculation")
        if st.sidebar.button("CALCULATE", type="primary", use_container_width=True):

            data = load_data()
            results = run_simulation(data)
            st.session_state["simulation_results"] = results
            st.switch_page("pages/2_📈_Results.py")
    elif page_type == "results_page":
        st.sidebar.subheader("Visualization Settings")
        plot_option = st.sidebar.radio(
                "Choose plot:",
                ["3D Orbit", "Velocity vs Time", "Altitude vs Time", "Mass vs Time", "Density vs Altitude", "Temperature profile", "Drag coefficient vs Mach"] 
            )
        return plot_option