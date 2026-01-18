import streamlit as st

def sidebar(page_type):

    if page_type == "input_page":
        st.sidebar.subheader("Calculation")
        if st.sidebar.button("CALCULATE", type="primary", use_container_width=True): 
            st.switch_page("pages/2_📈_Results.py")
    elif page_type == "results_page":
        st.sidebar.subheader("Visualization Settings")
        plot_option = st.sidebar.radio(
                "Choose plot:",
                ["3D Orbit", "Velocity Profile", "Altitude"]
            )
        return plot_option