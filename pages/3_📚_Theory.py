import streamlit as st

st.title("📚 Theory of RFM Analysis")
st.header("Workflow Overview")

# placeholder for workflow image


st.header("1. Data Collection")
st.write("""
The first step in RFM analysis is to gather data from users. This data typically includes:
""")
         
col1, col2, col3 = st.columns(3)
with col1:
    st.subheader("Rocket Parameters:")
    st.markdown("""
    - Payload mass
    - Theta angle
    - Number of stages
    - Number of boosters if applicable
    - Effective exhaust velocity of each stage (and boosters if applicable)
    - Mass flow rate of each stage (and boosters if applicable)
    - Start mass and propellant mass of each stage (and boosters if applicable) or construction mass ratio of each stage and payload mass ratio of entire rocket
    - Type of a rocket:
            - Optimal
            - Non-Optimal
    """)

with col2:
    st.subheader("Launch Parameters")
    st.markdown("""
    - Launch latitude
    - Launch longitude
    - Launch altitude
    - Launch date and time
    """)

with col3:
    st.subheader("Target Orbit Parameters")
    st.markdown("""
    - Semi-major axis
    - Eccentricity
    - Inclination
    - Right ascension of ascending node
    - Argument of perigee
    - True anomaly
    """)


st.header("2. Preprocessing")
st.write("""
         
         """)

























# tab1, tab2, tab3 = st.tabs(["🚀 Propulsion", "🌍 Trajectory", "📊 Optimization"])

# with tab1:
#     st.header("Propulsion & Tsiolkovsky Equation")
#     st.write("""
#     The core of rocket motion is based on momentum conservation. 
#     The change in velocity is calculated using the Tsiolkovsky equation:
#     """)
    
#     # Математика в LaTeX — мастхэв для курсовой
#     st.latex(r"\Delta v = v_e \ln \left( \frac{m_0}{m_f} \right)")
    
#     with st.expander("See implementation details"):
#         st.code("""
# def calculate_delta_v(v_e, m_start, m_final):
#     return v_e * np.log(m_start / m_final)
#         """, language="python")