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
         Before proceeding with the simulation based on user-defined parameters, it is essential to establish the mathematical and physical foundations and core functions that drive the algorithm.
         """)

tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs(["Optimal Rocket", "Non-optimal Rocket", "Aeorodynamic Force", "Gravity Force", "Propulsion", "Derivatives"])

with tab1:
    st.subheader('Optimal Rocket')
    st.write("""
    When user select "Optimal" type of a rocket, the algorithm will calculate the optimal rocket parameters based on the user-defined total payload mass ratio and construction mass ratio.
    As we know
    $$
    \lambda_{total} = \lambda_1 \cdot \lambda_2 \cdot \lambda_3 \cdots \lambda_N
    $$
    where $\lambda_{total}$ is the total payload mass ratio, $\lambda_1, \lambda_2, \lambda_3, \cdots, \lambda_N$ are the payload mass ratios of each stage.
    """)
    

with tab2:
    st.subheader("Non-optimal Rocket")
    st.write("""
    
    """)

with tab3:
    st.subheader('Aeorodynamic Force')

with tab4:
    st.subheader('Gravity Force')

with tab5:
    st.subheader('Propulsion')

with tab6:
    st.subheader('Derivatives')


st.header("3. Simulation")





















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