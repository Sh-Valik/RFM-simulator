import streamlit as st

st.title("📚 Theory of Multi-Stage Rocket Flight Mechanics")

st.markdown("""
This page documents the complete mathematical and physical theory underlying the rocket trajectory simulation algorithm.
All equations are sourced directly from the *Rocket Flight Mechanics* lecture notes (Stoil Ivanov, 2025).
""")

# ── TABLE OF CONTENTS ──────────────────────────────────────────────────────────
with st.expander("📋 Table of Contents", expanded=False):
    st.markdown("""
    1. [Propulsion Fundamentals](#1-propulsion-fundamentals)
    2. [Multi-Stage Rocket Theory](#2-multi-stage-rocket-theory)
    3. [Boosters](#3-boosters)
    4. [Flight Profile & Trajectory](#4-flight-profile-trajectory)
    5. [Aerodynamic Forces](#5-aerodynamic-forces)
    6. [Gravity Modeling](#6-gravity-modeling)
    7. [Reference Frames & Coordinate Transformations](#7-reference-frames-coordinate-transformations)
    8. [Orbital Elements from Cartesian State](#8-orbital-elements-from-cartesian-state)
    """)

st.divider()

# ══════════════════════════════════════════════════════════════════════════════
# 1. PROPULSION FUNDAMENTALS
# ══════════════════════════════════════════════════════════════════════════════
st.header("1. Propulsion Fundamentals")

st.subheader("1.1 Rocket Equation of Motion")
st.markdown("""
The fundamental principle of rocketry is the generation of propulsive force through expulsion of mass.
Starting from Newton's 2nd law applied to a system with changing mass, and considering that the rocket
expels mass $-dm$ with relative velocity $V_{ex}$ (effective exhaust velocity), the equation of motion becomes:
""")
st.latex(r"""
m \frac{d\mathbf{V}}{dt} = \dot{m}\mathbf{V}_{ex} + \mathbf{F}_u
""")
st.markdown("""
where $\\dot{m} = -dm/dt > 0$ is the propellant mass flow rate and $\\mathbf{F}_u$ is the sum of all external forces.
This is known as the **Principle of Solidification** — Newton's 2nd law for constant-mass bodies can be applied
to rocket motion by adding the thrust term $\\dot{m}V_e$ to the external force sum.
""")

st.subheader("1.2 Thrust")
st.markdown("""
The thrust $T$ of a rocket engine consists of two components: the **impulse thrust** and the **pressure thrust**:
""")
st.latex(r"T = \underbrace{\dot{m} V_{ex}}_{\text{impulse thrust}} + \underbrace{(p_e - p_0)A_e}_{\text{pressure thrust}} = \dot{m} V_e")
st.markdown("""
where $V_e$ is the **effective exhaust velocity** that absorbs both terms, $p_e$ is the nozzle exit pressure,
$p_0$ is the ambient pressure, and $A_e$ is the nozzle exit area.
Note that thrust **increases with altitude** because $p_0$ decreases — maximum thrust is achieved in vacuum:
""")
st.latex(r"T_{vac} = \dot{m} V_{ex} + p_e A_e")

st.subheader("1.3 Specific Impulse")
st.markdown("""
The **specific impulse** $I_{sp}$ (in seconds) is an engine-independent figure of merit:
""")
st.latex(r"I_{sp} = \frac{T}{\dot{m} g_0} = \frac{V_e}{g_0}")
st.markdown("where $g_0 = 9.80665 \\ \\mathrm{m/s^2}$ is the standard gravitational acceleration.")

st.subheader("1.4 Tsiolkovsky Equation")
st.markdown("""
In gravity-free vacuum, integrating the equation of motion from initial mass $m_0$ to final mass $m_f$ yields
the celebrated **Tsiolkovsky rocket equation**:
""")
st.latex(r"\Delta V = V_e \ln \frac{m_0}{m_f} = V_e \ln \Lambda")
st.markdown(r"""
where $\Lambda = m_0 / m_f$ is the **burn mass ratio**. The velocity increment $\Delta V$ is **independent of the 
burn program** (how propellant is expelled over time) — it depends only on the exhaust velocity and mass ratio.

The burn time for a **constant propellant mass flow** program is:
""")
st.latex(r"t_b = \frac{m_0 - m_f}{\dot{m}} = \frac{V_e}{g_0 \Psi_0} \left(1 - \frac{1}{\Lambda}\right)")
st.markdown(r"where $\Psi_0 = T / (m_0 g_0)$ is the initial **thrust load** (thrust-to-weight ratio).")

st.subheader("1.5 Mass Notation")
st.markdown(r"""
For a single-stage rocket the total mass is:
$$m_{wet} = m_{dry} + m_{prop}$$
The key mass ratios are defined as:

| Ratio | Symbol | Formula |
|---|---|---|
| Construction mass ratio | $\varepsilon$ | $m_c \, / \, (m_c + m_p)$ |
| Propellant mass ratio | $\varphi$ | $m_p \, / \, m_0 = (1-\varepsilon)(1-\lambda)$ |
| Payload mass ratio | $\lambda$ | $m_u \, / \, m_0$ |
| Burn mass ratio | $\Lambda$ | $m_0 \, / \, m_f = 1\,/\,(1-\varphi)$ |

From these definitions it follows that:
$$\frac{1}{\Lambda} = 1 - \varphi = \lambda(1-\varepsilon) + \varepsilon$$
""")

st.divider()

# ══════════════════════════════════════════════════════════════════════════════
# 2. MULTI-STAGE ROCKETS
# ══════════════════════════════════════════════════════════════════════════════
st.header("2. Multi-Stage Rocket Theory")

st.subheader("2.1 Motivation")
st.markdown(r"""
Reaching orbit with a Single-Stage-to-Orbit (SSTO) launcher is practically impossible with current chemical propulsion.
The payload mass ratio decreases exponentially with required $\Delta V$:
$$\lambda = \frac{e^{-V_{f,id}/V_e} - \varepsilon}{1 - \varepsilon}$$
Typical losses (gravity + drag) require an additional $\Delta V / V_e \approx 1.0$, leaving insufficient margin for
a reasonable payload fraction. **Multi-staging** solves this by jettisoning empty hardware between burns.
""")

st.subheader("2.2 Multi-Stage Mass Ratios")
st.markdown(r"""
For a rocket with $N$ stages, subscript $i$ denotes the stage/section number. The payload of stage $i$ is the 
total mass of the section above it: $m_{u,i} = m_{0,i+1}$.

The mass ratios for each stage are defined identically to the single-stage case:
$$\varepsilon_i = \frac{m_{c,i}}{m_{c,i} + m_{p,i}}, \quad \varphi_i = \frac{m_{p,i}}{m_{0,i}}, \quad \lambda_i = \frac{m_{0,i+1}}{m_{0,i}}, \quad \Lambda_i = \frac{m_{0,i}}{m_{0,i} - m_{p,i}}$$

The **total payload mass ratio** is the product of all stage payload mass ratios:
$$\lambda_{tot} = \frac{m_u}{m_{0,1}} = \prod_{i=1}^{N} \lambda_i$$

The **total ideal delta-V** of an $N$-stage rocket is the sum of each stage's contribution:
$$V_{f,id} = \sum_{i=1}^{N} V_{e,i} \ln \Lambda_i = -\sum_{i=1}^{N} V_{e,i} \ln\bigl[\lambda_i(1-\varepsilon_i) + \varepsilon_i\bigr]$$
""")

tab_opt, tab_nonopt = st.tabs(["Optimal Rocket", "Non-Optimal Rocket"])

with tab_opt:
    st.subheader("Optimal Rocket (Maximum $\\Delta V$ for given $\\lambda_{tot}$)")
    st.markdown(r"""
    To maximise $V_{f,id}$ subject to $\prod \lambda_i = \lambda_{tot} = \text{const}$, the **Lagrange multiplier method** is applied.
    The Lagrangian is:
    $$\mathcal{L}(\lambda_i, \mu) = \sum_{i=1}^{N} V_{e,i} \ln\Lambda_i + \mu \left(\sum_{i=1}^{N} \ln\lambda_i - \ln\lambda_{tot}\right)$$
    Setting $\partial \mathcal{L}/\partial \lambda_i = 0$ for each stage yields the **optimal payload mass ratio**:
    $$\lambda_i^* = \frac{\mu \varepsilon_i}{(V_{e,i} - \mu)(1 - \varepsilon_i)}$$
    where $\mu$ is determined numerically from the constraint equation:
    $$\prod_{i=1}^{N} \frac{\mu \varepsilon_i}{(V_{e,i} - \mu)(1 - \varepsilon_i)} = \lambda_{tot}$$
    Taking the natural logarithm, this is equivalent to solving:
    $$\sum_{i=1}^{N} \left[\ln(\mu\varepsilon_i) - \ln(V_{e,i} - \mu) - \ln(1 - \varepsilon_i)\right] = \ln(\lambda_{tot})$$
    The algorithm solves this scalar equation for $\mu$ using the Brent root-finding method (bounded search: $\mu \in (0,\, \min_i V_{e,i})$).

    Once $\mu$ is known, the stage parameters follow:
    $$\varphi_i = (1-\varepsilon_i)(1-\lambda_i^*), \quad \Lambda_i = \frac{1}{1-\varphi_i}, \quad V_{f,id,i} = V_{e,i}\ln\Lambda_i$$

    The section masses are built bottom-up from the payload:
    $$m_{0,N} = \frac{m_u}{\lambda_N^*}, \quad m_{0,i} = \frac{m_{0,i+1}}{\lambda_i^*} \quad (i = N-1,\ldots,1)$$
    """)

with tab_nonopt:
    st.subheader("Non-Optimal Rocket ($\\lambda_i = \\text{const for all stages}$)")
    st.markdown(r"""
    The non-optimal (equal payload mass ratio) configuration sets:
    $$\lambda_i = \lambda_{tot}^{1/N} \quad \forall i$$
    All stage propellant mass ratios and burn mass ratios follow directly:
    $$\varphi_i = (1-\varepsilon_i)(1-\lambda_{tot}^{1/N}), \quad \Lambda_i = \frac{1}{1-\varphi_i}, \quad V_{f,id,i} = V_{e,i}\ln\Lambda_i$$
    The non-optimal rocket requires **approximately 12% more total initial mass** than the optimal configuration
    for the same payload and velocity increment.
    """)

st.divider()

# ══════════════════════════════════════════════════════════════════════════════
# 3. BOOSTERS
# ══════════════════════════════════════════════════════════════════════════════
st.header("3. Boosters")

st.markdown(r"""
Boosters burn in parallel with the **first stage**, splitting its burn time into two sub-phases:
$$0 \le t \le t_{b,s} \quad \text{(1st stage + boosters)} \qquad t_{b,s} < t \le t_{b,1} \quad \text{(1st stage alone)}$$

**Goal:** increase payload capacity $\Delta m_u$ while keeping the same total rocket structure (same upper stages).

The addition of boosters increases the effective first-section mass:
$$m_{0,1}' = m_{0,1} + \Delta m_u + m_{c,s} + m_{p,s}$$

The **corrected propellant mass ratio** of the first stage (normalised to the new section mass) is:
$$\varphi_1^c = \frac{m_{p,1}}{m_{0,1} + \Delta m_u} = \frac{\varphi_1}{1 + \Delta m_u / m_{0,1}}$$

The **burn mass ratio during the booster phase** is:
$$\Lambda_1' = \frac{1}{1 - \varphi_1^c \left(\frac{m_{p,s}}{m_{p,1}} + \frac{t_{b,s}}{t_{b,1}}\right) \Big/ \left(1 + \frac{\Delta m_u}{m_{0,1}} + \frac{\varphi_1^c m_{p,s}}{1-\varepsilon_s}\frac{1}{m_{p,1}}\right)}$$

Because the first-stage and booster exhaust velocities generally differ, an **equivalent exhaust velocity** 
is used during the booster phase:
$$V_{e,eq} = V_{e,1} - \frac{1}{1 + \dot{m}_1/\dot{m}_s}(V_{e,1} - V_{e,s})$$

The **burn mass ratio for the remainder of the first stage** (after booster jettison) is:
$$\Lambda_1'' = \frac{1 - \varphi_1^c \cdot (t_{b,s}/t_{b,1})}{1 - \varphi_1^c}$$

The **total ideal velocity increment with boosters** is:
$$V_{f,id}' = \underbrace{V_{e,eq} \ln \Lambda_1'}_{\text{booster phase}} + \underbrace{V_{e,1} \ln \Lambda_1''}_{\text{post-booster phase}} + \sum_{i=2}^{N} V_{e,i} \ln \Lambda_i$$
""")

st.divider()

# ══════════════════════════════════════════════════════════════════════════════
# 4. FLIGHT PROFILE & TRAJECTORY
# ══════════════════════════════════════════════════════════════════════════════
st.header("4. Flight Profile & Trajectory")

st.subheader("4.1 General Equations of Motion in 3D")
st.markdown(r"""
The full 3D equations of motion in an Earth-Centered Inertial (ECI) frame are integrated numerically.
The state vector is $\mathbf{s} = [x, y, z, \dot{x}, \dot{y}, \dot{z}, m]^\top$ and its time derivative is:
$$\dot{\mathbf{s}} = \begin{bmatrix} \dot{x} \\ \dot{y} \\ \dot{z} \\ \ddot{x} \\ \ddot{y} \\ \ddot{z} \\ \dot{m} \end{bmatrix} = \begin{bmatrix} V_x \\ V_y \\ V_z \\ (F_{g,x} + F_{D,x} + T_x)/m \\ (F_{g,y} + F_{D,y} + T_y)/m \\ (F_{g,z} + F_{D,z} + T_z)/m \\ -\dot{m}_{prop} \end{bmatrix}$$
The three force contributions — gravity, aerodynamic drag, and thrust — are described below.
""")

st.subheader("4.2 Flight Phases")
st.markdown(r"""
The simulation distinguishes the following flight phases for each stage:

| Phase | Description |
|---|---|
| **Vertical flight** | $0 \le t \le t_{vert}$. Thrust aligned with the local vertical (radial direction $\hat{r}$). |
| **Pitch program (below ~80–95 km)** | Thrust directed at fixed pitch angle $\theta$ from local vertical toward the azimuth. |
| **Gravity turn (above ~80–95 km)** | Thrust aligned with the velocity vector (zero angle of attack). |
| **Ballistic / coasting** | Engine off; only gravity and drag act on the rocket. |
| **Circularisation burn** | For the final stage, a fraction of propellant is reserved and burned at apogee. |

During the pitch program, the thrust vector is assembled in the local East-North-Up (ENU) frame and 
then expressed in ECI:
$$\hat{T}_{ECI} = \cos\theta \, \hat{r} + \sin\theta \left(\sin Az \, \hat{e} + \cos Az \, \hat{n}\right)$$
where $\hat{r}$, $\hat{e}$, $\hat{n}$ are the local up, east and north unit vectors respectively, and
$Az$ is the (corrected) launch azimuth angle.
""")

st.subheader("4.3 2D Analytical Summary (vacuum, flat Earth)")
st.markdown(r"""
For educational reference, the analytical solutions for the simplified 2D flat-Earth vacuum case are 
summarised below. The pitch angle $\theta$ is fixed and thrust is constant ($\dot{m} = \text{const}$).

**Velocity components (propelled flight):**
$$V_x(t) = V_e \ln\!\frac{m_0}{m} \cos\theta, \qquad V_z(t) = V_e \ln\!\frac{m_0}{m} \sin\theta - g_0 t$$

**Flight-path angle** $\gamma$ (angle between velocity vector and local horizon):
$$\tan\gamma = \tan\theta - \frac{g_0 t}{V_e \ln(m_0/m)\cos\theta}$$

**Burnout velocity magnitude:**
$$V_b = \sqrt{V_{f,id}^2 - 2V_{f,id} g_0 t_b \sin\theta + g_0^2 t_b^2}$$

**Apogee:** reached when $V_z = 0$:
$$t_{apo} = \frac{V_{f,id}}{g_0} \sin\theta, \qquad V_{apo} = V_{f,id}\cos\theta$$

Note that both $t_{apo}$ and $V_{apo}$ are **independent of the burn program**.
""")

st.subheader("4.4 Gravity Turn")
st.markdown(r"""
A gravity turn (zero-lift trajectory, $\alpha = 0$) keeps the angle of attack at zero by aligning the
thrust vector with the velocity vector at all times ($\gamma = \theta$). The equations of motion reduce to:
$$m \frac{dV}{dt} = T - mg_0 \sin\gamma \qquad V\frac{d\gamma}{dt} = -g_0 \cos\gamma$$
The second equation shows that **gravity steers the rocket**, eliminating aerodynamic loads.

In the numerical model, the gravity turn is activated above a threshold altitude (≈80–95 km) where the 
atmosphere is sufficiently rarefied. Below this altitude, a fixed pitch angle is maintained to establish 
the correct trajectory angle for the gravity turn to take over.
""")

st.divider()

# ══════════════════════════════════════════════════════════════════════════════
# 5. AERODYNAMIC FORCES
# ══════════════════════════════════════════════════════════════════════════════
st.header("5. Aerodynamic Forces")

tab_drag, tab_density, tab_mach, tab_dynpres = st.tabs(["Drag Force", "Atmospheric Density", "Mach Number & Speed of Sound", "Dynamic Pressure"])

with tab_drag:
    st.subheader("Aerodynamic Drag")
    st.markdown(r"""
    The drag force vector (opposing velocity) is:
    $$\mathbf{F}_D = -\frac{1}{2} \rho(h)\, C_D(M_a)\, A\, |\mathbf{V}|^2 \,\hat{V} = -\frac{1}{2}\rho(h)\, C_D(M_a)\, A \,|\mathbf{V}|\,\mathbf{V}$$
    where:
    - $\rho(h)$ — atmospheric density at altitude $h$
    - $C_D(M_a)$ — drag coefficient as a function of Mach number (tabulated data, linearly interpolated)
    - $A$ — cross-sectional area along the velocity vector
    - $\mathbf{V}$ — velocity vector of the rocket relative to the atmosphere

    **Cross-sectional area switching:** When the radial velocity component becomes negative  
    ($\dot{r} = \mathbf{r} \cdot \mathbf{V}/|\mathbf{r}| < 0$, i.e. the rocket is descending past apogee),  
    the falling stage is modelled as a cylinder falling **sideways**, so the frontal area switches from  
    $A_{pf}$ (along-flight) to $A_{bf}$ (broadside), and the drag coefficient switches to that of a cross-flow cylinder.

    Drag is set to zero when $\rho < 10^{-6}$ kg/m³ (effectively vacuum).
    """)

with tab_density:
    st.subheader("Atmospheric Density Model")
    st.markdown(r"""
    The atmospheric density decreases exponentially with altitude:
    $$\rho(h) = \rho_0 \exp\!\left(-\frac{h}{H_0}\right)$$
    where $\rho_0 = 1.293 \ \mathrm{kg/m^3}$ is the sea-level air density and $H_0$ is the **density scale height**:
    $$H_0 = -\frac{R \, T_0}{g_0 \, M_0 \, L} \approx 8432.56 \ \mathrm{m}$$
    with universal gas constant $R = 8.314 \ \mathrm{J/(mol \cdot K)}$, sea-level temperature $T_0 = 288.15 \ \mathrm{K}$,
    molar mass of air $M_0 = 0.02897 \ \mathrm{kg/mol}$, and temperature lapse rate $L = 0.0065 \ \mathrm{K/m}$.

    For the temperature profile used in computing the speed of sound, a tabulated atmospheric profile  
    (0–600 km altitude) is loaded and linearly interpolated at each integration step:
    $$T = T_{interp}(h)$$
    """)

with tab_mach:
    st.subheader("Mach Number & Local Speed of Sound")
    st.markdown(r"""
    The **local speed of sound** in air depends on the local temperature $T$:
    $$v_s = \sqrt{\gamma_a \, R_{\!f} \, T}$$
    where $\gamma_a = c_p/c_v \approx 1.4$ is the specific heat ratio of air and $R_f = 287.05 \ \mathrm{J/(kg \cdot K)}$ 
    is the specific gas constant of dry air.

    The **Mach number** is then:
    $$M_a = \frac{|\mathbf{V}|}{v_s}$$
    The drag coefficient $C_D(M_a)$ is read from a tabulated Mach–$C_D$ relation (loaded at startup) 
    and interpolated for the current Mach number at each integration step. This captures the 
    sharp transonic drag rise near $M_a \approx 1$.
    """)

with tab_dynpres:
    st.subheader("Dynamic Pressure")
    st.markdown(r"""
    Dynamic pressure represents the kinetic energy per unit volume of the fluid flowing around the rocket:
    $$q = \frac{1}{2}\rho(h)\,|\mathbf{V}|^2$$
    This quantity peaks at **MAX-Q** — the point of maximum aerodynamic structural loading during ascent.
    At launch, velocity is low despite high density; at high altitude, density collapses faster than velocity
    rises, so $q$ passes through a maximum and then decreases.

    Dynamic pressure is used internally to weight the drag force and is available as a simulation output.
    """)

st.divider()

# ══════════════════════════════════════════════════════════════════════════════
# 6. GRAVITY MODELING
# ══════════════════════════════════════════════════════════════════════════════
st.header("6. Gravity Modeling")

tab_j2, tab_point = st.tabs(["J2 Perturbation (used in simulation)", "Point-mass reference"])

with tab_j2:
    st.subheader("J2 Gravitational Acceleration")
    st.markdown(r"""
    Earth is not a perfect sphere — it bulges at the equator. The dominant non-spherical effect is 
    captured by the **J2 zonal harmonic** ($J_2 = 0.00108263$). The gravitational potential including J2 is:
    $$V_{J_2} = -\frac{GM}{r} - \frac{GM}{r}\frac{J_2 R_e^2}{r^2} P_2(\sin\varphi)$$
    where $P_2(\sin\varphi) = \tfrac{1}{2}(3\sin^2\varphi - 1)$ is the second Legendre polynomial and
    $\sin\varphi = z/r$ in Cartesian coordinates.

    Differentiating this potential, the J2 gravitational acceleration components are:
    """)
    st.latex(r"""
    g_x = -\frac{GM}{r^3}x + \frac{3 J_2 GM R_e^2}{2} \cdot \frac{x}{r^5}\left(1 - \frac{5z^2}{r^2}\right)
    """)
    st.latex(r"""
    g_y = -\frac{GM}{r^3}y + \frac{3 J_2 GM R_e^2}{2} \cdot \frac{y}{r^5}\left(1 - \frac{5z^2}{r^2}\right)
    """)
    st.latex(r"""
    g_z = -\frac{GM}{r^3}z + \frac{3 J_2 GM R_e^2}{2} \cdot \frac{z}{r^5}\left(3 - \frac{5z^2}{r^2}\right)
    """)
    st.markdown(r"""
    where $r = \sqrt{x^2 + y^2 + z^2}$, $G = 6.6742\times10^{-11}$ N·m²/kg² and $M = 5.97219\times10^{24}$ kg.
    The gravity force on the rocket is $\mathbf{F}_g = m\,[g_x,\, g_y,\, g_z]^\top$.

    Gravity is set to zero if the rocket descends below Earth's surface ($r < R_E$).
    """)

with tab_point:
    st.subheader("Point-Mass Gravity (simplified reference)")
    st.markdown(r"""
    The simplified point-mass gravitational model (used for analytical derivations) gives:
    $$\mathbf{g} = -\frac{GM}{r^2}\hat{r} = -\frac{GM}{r^3}\mathbf{r}$$
    At sea level this equals the standard acceleration $g_0 = GM/R_E^2 = 9.80665 \ \mathrm{m/s^2}$.

    For most analytical estimates, gravity is treated as constant ($g = g_0$) because typical flight altitudes 
    are much smaller than Earth's radius (at 200 km, $g \approx 0.94\, g_0$).
    """)

st.divider()

# ══════════════════════════════════════════════════════════════════════════════
# 7. REFERENCE FRAMES & TRANSFORMATIONS
# ══════════════════════════════════════════════════════════════════════════════
st.header("7. Reference Frames & Coordinate Transformations")

st.markdown(r"""
The simulation pipeline traverses the following reference frame chain:
$$\underbrace{\text{Geodetic (WGS84)}}_{\text{launch pad}} \xrightarrow{(1)} \underbrace{\text{ECI (Greenwich-based)}}_{\text{integration frame}} \xrightarrow{(2)} \underbrace{\text{J2000 Inertial}}_{\text{results}} \xrightarrow{(3)} \underbrace{\text{Keplerian Elements}}_{\text{orbit}}$$
""")

step1, step2, step3, step4 = st.tabs([
    "① Geodetic → Cartesian",
    "② Julian Date & GMST",
    "③ ECI → J2000",
    "④ Azimuth Correction"
])

with step1:
    st.subheader("Geodetic (WGS84) → Cartesian ECI")
    st.markdown(r"""
    The **World Geodetic System 1984 (WGS84)** ellipsoid describes Earth's shape with semi-major axis
    $a = 6378137$ m and first eccentricity squared $e^2 = 0.00669437999014$.

    The **radius of the prime vertical** at geodetic latitude $La$ is:
    $$N = \frac{a}{\sqrt{1 - e^2 \sin^2(La)}}$$

    The Cartesian ECEF coordinates of the launch pad at geodetic $(La, Lon, h)$ are:
    $$x = (N + h)\cos(La)\cos(Lon)$$
    $$y = (N + h)\cos(La)\sin(Lon)$$
    $$z = \bigl[N(1-e^2) + h\bigr]\sin(La)$$

    These Cartesian coordinates become the **initial position** $\mathbf{r}_0$ of the rocket.

    **Initial velocity** from Earth's rotation: The launch pad moves eastward with speed
    $V_{pad} = V_{eq}\cos(La)$ where $V_{eq} = 2\pi R_E / t_{day}$.
    In ECI coordinates this contributes a velocity perpendicular to both the spin axis and the position vector.
    """)

with step2:
    st.subheader("UTC → Julian Date → GMST")
    st.markdown(r"""
    To align the ECI Greenwich-based frame with the inertial J2000 frame, the **Greenwich Mean Sidereal Time (GMST)**
    at the moment of launch must be known.

    **Step 1 — Julian Date from UTC:**
    $$D_f = D + \frac{h}{24} + \frac{m}{1440} + \frac{s}{86400}$$
    If month $M \le 2$: set $Y' = Y-1$, $M' = M+12$, else $Y'=Y$, $M'=M$.
    $$A = \left\lfloor\frac{Y'}{100}\right\rfloor, \quad B = 2 - A + \left\lfloor\frac{A}{4}\right\rfloor$$
    $$J_D = \lfloor 365.25(Y'+4716)\rfloor + \lfloor 30.6001(M'+1)\rfloor + D_f + B - 1524.5$$

    **Step 2 — Centuries since J2000:**
    $$T = \frac{J_D - 2451545.0}{36525}$$

    **Step 3 — GMST (in degrees, reduced to $[0°,360°)$):**
    $$\theta_{GMST} = 280.46061837 + 360.98564736629\,(J_D - 2451545.0) + 0.000387933\,T^2 - \frac{T^3}{38710000}$$
    """)

with step3:
    st.subheader("ECI (Greenwich-based) → J2000 Cartesian")
    st.markdown(r"""
    The angle $\theta_{GMST}$ is the rotation from the J2000 x-axis (vernal equinox) to the Greenwich meridian, 
    measured about Earth's spin axis (z-axis, common to both frames).

    The rotation matrix from J2000 to ECI Greenwich is $R_z(\theta_{GMST})$; therefore the **inverse rotation** 
    $R_z(-\theta_{GMST})$ transforms from ECI to J2000:
    $$\mathbf{r}_{J2000} = R_z(-\theta_{GMST})\,\mathbf{r}_{ECI}$$
    $$\mathbf{V}_{J2000} = R_z(-\theta_{GMST})\,\mathbf{V}_{ECI}$$
    where the rotation matrix is:
    $$R_z(-\theta) = \begin{pmatrix}\cos\theta & \sin\theta & 0 \\ -\sin\theta & \cos\theta & 0 \\ 0 & 0 & 1\end{pmatrix}$$
    Component-wise:
    $$x_{J} = x\cos\theta + y\sin\theta, \quad y_{J} = -x\sin\theta + y\cos\theta, \quad z_{J} = z$$
    and identically for velocity components.
    """)

with step4:
    st.subheader("Azimuth Correction for Earth's Rotation")
    st.markdown(r"""
    Without correction, the naive azimuth $Az_0 = \arcsin(\cos i / \cos La)$ (from the inclination $i$ and 
    launch latitude $La$) ignores Earth's rotation. The **corrected azimuth** accounts for the fact that the 
    launch pad already contributes an eastward velocity $V_{pad} = V_{eq}\cos(La)$.

    For a circular target orbit the procedure is:

    1. Compute orbital speed: $V_{orb} = \sqrt{\mu/(R_E + h_{orb})}$
    2. Break orbital velocity into North/East components using the **uncorrected** azimuth:
       $$V_{orb,N} = V_{orb}\cos(Az_0), \qquad V_{orb,E} = V_{orb}\sin(Az_0)$$
    3. Subtract the launch pad velocity from the eastward component to get the **required launch velocity**:
       $$V_{L,N} = V_{orb,N}, \qquad V_{L,E} = V_{orb,E} - V_{pad}$$
    4. The **corrected azimuth** is then:
       $$Az_{corr} = \arctan_2(V_{L,E},\, V_{L,N})$$
    5. The required launch velocity magnitude:
       $$V_L = \sqrt{V_{L,N}^2 + V_{L,E}^2}$$
    6. The velocity gain from Earth's rotation:
       $$\Delta V_{rot} = V_{pad}\sin(Az_{corr}) = V_{eq}\cos(La)\sin(Az_{corr})$$
    """)

st.divider()

# ══════════════════════════════════════════════════════════════════════════════
# 8. ORBITAL ELEMENTS FROM CARTESIAN STATE
# ══════════════════════════════════════════════════════════════════════════════
st.header("8. Orbital Elements from Cartesian State")

st.markdown(r"""
Given the J2000 Cartesian state $(\mathbf{r},\mathbf{V})$ at the end of powered flight, the **Keplerian orbital elements**
are derived as follows. Let $\mu = GM$.
""")

col1, col2 = st.columns(2)

with col1:
    st.markdown("**Specific angular momentum:**")
    st.latex(r"\mathbf{h} = \mathbf{r} \times \mathbf{V}")

    st.markdown("**Node vector:**")
    st.latex(r"\mathbf{n} = \hat{k} \times \mathbf{h}, \quad \hat{k} = (0,0,1)^\top")

    st.markdown("**Eccentricity vector:**")
    st.latex(r"\mathbf{e} = \frac{1}{\mu}\left(\mathbf{V}\times\mathbf{h} - \mu\frac{\mathbf{r}}{r}\right)")

    st.markdown("**Semi-major axis (vis-viva):**")
    st.latex(r"a = \left(\frac{2}{r} - \frac{V^2}{\mu}\right)^{-1}")

with col2:
    st.markdown("**Inclination:**")
    st.latex(r"i = \arccos\!\left(\frac{h_z}{h}\right)")

    st.markdown("**Right Ascension of Ascending Node (RAAN):**")
    st.latex(r"\Omega = \begin{cases}\arccos(n_x/n) & n_y \ge 0 \\ 2\pi - \arccos(n_x/n) & n_y < 0\end{cases}")

    st.markdown("**Argument of Perigee:**")
    st.latex(r"\omega = \begin{cases}\arccos\!\left(\frac{\mathbf{n}\cdot\mathbf{e}}{ne}\right) & e_z \ge 0 \\ 2\pi - \arccos\!\left(\frac{\mathbf{n}\cdot\mathbf{e}}{ne}\right) & e_z < 0\end{cases}")

    st.markdown("**True Anomaly:**")
    st.latex(r"\nu = \begin{cases}\arccos\!\left(\frac{\mathbf{e}\cdot\mathbf{r}}{er}\right) & \mathbf{r}\cdot\mathbf{V} \ge 0 \\ 2\pi - \arccos\!\left(\frac{\mathbf{e}\cdot\mathbf{r}}{er}\right) & \mathbf{r}\cdot\mathbf{V} < 0\end{cases}")

st.divider()

# ── VELOCITY BUDGET SUMMARY ────────────────────────────────────────────────────
st.header("9. Velocity Budget")

st.markdown(r"""
The total required delta-V to reach a target orbit is:
$$\Delta V_{req} = \Delta V_{me} + \Delta V_{drag} + \Delta V_{grav} - \Delta V_{rot}$$

| Component | Equation | Description |
|---|---|---|
| **Mechanical energy** | $\Delta V_{me} = \sqrt{\frac{2\mu}{R_E+h_{pad}} - \frac{\mu}{R_E+h_{orb}}}$ | $\Delta V$ to change orbital energy from pad to orbit |
| **Drag loss** | $\Delta V_{drag} = \int_0^{t_{drag}} a_{drag}\,dt$ | Velocity lost fighting atmospheric drag |
| **Gravity loss** | $\Delta V_{grav} = \int_0^{t_{flight}} g(h)\sin\gamma(t)\,dt$ | Velocity lost fighting gravity in radial direction |
| **Earth rotation gain** | $\Delta V_{rot} = V_{eq}\cos(La)\sin(Az)$ | Free velocity boost from Earth's rotation |

Typical values for LEO launch: $\Delta V_{drag} \approx 100$ m/s, $\Delta V_{grav} \approx 1.5$–$2.2$ km/s.
""")

st.divider()

st.caption("""
All theory derived from: Stoil Ivanov, *Rocket Flight Mechanics*, Lecture Notes 1–11, Sofia, Bulgaria, 2025.
Reference sources: Wittenburg (TU Delft, 2014), Miele (1962), Bate–Mueller–White (1971), US Standard Atmosphere (1976).
""")