# 🚀 RFM Simulator - Rocket Flight Mechanics Simulator

A comprehensive, interactive Python application for simulating and analyzing rocket trajectories, rocket simulation, and flight dynamics. This simulator uses advanced numerical methods to model multi-stage rockets with optional payload mass distribution, accounting for atmospheric drag, gravitational forces, rotation of the Earth, and thrust dynamics.

---

## 📋 Table of Contents

1. [Overview](#overview)
2. [Features](#features)
3. [Project Structure](#project-structure)
4. [Installation & Setup](#installation--setup)
5. [System Requirements](#system-requirements)
6. [Usage Guide](#usage-guide)
7. [Core Concepts](#core-concepts)
8. [Module Documentation](#module-documentation)
9. [Physics & Mathematics](#physics--mathematics)
10. [Configuration Guide](#configuration-guide)
11. [Visualization & Results](#visualization--results)
12. [Advanced Features](#advanced-features)
13. [Troubleshooting](#troubleshooting)
14. [Contributing](#contributing)

---

## Overview

The RFM (Rocket Flight Mechanics) Simulator is an advanced educational and research tool designed to simulate rocket launches and predict orbital insertion trajectories. The simulator integrates real-world physical models including:

- **Multi-stage rocket propulsion systems** with configurable mass flow rates and thruster parameters
- **Booster systems** that detached during flight
- **Atmospheric drag** using Mach-dependent drag coefficient data
- **Gravitational modeling** based on Earth's actual mass and radius
- **Atmospheric density** and temperature profiles based on real atmospheric data
- **Orbital mechanics** with Keplerian elements for target orbit specification

The application features a user-friendly web interface built with **Streamlit**, allowing users to:
- Configure launch parameters and rocket specifications
- Define target orbital characteristics
- Run simulations and visualize results
- Analyze flight performance across multiple metrics
- Read a fairly detailed theory that the algorithm uses

---

## Features

### 🎯 Core Simulation Features

- **Multi-Stage Rocket Support**: Model rockets with more than one stage with independent thrust profiles
- **Booster Systems**: Boosters that can be configured separately
- **Dual Input Modes**:
  - Direct mass input (Start mass & Propellant mass)
  - Construction and payload mass ratio input (EPS & Payload mass ratio)
- **Flexible Rocket Types**: Support for different rocket configuration models
- **Real Atmospheric Data**: Temperature and density profiles from sea level to 600 km altitude
- **Drag Coefficient Database**: Mach-dependent drag coefficients for accurate aerodynamic modeling
- **Time-Dependent Physics**: 
  - Atmospheric density decay
  - Temperature profile changes
  - Mach number-dependent drag

### 🌍 Launch & Orbital Mechanics

- **Geographic Launch Parameters**:
  - Latitude/Longitude specification
  - Launch altitude above sea level
  - Custom launch date and time
- **Keplerian Orbital Elements**:
  - Semi-major axis (a)
  - Eccentricity (e)
  - Inclination (i)
  - Right Ascension of Ascending Node (Ω)
  - Argument of Periapsis (ω)
  - True Anomaly (ν)
- **3D Orbit Visualization**: Real-time visualization of target orbits in 3D space
- **Rocket Position Tracking**: Display rocket position relative to target orbit

### 📊 Visualization & Analysis

- **3D Trajectory Plotting**: Complete rocket flight path in 3D space
- **Velocity Analysis**: Velocity vs. time throughout flight phases
- **Altitude Tracking**: Altitude changes during ascent and coasting phases
- **Mass Depletion**: Rocket mass reduction over time due to propellant consumption
- **Atmospheric Profiles**: 
  - Air density vs. altitude
  - Temperature profile visualization
  - Mach-dependent drag coefficient curves
- **Interactive Plots**: Fully interactive Plotly-based visualizations

### ⚙️ Advanced Configuration

- **Flexible Input Formats**:
  - Start mass and propellant mass per stage
  - Energy ratios (EPS) with mass flow rates
- **Payload Management**: Configurable payload mass and mass ratios
- **Propulsion Parameters**:
  - Exhaust velocity (Ve) per stage
  - Mass flow rates
  - Thrust vector control angles
- **Session Management**: Automatic temporary file handling for multi-session workflows

---

## Project Structure

```
RFM-simulator/
├── 🚀 Rocket_parameters.py          # Main entry point and rocket configuration
├── data_manager.py                  # Session data persistence and management
├── sidebar.py                       # Streamlit UI sidebar navigation
├── requirements.txt                 # Python dependencies
├── Algorithm/
│   ├── __init__.py                 # Package initialization
│   ├── main.py                     # Main simulation execution logic
│   └── functions.py                # Core physics and numerical functions
├── pages/
│   ├── 1_🌍_Launche_&_Target_orbit.py   # Launch parameter input interface
│   ├── 2_📈_Results.py               # Results visualization and analysis
│   └── 3_📚_Theory.py               # Educational theory and reference material
├── resources/
│   ├── temperature_profile_0_to_600km.txt  # Atmospheric temperature data
│   └── CD-Mach_relation.txt               # Drag coefficient vs. Mach data
└── temp_sessions/                   # Temporary session storage (auto-generated)
```

### File Descriptions

#### **🚀 Rocket_parameters.py**
- **Purpose**: Main configuration and initialization file
- **Contains**:
  - Physical constants (gravitational constant, Earth radius, Earth mass, etc.)
  - Rocket and mission parameters loading
  - Stage and booster configuration processing
  - Launch location and orbital target parameters
  - Initial state setup for simulation

#### **data_manager.py**
- **Purpose**: Manage application state and user input persistence
- **Key Features**:
  - Session-based data storage (JSON format)
  - Automatic cleanup of old session files
  - Default configuration templates
  - Real-time data updates from UI input
  - Multi-user support through session isolation

#### **Algorithm/functions.py**
- **Purpose**: Core physics calculations and numerical methods
- **Key Functions**:
  - `gravity()`: Gravitational acceleration computation
  - `density()`: Air density at given altitude
  - `temperature_by_altitude()`: Temperature interpolation from data
  - `drag_interp()`: Drag coefficient interpolation
  - `propulsion()`: Thrust and mass flow calculation
  - `Derivatives()`: ODE system derivatives computation
  - `parameters_of_stages()`: Stage parameter processing

#### **Algorithm/main.py**
- **Purpose**: Main simulation execution
- **Contains**:
  - Numerical integration setup
  - Initial condition definition
  - ODE solver execution
  - Results post-processing

#### **pages/1_🌍_Launche_&_Target_orbit.py**
- **Purpose**: User interface for input parameters
- **Features**:
  - Launch parameter form (location, date, time)
  - Target orbit Keplerian element input
  - Interactive 3D orbit visualization
  - Real-time validation

#### **pages/2_📈_Results.py**
- **Purpose**: Results visualization and analysis
- **Features**:
  - Multiple visualization options
  - Interactive plot selection
  - Data analysis tools
  - Export capabilities

#### **pages/3_📚_Theory.py**
- **Purpose**: Educational reference material
- **Contains**: Physics background, equations, and explanations

#### **resources/**
- **temperature_profile_0_to_600km.txt**: Altitude-temperature lookup table
- **CD-Mach_relation.txt**: Mach number-drag coefficient lookup table

---

## Version History

- **v1.0.0** (Current): Initial release with multi-stage rocket simulation, atmospheric modeling, and orbital mechanics

---

### Common Abbreviations

| Abbrev | Definition |
|--------|-----------|
| LEO | Low Earth Orbit |
| MEO | Medium Earth Orbit |
| GEO | Geostationary Orbit |
| ECI | Earth-Centered Inertial |
| RAAN | Right Ascension of Ascending Node |
| Ve | Exhaust Velocity |
| Isp | Specific Impulse |
| Cd | Drag Coefficient |
| ODE | Ordinary Differential Equation |
| SRB | Solid Rocket Booster |
| STS | Space Transportation System |

---

**Last Updated**: January 2026
**Maintainer**: Development Team
**Status**: Active Development

