AirDensityDrive: F1 Power & Top Speed Analysis

A MATLAB-based simulation project exploring how air density affects the engine power output, drag force, and theoretical top speed of Formula 1 cars across all 24 tracks on the 2024 calendar. The model compares regulatory changes between 2025 and 2026 power unit configurations, including the switch to 100% sustainable fuel.

## 📌 Objectives

- Model the effect of ambient air density (influenced by elevation and temperature) on engine performance and vehicle top speed.
- Quantify power output and maximum velocity under:
  - 2025 regulations
  - 2026 regulations with increased MGU-K contribution
  - 2026 using ethanol-based sustainable fuel
## ⚙️ Methodology

1. **Track Data Collection**
   - Elevation and average temperature data for each circuit.
2. **Air Density Estimation**
   - Barometric pressure formula used with ideal gas law.
3. **Engine Model**
   - Choked flow assumption for maximum intake airflow.
   - Combustion power calculated using LHV and efficiency.
4. **Top Speed Solver**
   - Power balance solved numerically using `fzero`.
5. **Validation**
   - Comparison with publicly available FIA-reported top speeds (limited to available circuits).

---

## 📊 Key Results

- Overall model accuracy: **~95.5%** compared to known FIA top speeds.
- Demonstrated that **high-altitude circuits** (e.g., Mexico City, São Paulo) show reduced ICE power but still high top speeds due to lower aerodynamic drag.
- 2026 hybrid upgrades (MGU-K power boost from 120 kW to 350 kW) show **clear top-speed gains**, especially at low-density tracks.
- Ethanol fuel results in **only minor ICE power drop (~4 kW)** due to lower LHV but allows more fuel flow due to lower AFR.
