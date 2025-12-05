# Marcos Esparza – Engineering Portfolio  
Mechanical Engineering • Energy & Thermo-Fluids • High-Power Rocketry  

[About](#about-me) • [Resume](#resume) • [Featured Project](#featured-project--irec-2026) • [Projects](#projects) • [Skills](#skills) • [Contact](#contact)

---

## About Me

I’m a Mechanical Engineering student at the University of Texas Permian Basin with a focus on energy systems, thermo-fluids, and experimental testing. I work at the intersection of high-power rocketry and lab-scale energy hardware, where I care about not just modeling performance but also how hardware is built, instrumented, and verified.

I am part of the Falcon Aeronautics & Space Team (FAST), where I lead aerodynamic design and fin-can development for our 10,000-ft COTS rocket competing at the Intercollegiate Rocket Engineering Competition (IREC). My work spans fin geometry, flutter and static-stability analysis, and integrated airframe modeling in SOLIDWORKS, MATLAB, and OpenRocket.

In the lab, I have run experiments on cross-flow heat exchangers, fluid friction and bend losses, flow-measurement devices, transient cooling, and velocity profiles. These projects have given me hands-on experience with piping systems, pressure drop, heat transfer, and DAQ/LabVIEW-style measurement setups—the same fundamentals that show up in real energy and pipeline infrastructure.

*Currently seeking Summer 2026 internships in energy systems, turbomachinery, or mechanical design and testing.*

---

## Resume

- **[📄 View Resume (PDF)](assets/img/MarcosEsparza_Resume.pdf)**
- **[🔗 View IREC GitHub Repo](https://github.com/MarcosEsparza/Falcon-Aeronautics-and-Space-Team)**

---

## Featured Project — IREC 2026  

**Role:** Aerodynamic Design Lead  
**Team:** Falcon Aeronautics & Space Team (FAST), UTPB  
**Tools:** SOLIDWORKS, MATLAB, OpenRocket, Excel  

### Overview

- Lead aerodynamic design and fin-can development for FAST’s 10,000-ft COTS rocket for IREC.  
- Design a modular fin-can and airframe to safely reach ~10,000 ft on an L-class motor while maintaining ≥1.5× flutter margin and competition-typical stability ranges.  
- Treat the rocket as a small field system: define loads and margins, keep interfaces clean in CAD, and document decisions so the hardware is safe, repeatable, and buildable by the team.

### Technical Scope

- Designed fin geometry and performed flutter and static-stability analysis using MATLAB and OpenRocket.  
- Modeled a 4.02-in modular fiberglass airframe with an aluminum fin-can in SOLIDWORKS, including internal structure.  
- Verified ≥1.5× V_max flutter margin and maintained stability margin within competition-typical ranges.  
- Integrated internal structure, recovery system, and avionics into a unified CAD model to align mass properties and mechanical interfaces.  
- Established a documentation workflow for the IREC 2026 technical report and team GitHub repository.  
- **MATLAB Code:** [Fin Flutter Solver](https://github.com/MarcosEsparza/Falcon-Aeronautics-and-Space-Team/blob/main/simulations/matlab/fin_flutter/AluminumFinCan.m)  

### Flight Simulation — L1500T

- Apogee: **10,990 ft**  
- Max velocity: **1,191 ft/s (Mach 1.077)**  
- Max acceleration: **431 ft/s²**  
- Stability margin: **2.76 cal**  
- Simulation tool: OpenRocket with custom MATLAB post-processing  

![FAST 10k COTS rocket with aluminum fin can](assets/img/IRECRocket.PNG)

### Aluminum Fin Can – SOLIDWORKS

- **Material:** CNC/lathe-machined 6061-T6 aluminum fin can, bolted to a 4" fiberglass airframe  
- **Design target:** ≥1.5× V_max flutter margin  
- **Status:** Preliminary CAD model for the IREC 2026 fin-can design  

**Views:** Planform • Side • Isometric  

Planform view:  
![Aluminum fin can planform view](assets/img/AlFinCan3.PNG)

Side view:  
![Aluminum fin can side view](assets/img/AlFinCan2.PNG)

Isometric view:  
![Aluminum fin can isometric view](assets/img/AlFinCan.PNG)

---

## Falcon Aeronautics & Space Team

- FAST engages students at UTPB in aerospace and space-exploration projects through hands-on design, testing, and competition work.  
- The team focuses on building practical experience in rocketry, avionics, and aerospace systems while growing aerospace opportunities in the Permian Basin.  
- **FAST Instagram:** [fast.utpb](https://www.instagram.com/fast.utpb/)  

---

## Projects

### Level 1 High-Power Certification

**Role:** Designer & Builder  
**Tools:** OpenRocket, hand calculations, shop tools  
**Launch Site:** San Angelo, TX  

- Designed, built, and flew an L1 certification rocket.  
- Achieved altitude: **2,148 ft** (≈2.4% error vs. prediction), validating the simulation and modeling approach.  
- Motor: **H219T-14**  

<video width="480" controls>
  <source src="/engineering-portfolio/assets/img/IMG_4424.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>

---

### Level 2 Certification — In Progress

**Role:** Designer & Builder  
**Tools:** OpenRocket, CAD, testing hardware  

- Designing a 4" dual-deploy Level 2 rocket with an avionics bay and separate drogue/main recovery.  
- Expected altitude: **~4,200 ft** on a **J425R-14** motor.  
- Target stability margin: **~1.9 cal**.  
- Focus on robust airframe design, recovery integration, and repeatable preparation and launch procedures.

![L2 Certification Rocket](/engineering-portfolio/assets/img/L2Cert.PNG)

---

## Featured Lab Reports

A selection of lab work most relevant to thermo-fluids, energy systems, and experimental testing.

### Pelton Turbine – Performance & Efficiency

**Role:** Lab Team Member  
**Tools:** Pelton turbine rig, load cell, tachometer, Excel  

- Characterized Pelton turbine performance under varying loads and spear-valve positions.  
- Generated torque–speed curves, computed hydraulic power, and evaluated overall turbine efficiency.  
- Connected experimental results to turbomachinery performance and energy-system applications.  

**Report:** [PDF](/engineering-portfolio/assets/img/Lab2%20Pelton%20Turbine.pdf)

<p align="center">
  <img src="/engineering-portfolio/assets/img/PowerVsSpeed.png" alt="Power vs Speed for Pelton turbine" width="450">
</p>

*Figure: Measured power vs. rotational speed for the Pelton turbine at 50% and 100% spear-valve openings, showing peak power at intermediate speeds.*

---

### Performance Analysis of a Turbojet Engine (SR-30 Mini-Lab)

**Role:** Lab Team Member  
**Tools:** SR-30 Mini-Lab test stand, DAQ/LabVIEW, MATLAB, Excel  

- Processed SR-30 turbojet Mini-Lab data at a representative high-power condition using LabVIEW and Excel.  
- Computed density, velocity, and Mach number at the compressor inlet (station 1) and nozzle exit (station 5), confirming a high-speed but subsonic exhaust jet.  
- Estimated inlet and nozzle mass flow rates and checked approximate conservation of mass across the engine.  
- Built a 1-D control-volume momentum model to predict static thrust and compared it to load-cell thrust measurements (prediction ≈9.5 lbf vs. ≈11 lbf measured, about 14% low).  
- Interpreted discrepancies in terms of probe placement, non-uniform exhaust flow, and departures from ideal Brayton-cycle assumptions.  

**Report:** [PDF](/engineering-portfolio/assets/img/Performance%20Analysis%20of%20a%20Turbojet%20Engine.pdf)

<p align="center">
  <img src="/engineering-portfolio/assets/img/Gas%20Turbine.jpeg"
       alt="SR-30 turbojet engine mounted in Mini-Lab test stand"
       width="450">
</p>

*Figure: SR-30 turbojet engine installed in the Mini-Lab test stand with enclosed test cell, control gauges, and throttle lever.*

---

### Cross-Flow Heat Exchanger – TE93

**Role:** Lab Team Member  
**Tools:** TE93 cross-flow rig, pitot-static tube, DAQ, Excel  

- Measured velocity, pressure drop, and dynamic pressure across a cross-flow rod bank.  
- Used correlations to evaluate heat-transfer behavior and compared experimental and theoretical trends.  
- Demonstrated understanding of convective heat transfer, pressure loss, and airflow behavior.  

**Report:** [PDF](/engineering-portfolio/assets/img/Lab6%20Cross%20Flow%20Heat%20Exchanger%20(1).pdf)

<p align="center">
  <img src="/engineering-portfolio/assets/img/IdealVsActual.PNG"
       alt="Ideal vs actual velocity for different valve openings"
       width="450">
</p>

*Figure: Ideal vs. actual air velocity for multiple valve openings, with deviations typically within about 10% of theory.*

---

### Cooling Rate / Transient Convection

**Role:** Lab Team Member  
**Tools:** TE93 cross-flow heat exchanger, thermocouples, DAQ, Excel  

- Used the TE93 cross-flow heat exchanger to measure transient cooling of a heated copper rod.  
- Applied lumped-capacitance methods to estimate convective heat-transfer coefficients from temperature–time data.  
- Compared experimental \(h\)-values to textbook correlations and analyzed sources of error and deviation.  

**Report:** [PDF](/engineering-portfolio/assets/img/Lab7%20Cooling%20Rate.pdf)

<p align="center">
  <img src="/engineering-portfolio/assets/img/CoolingCurve_SingleRod1.PNG"
       alt="Cooling curves for heated rod at different airflow settings"
       width="450">
</p>

*Figure: Cooling curves for a single heated rod at 20–80% airflow; steeper slopes at higher airflow correspond to larger convective heat-transfer coefficients.*

---

### Fluid Friction & Bend Loss – H408

**Role:** Lab Team Member  
**Tools:** H408 apparatus, pressure gauges, Excel  

- Measured friction factors for smooth and rough pipes and multiple fittings.  
- Calculated Reynolds number, head loss, and minor loss coefficients (k, k_l) and compared them with Moody/Blasius correlations.  
- Related findings to real piping networks, pump sizing, and industrial fluid systems.  

**Report:** [PDF](/engineering-portfolio/assets/img/Lab3%20Friction%20Fluid.pdf)

<p align="center">
  <img src="/engineering-portfolio/assets/img/FrictionFactorAnalysis.png" alt="Friction factor vs Reynolds number compared to Blasius correlation" width="450">
</p>

*Figure: Experimental friction factor vs. Reynolds number for a smooth pipe, compared with the Blasius correlation, showing good agreement over the tested range.*

---

### Other Lab Reports

- **Flow Measuring Devices**  
  Compared orifice plates, Venturi meters, and rotameters; evaluated discharge coefficients and measurement uncertainties.  
  **Report:** [PDF](/engineering-portfolio/assets/img/Lab5%20Flow%20Measuring%20Devices.pdf)

- **Bend Loss in Pipe Fittings**  
  Measured head loss across elbows and other fittings; extracted minor loss coefficients (\(k\)) and compared them to standard handbook values.  
  **Report:** [PDF](/engineering-portfolio/assets/img/Lab%204%20Bend%20Loss.pdf)

- **Velocity Profile in Cross-Flow**  
  Mapped the velocity profile downstream of a rod bank using a pitot tube, resolving wake behavior and recirculation regions.  
  **Report:** [PDF](/engineering-portfolio/assets/img/Lab8-Velocity%20Profile.pdf)

- **Impact of a Jet**  
  Investigated momentum transfer and force from a water jet striking targets; compared measured forces with theoretical predictions.  
  **Report:** [PDF](/engineering-portfolio/assets/img/Lab1%20Impact%20of%20a%20Jet-1.pdf)

- **Impact: Charpy & Izod**  
  Performed impact testing on multiple 3D-printed polymers to compare absorbed energy and fracture behavior.  
  **Report:** [PDF](/engineering-portfolio/assets/img/Lab9-Impact%20Charpy%20and%20Izod.pdf)

- **Torsion Test**  
  Obtained shear modulus and torque–angle behavior for circular shafts; compared elastic and plastic response to theoretical models.  
  **Report:** [PDF](/engineering-portfolio/assets/img/lab10-%20Torsion%20Lab%20Report.pdf)

---

## Skills

### Software

- SOLIDWORKS  
- MATLAB  
- OpenRocket  
- LabVIEW / DAQ systems  
- Excel (engineering / data analysis)  

### Engineering & Domain

- Thermo-Fluids & Energy Systems  
- Heat Transfer (cross-flow, transient)  
- Piping, Pressure Losses & Head Loss  
- Mechanical / Structural Design & CAD  
- Experimental Methods, DAQ & Data Analysis  
- High-Power Rocketry, Aerodynamic Design & Stability  
- Recovery Systems (dual-deploy, parachute sizing & rigging)

### Certifications

- Tripoli / NAR Level 1 (Level 2 in progress)  

---

## Organizations

- ASME — Student Member  
- NAR — Student Member  
- TRA — Student Member  
- AIAA — Student Member  

---

## Contact

- **Email:** esparza_m58311@utpb.edu  
- **LinkedIn:** [linkedin.com/in/marcos-v-esparza](https://www.linkedin.com/in/marcos-v-esparza/)

[Back to Top](#marcos-esparza--engineering-portfolio)
