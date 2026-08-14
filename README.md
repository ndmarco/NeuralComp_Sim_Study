# NeuralComp_Sim_Study
Simulation Studies for "Modeling Neural Switching via Drift-Diffusion Models"


**Required Packages**: NeuralComp, MASS, truncnorm, ggplot2, latex2exp, future.apply, ggpubr, scales, transport, ggallin, statmod, gridExtra, grid

**How to Run**: 
1. Run WAIC Simulation Study 1 (Section 3.2 of the Supplementary Materials) using `WAIC_Simulations/Sim1.R`.
    - Load packages by running lines 1-13.
    - Set save directory (`save_dir`) in line 15 and run.
    - Run competition model part of simulation study by running lines 17-172.
    - Run lines 174-273 to obtain Figure 5 of the Supplementary Materials.
    - Run IIGPP model part of simulation study by running lines 279-402.
    - Run lines 407-476 to obtain Figure 6 of the Supplementary Materials.
2. Run WAIC Simulation Study 2 (Section 3.3 of the Supplementary Materials) using `WAIC_Simulations/Sim2.R`.
    - Load packages by running lines 1-14.
    - Set save directory (`save_dir`) in line 16 and run.
    - Run simulation study by running lines 20-249.
    - Obtain figures by running lines 252-277.
3. Run WAIC Simulation Study 3 (Section 3.4 of the Supplementary Materials) using `WAIC_Simulations/Sim3.R`.
    - Load packages by running lines 1-16.
    - Set save directory (`save_dir`) in line 18 and run.
    - Run simulation study by running lines 20-170.
    - Obtain figures by running lines 172-235.
4. Run WAIC Simulation Study 4 (Section 3.5 of the Supplementary Materials) using `WAIC_Simulations/Sim4.R`.
    - Load packages by running lines 1-14.
    - Set save directory (`save_dir`) in line 16 and run.
    - Run winner take all part of simulation study by running lines 20-187.
    - Obtain Table 1 of the Supplementary Materials by running lines 191-212.
    - Run the competition model part of the simulation study by running lines 220-382.
    - Obtain Table 2 of the Supplementary Materials by running lines 385-404.
5. Run Recovery of Scientific Quantities simulation study (Section 4.1 of the Supplementary Materials) using `Sim_Recovery.R`.
    - Load packages by running lines 1-14.
    - Set save directory (`save_dir`) in line 16 and run.
    - Run simulation study by running lines 18-226.
    - Obtain figures by running lines 230-385.
    - Obtain Table 3 by running lines 387-414.
6. Run Prior Sensitivity Analysis (Section 4.2 of the Supplementary Materials) using `Prior_Sensitivity/Prior_Sensitivity.r`.
    - Load packages by running lines 1-14.
    - Set save directory (`save_dir`) in line 16 and run.
    - Run simulation study by running lines 18-233.
    - Obtain figures by running lines 237-386.
7. Run Simulation Study on Posterior Predictive P-Values (Section 5.2 of the Supplementary Materials) using `Posterior_Pval_Simulations/Sim.R`.
    - Load packages by running lines 1-15.
    - Set save directory (`save_dir`) in line 17 and run.
    - Run simulation study by running lines 21-166.
    - Obtain figures by running lines 169-251.
