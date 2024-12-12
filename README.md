# PINN-based Temperature Field Prediction for 2D Beam

## Project Overview

This project aims to apply **Physics-Informed Neural Networks (PINN)** to predict the temperature field of a 2D beam. Specifically, the model inputs are the temperature data measured by **RTDs (Resistance Temperature Detectors)** at different positions during experiments, and the output is the temperature prediction for the entire beam. The project consists of two main parts:

1. **FEM Simulation Data as PINN Input**: Before conducting actual experiments, we first solve the heat conduction model using **Finite Element Method (FEM)** in MATLAB and extract temperature data at specific locations as the input for PINN. This data can be considered "synthetic experimental data" and used to test the feasibility of the PINN model.
2. **Real Experimental Data Comparison**: Real experimental data is collected from RTDs and compared with the temperature predictions made by the PINN model to evaluate its accuracy using thermocouples' data, which are attached at the different positions of the beam.

## Project Goals

- Generate temperature data from **FEM simulations** & **Experimental data** as training data for the PINN model
- Compute and evaluate the performance of the PINN model using metrics such as R² and RMSE for selected test points along with time

## Run the code

- Beam_Sensor_Position_Plot.m visualizes the geometry of a 2D beam and the positions of sensors (RTD & TC) located on it.
  ![RTD & TC locations](https://github.com/The-Alchemist-0517/Thermal_PINN/blob/main/Example_Figures/Beam_Sensor_Position_Plot/beam_sensor.png)
- Thermal_FEM_2D.m applies FEM to a 2D transient thermal model for a rectangular beam using the MATLAB PDE Toolbox. The result is visualized using video. After generating FEM results, downsample the result for PINN training and further PINN vs. FEM comparison.
- PINN_vs_FEM_comparison.m extracts temperature predictions from the PINN model at specific x-positions and a fixed y-position and calculates statistical metrics (R² and RMSE) to quantify the agreement between the PINN predictions and FEM model.
- Thermal_exp_2D.m extracts and downsamples the temperature data from multiple RTDs, organizes them for PINN training
- PINN_vs_TC_comparison.m extracts temperature predictions from the PINN model at specific x-positions and a fixed y-position corresponding to TC sensor locations. It then calculates statistical metrics (R² and RMSE) to quantify the agreement between the PINN predictions and TC measurements over time.



