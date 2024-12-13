# PINN-based Temperature Field Prediction for 2D Beam

## Project Overview

This project aims to apply **Physics-Informed Neural Networks (PINN)** to predict the temperature field of a 2D beam. Specifically, the model inputs are the temperature data measured by **RTDs (Resistance Temperature Detectors)** at different positions during experiments, and the output is the temperature prediction for the entire beam. The project consists of two main parts:

1. **FEM Simulation Data as PINN Input**: Before conducting actual experiments, we first solve the heat conduction model using **Finite Element Method (FEM)** in MATLAB and extract temperature data at specific locations as the input for PINN. This data can be considered "synthetic experimental data" and used to test the feasibility of the PINN model.
2. **Real Experimental Data Comparison**: Real experimental data is collected from RTDs and compared with the temperature predictions made by the PINN model to evaluate its accuracy using thermocouples' data, which are attached at the different positions of the beam.

## Project Goals

- Generate temperature data from **FEM simulations** & **Experimental data** as training data for the PINN model
- Compute and evaluate the performance of the PINN model using metrics such as R² and RMSE for selected test points along with time
<table align="center">
  <tr>
    <td align="center">
      <img src="https://github.com/The-Alchemist-0517/Thermal_PINN/blob/main/Example_Figures/MS4_1.png" 
           alt="RTD 2(2 layers)" width="300">
      <p>Experiment set up</p>
    </td>
    <td align="center">
      <img src="https://github.com/The-Alchemist-0517/Thermal_PINN/blob/main/Example_Figures/MS4_2.png" 
           alt="RTD 3(2 layers)" width="300">
      <p>PINN working process</p>
    </td>
    <td align="center">
      <img src="https://github.com/The-Alchemist-0517/Thermal_PINN/blob/main/Example_Figures/MS4_3.png" 
           alt="RTD 4(2 layers)" width="300">
      <p>Prediction results</p>
    </td>
  </tr>
</table>

## Run the code
- **Beam_Sensor_Position_Plot.m**: This script visualizes the geometry of a 2D beam and the positions of sensors (RTD & TC) located on it.
<table align="center">
  <tr>
    <td align="center">
      <img src="https://github.com/The-Alchemist-0517/Thermal_PINN/blob/main/Example_Figures/Beam_Sensor_Position_Plot/beam_sensor.png" 
           alt="FEM Temperature Field Video" width="600">
      <p>RTD & TC locations </p>
    </td>
  </tr>
</table>

- **Thermal_FEM_2D.m**: This script applies FEM to a 2D transient thermal model for a rectangular beam using the MATLAB PDE Toolbox. The result is visualized as a video. After generating FEM results, downsample the result for PINN training and further comparison.

<table align="center">
  <tr>
    <td align="center">
      <img src="https://github.com/The-Alchemist-0517/Thermal_PINN/blob/main/Example_Figures/Thermal_FEM_2D/FEM_temperature_field-ezgif.com-video-to-gif-converter.gif" 
           alt="FEM Temperature Field Video" width="600">
      <p>FEM Temperature Field </p>
    </td>
  </tr>
</table>

- PINN_vs_FEM_comparison.m extracts temperature predictions from the PINN model at specific x-positions and a fixed y-position and calculates statistical metrics (R² and RMSE) to quantify the agreement between the PINN predictions and FEM model.
  
<table align="center">
  <tr>
    <td align="center">
      <img src="https://github.com/The-Alchemist-0517/Thermal_PINN/blob/main/Example_Figures/PINN_vs_FEM_comparison/PINN_baseline_200s_3_points_Transient-ezgif.com-video-to-gif-converter.gif" 
           alt="PINN Baseline Transient" width="400">
      <p>PINN Temperature Filed Prediction </p>
    </td>
    <td align="center">
      <img src="https://github.com/The-Alchemist-0517/Thermal_PINN/blob/main/Example_Figures/PINN_vs_FEM_comparison/Difference_baseline_0_200s_3_points_Trasient_input-ezgif.com-video-to-gif-converter.gif" 
           alt="Difference Between PINN and FEM" width="400">
      <p>Difference Baseline Transient</p>
    </td>
  </tr>
</table>

- Thermal_exp_2D.m extracts and downsamples the temperature data from multiple RTDs and organizes them for PINN training.2 layers means we have 2 layers of embedded RTDs in the beam, but for this 2D case, we only use 1 layer of RTDs' data for training. 

<table align="center">
  <tr>
    <td align="center">
      <img src="https://github.com/The-Alchemist-0517/Thermal_PINN/blob/main/Example_Figures/Thermal_exp_2D/RTD_2.png" 
           alt="RTD 2(2 layers)" width="300">
      <p>RTD 2(2 layers)</p>
    </td>
    <td align="center">
      <img src="https://github.com/The-Alchemist-0517/Thermal_PINN/blob/main/Example_Figures/Thermal_exp_2D/RTD_3.png" 
           alt="RTD 3(2 layers)" width="300">
      <p>RTD 3(2 layers)</p>
    </td>
    <td align="center">
      <img src="https://github.com/The-Alchemist-0517/Thermal_PINN/blob/main/Example_Figures/Thermal_exp_2D/RTD_4.png" 
           alt="RTD 4(2 layers)" width="300">
      <p>RTD 4(2 layers)</p>
    </td>
  </tr>
</table>

- PINN_vs_TC_comparison.m extracts temperature predictions from the PINN model at specific x-positions and a fixed y-position corresponding to TC sensor locations. It then calculates statistical metrics (R² and RMSE) to quantify the agreement between the PINN predictions and TC measurements over time.

<table align="center">
  <tr>
    <td align="center">
      <img src="https://github.com/The-Alchemist-0517/Thermal_PINN/blob/main/Example_Figures/PINN_vs_TC_Comparison/TC_error.png" 
           alt="FEM Temperature Field Video" width="600">
      <p>PINN temperature prediction vs. TC measurement  </p>
    </td>
  </tr>
</table>



