# PINN-based Temperature Field Prediction for 2D Beam

## Project Overview

This project aims to apply **Physics-Informed Neural Networks (PINN)** to predict the temperature field of a 2D beam. Specifically, the model inputs are the temperature data measured by **RTDs (Resistance Temperature Detectors)** at different positions during experiments, and the output is the temperature prediction for the entire beam. The project consists of two main parts:

1. **FEM Simulation Data as PINN Input**: Before conducting actual experiments, we first solve the heat conduction model using **Finite Element Method (FEM)** in MATLAB and extract temperature data at specific locations as the input for PINN. This data can be considered "synthetic experimental data" and used to test the feasibility of the PINN model.
2. **Real Experimental Data Comparison**: Real experimental data is collected from RTDs and compared with the temperature predictions made by the PINN model to evaluate its accuracy.

## Project Goals

- Generate temperature data from **FEM simulations** & **Experimental data** as training data for the PINN model
- Compute and evaluate the performance of the PINN model using metrics such as R² and RMSE


