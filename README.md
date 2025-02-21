# Future forest disturbance scenarios for Europe: Code Repository

This repository contains the code and resources to create and analyze future disturbance scenarios for Europe.



The code is organised in several key steps:

## Project Workflow

### 1. Creating Training Data for a Deep Neural Network (DNN)
Prepare datasets to train the DNN for predicting forest state transitions. This step relies heavily on the **European Forest Simulation Database**, available [here](https://zenodo.org/records/10730807).

### 2. Training the Deep Neural Network
Develop and fine-tune a neural network to model forest state transitions based on training data. Details are explained in Grünig et al. (in review), with code published in this repository: [europe_tree_competition_publi](https://github.com/magrueni/europe_tree_competition_publi).

### 3. Creating Disturbance Modules
Generate disturbance modules for **wind**, **wildfires**, **bark beetles**, and **forest management**. These modules are based on process understanding and statistical models calibrated on historical disturbance observations published by Senf & Seidl (2021).

### 4. Creating the Initial Forest Landscape
Set up the baseline forest landscape for simulations across Europe for the year 2020. This involves creating a landscape where each pixel represents a forest state recognized by the DNN, allowing the model to predict state transitions. Published datasets were used as inputs for this step.

### 5. Using the Scaling Vegetation Dynamics (SVD) Modelling Framework
Implement the **SVD framework** to simulate vegetation dynamics under various scenarios. The framework was introduced by Rammer & Seidl (2019) and is available [here](https://svdmodel.github.io/SVD/#/?id=svd-documentation).

### 6. Processing SVD Outputs
Process, analyze, and refine the outputs generated by the SVD simulations to enable further insights and visualizations.

### 7. Analyzing Disturbance Scenarios
Evaluate the disturbance scenarios to assess their impacts, with a focus on the dynamics of **young** and **old forests**.

## How to Use

1. Clone the repository:
```bash
git clone https://github.com/magrueni/disturbance_scenarios_europe.git
```

## Dependencies

- **Python**: Version 3.9 or higher  
- **R**: Version 4.4.1 or higher  


## Acknowledgments

This project utilizes the **Scaling Vegetation Dynamics (SVD)** framework for simulating vegetation dynamics under future scenarios (see https://github.com/edfm-tum/SVD). 
