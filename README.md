# Teukolsky Master Equation: Mathematical Analysis and Numerical Solution
# By:vNataly Phawllyn Neira Parra 

## Overview

This repository contains the computational framework and methodology for the study of the **Teukolsky master equation**, applied to gravitational perturbations in Kerr black holes. It focuses on understanding the perturbations in black holes, with particular attention given to Schwarzschild and Kerr geometries, their stability, and the quasinormal modes (QNM’s) associated with these solutions. The methodology used is based on **Newman-Penrose formalism**, leveraging **numerical techniques** such as **Newton iterations** and **Leaver’s method** to solve the equation.

## Objective

The primary goal of this research is to **study the mathematical structure of the Teukolsky master equation** for gravitational perturbations in Kerr black holes. This objective is pursued through a detailed examination of the equation's spectral properties and its implications for the stability of Schwarzschild and Kerr black hole solutions.

### Specific Objectives:
- **Derive the Schwarzschild and Kerr metrics** rigorously as solutions to Einstein's field equations.
- **Explore the Newman-Penrose formalism** and its application in deriving the Teukolsky equation.
- **Analyze the asymptotic behavior of solutions** to the Teukolsky equation, focusing on QNMs.
- **Investigate the mathematical and physical implications** of spectral analysis for the stability of Kerr-type solutions.

## Methods

1. **Theoretical Framework**:
    - The Schwarzschild and Kerr metrics are derived using spherical and rotational symmetries, respectively.
    - The **Newman-Penrose formalism** is employed to transform Einstein's field equations, facilitating the decoupling of variables in the Teukolsky equation.
  
2. **Numerical Techniques**:
    - **Leaver’s Method**: A continued fraction approach to solve the radial equation for perturbations.
    - **Newton Iterations**: Used to solve the angular part of the equation and calculate QNM frequencies.

3. **Computational Approach**:
    - **Python and Julia Scripts**: The Teukolsky equation is solved numerically in Julia using the **Galerkin method** and Newton’s iterations, with the **Schwarzschild limit** used for validation against known theoretical results for fields with spin \( s = 0 \) and \( s = 2 \).

## Key Results

- **Schwarzschild Solution**: The Schwarzschild solution is first validated as a special case for \( a = 0 \), where the perturbations are examined for both scalar ( \( s = 0 \) ) and gravitational ( \( s = 2 \) ) fields.
  
- **Kerr Solution**: The approach is extended to the Kerr metric, where the perturbations are more complex due to the rotation of the black hole. This necessitated the use of complex tetrads for efficient calculation.
  
- **Quasinormal Modes (QNM's)**: Numerical results show the stability of the Schwarzschild solution under linear perturbations. The QNM frequencies were computed for both Schwarzschild and Kerr geometries, and the results agree with known theoretical values for these black holes.

## Conclusion

The study of the Teukolsky master equation reveals the critical role of **quasinormal modes** in understanding the stability of black hole solutions, particularly in the Kerr geometry. The findings confirm the stability of the Schwarzschild metric under linear perturbations, while the Kerr black hole's perturbations exhibit similar behavior.

However, the analysis in this work was restricted to linear perturbations in the Kerr geometry. Future work could involve:
- Extending the analysis to **non-linear perturbations**.
- Implementing additional numerical methods to improve the **precision and stability** of the QNM calculations.
- Further exploring the **Kerr-Newman geometry** under perturbations, which presents more complex, coupled equations for perturbation analysis.

## Repository Structure

- **`Tesis.pdf`**: Provides the  complite  document of the project.
- **`notebook_calculos.ipynb`**: Contains some calculations performed in Sage.
- **`paquetes_juli.jl`**: Julia packages required for the environment setup.
- **`utilidades_y_funciones.jl`**: Julia functions used for solving the Teukolsky equation.
- **`sol_teukolsky.jl`**: Implements the numerical solution of the Teukolsky equation.



## Running the Code

This project involves solving the Teukolsky equation using multiple programming environments. The code is organized into different files for **Julia**, **Python**, and **SageMath**. Below are the instructions for running each section:

### 1. **Julia Setup**

To run the code in Julia, follow these steps:

1. **Set up the environment**:
    - Open the terminal or console and navigate to the directory where the `paquetes_juli.jl` file is located.
    - Run the following command to install necessary Julia packages:
      ```julia
      include("paquetes_juli.jl")
      ```
      This will install the required packages for solving the Teukolsky equation.

2. **Run the Utility File**:
    - Once the necessary packages are installed, run the utility functions by executing:
      ```julia
      include("utilidades_y_funciones.jl")
      ```
    - This will load all the helper functions needed for solving the Teukolsky equation.

3. **Solve the Teukolsky Equation**:
    - After loading the utilities, you can solve the Teukolsky equation by running the main solver:
      ```julia
      include("sol_teukolsky.jl")
      ```
    - This file contains the actual implementation for solving the Teukolsky equation numerically.

4. **Important Note**:
    - The `paquetes_juli.jl` file needs to be run only **once** to install the necessary packages. After the first setup, you can skip this step and go directly to running the utility and solver files.

---

### 2. **Python Setup**

To run the Python code, ensure you have the necessary Python environment set up:

1. **Install Dependencies**:
    - First, ensure that you have all the required Python libraries installed. You can do this by running:
      ```bash
      pip install -r requirements.txt
      ```
    - This will install all dependencies, including libraries for numerical computation and plotting.

2. **Run the Python Scripts**:
    - The primary Python script for solving the Teukolsky equation is located in `teukolsky_solver.py`. You can execute it by running:
      ```bash
      python teukolsky_solver.py
      ```
    - This will compute the numerical solutions for the perturbations and calculate the QNM frequencies.

---

### 3. **SageMath Setup**

To run the SageMath code, follow these steps:

1. **Install SageMath**:
    - If you don't have SageMath installed, you can download it from [SageMath's official website](http://www.sagemath.org/).

2. **Run the SageMath Notebook**:
    - The SageMath calculations are contained in the `notebook_calculos.ipynb` file.
    - You can run the notebook directly from the SageMath environment by opening it with:
      ```bash
      sage notebook notebook_calculos.ipynb
      ```
    - This will launch a SageMath notebook where you can execute the cells and view the results of the calculations.

---

### Important Notes

- For **Julia**, ensure that the first time you run the setup (`paquetes_juli.jl`), it installs the required packages. After that, you can skip it in subsequent runs.
- In **Python**, make sure to install the required libraries by running `pip install -r requirements.txt` before executing the scripts.
- **SageMath** requires a specific environment and can be run directly from its notebook interface.

By following these steps, you can set up and execute the Teukolsky equation solver in the respective environments. Each section of the project is modular, allowing you to test different methods and frameworks independently.


## References

1. Chandrasekhar, S. (1983). *The Mathematical Theory of Black Holes*. Oxford University Press.
2. Cook, G. B., & Zalutskiy, M. (2014). *Numerical Solutions of the Teukolsky Equation*. Astrophysical Journal.
3. Teukolsky, S. A. (1973). *Perturbations of a Rotating Black Hole. I. Fundamental Equations for Gravitational, Electromagnetic, and Neutrino-Field Perturbations*. Astrophysical Journal.
