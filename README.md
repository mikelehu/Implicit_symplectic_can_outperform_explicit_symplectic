# SIMD-vectorized implicit symplectic integrators can outperform explicit symplectic ones

## abstract

The main purpose of this work is to present a SIMD-vectorized implementation of the symplectic 16th-order 8-stage implicit Runge-Kutta integrator based on collocation with Gauss-Legendre nodes (IRKGL16), and to show that it can outperform state-of-the-art symplectic explicit integrators for high-precision numerical integrations (in double-precision floating-point arithmetic) of non-stiff Hamiltonian ODE systems.  Our IRKGL16 integrator, implemented in Julia language, leverages Single Instruction Multiple Data (SIMD) based parallelism (in a way that is transparent to the user) to significantly enhance the performance of the sequential IRKGL16 implementation.
We present numerical experiments comparing IRKGL16 with state-of-the-art high-order explicit symplectic methods for the numerical integration of several Hamiltonian systems in double-precision floating-point arithmetic.


## Repository Contents

### 1) Splitting Solvers Implementation  
This folder contains the code and experiments associated with the article.

- **Julia implementation of splitting solvers**  
  Core algorithms developed and used in the study.

- **Jupyter notebooks reproducing the experiments from the article**  
  Located in the **`Experiments_Article`** folder. 
  Running these notebooks generates the figures and numerical results presented in the article.

- **Examples**  
  The **`Examples`** folder provides illustrative notebooks demonstrating how to use the solvers on simple test problems.

### 2) Other Benchmarks  
- Additional benchmark results for the **IRKGL16** method.

## Contact

If you have any questions or suggestions, feel free to open an issue or contact us at mikel.antonana@ehu.eus.

Updated November 28, 2025
