# SIMD-vectorized implicit symplectic integrators can outperform explicit symplectic ones


## Overview

This repository contains the Julia implementation and numerical experiments
supporting the article: 

> *SIMD-vectorized implicit symplectic integrators can outperform explicit symplectic ones*

The goal is to demonstrate that a SIMD-vectorized implementation of the
16th-order Gauss–Legendre implicit Runge–Kutta method (IRKGL16) can outperform
state-of-the-art explicit symplectic integrators for high-precision integration
of non-stiff Hamiltonian systems.

## How to reproduce the results

This section describes how to reproduce the numerical experiments and figures
presented in the preprint.


### Environment
- Julia ≥ 1.x
- All dependencies are specified in `Project.toml` and `Manifest.toml`


### Setup
```julia
julia --project
] instantiate
```

### Running the experiments

The Jupyter notebooks in the `Experiments_Article/` directory reproduce the
figures and numerical results from the article.

- Each notebook corresponds to a specific experiment or figure.

- Notebooks can be run independently unless otherwise stated.

- Figures are generated automatically when the notebooks are executed.

## Repository Contents
 
This folder contains the code and experiments associated with the article.

- **Splitting solvers (Julia implementation)**  
  Core algorithms developed and used in the study.

- **`Experiments_Article/`**  
  Jupyter notebooks reproducing the numerical experiments and figures from the article.


- **`Examples/`**   
  Simple illustrative notebooks demonstrating how to use the solvers on test problems.

- **`Other Benchmarks/`**  
  Additional performance results for the IRKGL16 method.

## IRKGL16 Solver

This repository 
[IRKGaussLegendre.jl](https://github.com/SciML/IRKGaussLegendre.jl) contains the IRKGL16 implementation

## Citation

If you use this code or reproduce results from this repository, please cite
the associated preprint:

```bibtex
@article{Antonana2025IRKGL16,
  title   = {SIMD-vectorized implicit symplectic integrators can outperform explicit symplectic ones},
  author  = {Antonana, M., Makazaga, J., and Murua, A.},
  journal = {arXiv preprint arXiv:2511.03655},
  year    = {2025}
}
```

This reference will be updated once the article is published.



## Contact

If you have any questions or suggestions, feel free to open an issue or contact us at mikel.antonana@ehu.eus.

Updated February 23, 2026
