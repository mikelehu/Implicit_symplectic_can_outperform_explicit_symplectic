# SIMD-vectorized implicit symplectic integrators can outperform explicit symplectic ones

## abstract

The main purpose of this work is to present a SIMD-vectorized implementation of the symplectic 16th-order 8-stage implicit Runge-Kutta integration based on collocation with Gauss-Legendre nodes (SIMD-IRKGL16), and to show that it can outperform state-of-the-art symplectic explicit integrators for high-precision numerical integrations (in double-precision floating-point arithmetic) of non-stiff Hamiltonian ODE systems.  Our SIMD-IRKGL16 integrator leverages Single Instruction Multiple Data (SIMD) based parallelism (in a way that is transparent to the user) to significantly enhance the performance of the sequential IRKGL16 implementation.
We present numerical experiments comparing SIMD-IRKGL16 with state-of-the-art high-order explicit symplectic methods for the numerical integration of several Hamiltonian systems in double-precision floating-point arithmetic.


