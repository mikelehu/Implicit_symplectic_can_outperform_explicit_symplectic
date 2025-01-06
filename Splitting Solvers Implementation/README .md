# Splitting methods



## Description

- Only constant step size

## Installation


## Solver options


### Available common arguments

- dt:  is the step size for the integration
- save_everystep: specifies whether intermediate solutions are saved (default is true)

### No-common arguments

- mstep: output saved at every 'mstep' steps (default mstep=1)


### Splitting Algorithm 

 - r : order of the method
 - rkn : =true for Runge-Kutta-Nystrom methods; false= otherwise


### Splitting Problem

  - flows: array of functions
  - u0:
  - tspan:
  - p: