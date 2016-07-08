# DiffEQ.jl

The `DiffEQ` package is for solving differential equations numerically in Julia. Currently
this is in very early development and really only supports using the Runge-Kutta solver
which uses the Dormand-Prince (4, 5) method. Ultimately this package aims to have the
basic set of general purpose solvers needed to solve ODE, DAE, and BVP's all coded in
pure Julia.

If you are looking for a more mature package now, excellent solvers are found in the
`ODE.jl` package (which this code is heavily based on), the `Sundials.jl` package, and
finally the `DASSL.jl` package.

To install the package run the following from the Julia prompt:
```julia
Pkg.clone("https://github.com/gabrielgellner/DiffEQ.jl")
```

## Example Usage
The solvers in `DiffEQ` are based around using core top level functions for array outputs
`aode`, and continuous extension (interpolating function) `dode`. Each of these top level
functions expects a `ODESystem` type as its first argument. Currently the only system
supported is `Dopri54`, which represents the memory and functions needed for the solver.

```julia
using DiffEQ
using PyPlot

function f(t, y)
    ydot = similar(y)
    ydot[1] = y[2]*y[3]
    ydot[2] = -y[1]*y[3]
    ydot[3] = -0.51*y[1]*y[2]
    ydot
end

tout = linspace(0.0, 12.0, 10)
# Here we define the type in the call, but you could also save the type as a variable
# and then only incur the memory allocation once for all subsequent calls to the solver.
sol = aode(Dopri54(f, [0.0, 1.0, 1.0]), tout; reltol = 1e-4, abstol = 1e-4)

plot(sol.x, sol.y)

# or if you want to have an interpolating function
csol = dode(Dopri54(f, [0.0, 1.0, 1.0]), [0.0, 1.0]; reltol = 1e-4, abstol = 1e-4))

# which can be called at any time between tstart -> tend
csol(0.5)
```
