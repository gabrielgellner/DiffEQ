module DiffEQ

# So to start I will just port over the ode45 from `ODE.jl` so that I can work
# on the interface that I want to use.

export ode45

#TODO: So I am doing all the package import logic in this file, which feels
# super ugly me, need to see if there is a better way
include("types.jl")
include("runge_kutta_tableaus.jl")
include("utils.jl")
include("runge_kutta_minimal.jl")

end # module