using DiffEQ
using ODE

function f(t, y)
    ydot = similar(y)
    ydot[1] = y[2]*y[3]
    ydot[2] = -y[1]*y[3]
    ydot[3] = -0.51*y[1]*y[2]
    ydot
end

# warm up the solvers
aode(Dopri5(f, [0.0, 1.0, 1.0]), linspace(0.0, 10, 10); reltol = 1e-4, abstol = 1e-4)
ode45(f, [0.0, 1.0, 1.0], linspace(0.0, 10, 10); reltol = 1e-4, abstol = 1e-4)

## Strangely this problem has aode faster, whereas the `vsODE.jl` version has `ODE.jl` much faster!
# matlab speed: approx 0.067503sec vs aode: 0.001025sec -- ie crazy faster
tout = linspace(0.0, 12.0, 10000)
println("aode:")
@time sol = aode(Dopri5(f, [0.0, 1.0, 1.0]), tout; reltol = 1e-4, abstol = 1e-4)
println("ode45:")
@time t, y = ode45(f, [0.0, 1.0, 1.0], tout; reltol = 1e-4, abstol = 1e-4, points = :specified)

# two ref runs from matlab 2015b ... notice the error
#refsol = [-0.705810884230406, -0.708700436943760, 0.863899310414768]
refsol = [-0.705653253786900, -0.708731728427065, 0.863900202621616]

println("solutions")
println(sol.y[end - 1, :])
println(y[end - 1])
println("ref")
println(refsol)
println(sol.y[end, :])
println(y[end])
