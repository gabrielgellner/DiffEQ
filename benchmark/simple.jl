using DiffEQ
using ODE

function f(t, y)
    global i
    i += 1
    ydot = similar(y)
    ydot[1] = y[2]*y[3]
    ydot[2] = -y[1]*y[3]
    ydot[3] = -0.51*y[1]*y[2]
    ydot
end

# warm up the solvers
i = 0
aode(Dopri54(f, [0.0, 1.0, 1.0]), linspace(0.0, 10, 10); reltol = 1e-4, abstol = 1e-4)
ode45(f, [0.0, 1.0, 1.0], linspace(0.0, 10, 10); reltol = 1e-4, abstol = 1e-4)

## Strangely this problem has aode faster, whereas the `vsODE.jl` version has `ODE.jl` much faster!
# matlab speed: approx 0.067503sec vs aode: 0.001025sec -- ie crazy faster
tout = linspace(0.0, 12.0, 10000)
i = 0
println("aode:")
@time sol = aode(Dopri54(f, [0.0, 1.0, 1.0]), tout; reltol = 1e-4, abstol = 1e-4)
println("aode: $i")

println("ode45:")
i = 0
@time t, y = ode45(f, [0.0, 1.0, 1.0], tout; reltol = 1e-4, abstol = 1e-4, points = :specified)
println("ode45: $i")

# two ref runs from matlab 2015b ... notice the error
#refsol = [-0.705810884230406, -0.708700436943760, 0.863899310414768]
reft = 11.993999399939995
refsol = [-0.701710883180149, -0.712461813852860, 0.865376759827586]

println("Compare:")
println("ref:   ", reft, refsol)
println("aode:  ", sol.x[end - 5], sol.y[end - 5, :])
println("ode45: ", t[end - 5], y[end - 5])
