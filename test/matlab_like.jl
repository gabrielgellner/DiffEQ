using DiffEQ
using Base.Test

function f(t, y)
    ydot = similar(y)
    ydot[1] = y[2]*y[3]
    ydot[2] = -y[1]*y[3]
    ydot[3] = -0.51*y[1]*y[2]
    ydot
end

t = linspace(0.0, 12.0, 10)
sol = ode45(f, [0.0, 1.0, 1.0], t; reltol = 1e-4, abstol = 1e-4)

@test size(sol.y) == (10, 3)
#@test sol.y == 1
