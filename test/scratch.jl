using DiffEQ

function f(t, y)
    ydot = similar(y)
    ydot[1] = y[2]*y[3]
    ydot[2] = -y[1]*y[3]
    ydot[3] = -0.51*y[1]*y[2]
    ydot
end

t = linspace(0.0, 12.0, 10)
t, y = ode45(f, [0.0, 1.0, 1.0], t; reltol = 1e-4, abstol = 1e-4)

println(size(y))
println(y)
