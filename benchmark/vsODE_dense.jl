using DiffEQ
using ODE

function f(t, y)
    ydot = similar(y)
    ydot[1] = -0.04*y[1] + 1.0e4*y[2]*y[3]
    ydot[3] = 3.0e7*y[2]*y[2]
    ydot[2] = -ydot[1] - ydot[3]
    ydot
end

t = linspace(0.0, 1e3, 1000)

# warm up both solvers
dode(Dopri54(f, [1.0, 0.0, 0.0]), [0.0, 1.0])
ode45(f, [1.0, 0.0, 0.0], [0.0, 1.0])

#TODO: so currently dode is much slower than ode45. I need to figure out why.
@time dsol = dode(Dopri54(f, [1.0, 0.0, 0.0]), [0.0, 1e3]; abstol=1e-8, reltol=1e-8, maxstep=1e11/10, minstep=1e11/1e18);
@time t2, sol2 = ode45(f, [1.0, 0.0, 0.0], [0.0, 1e3]; abstol=1e-8, reltol=1e-8, maxstep=1e11/10, minstep=1e11/1e18);

#TODO: I also need to think about testing the speed of the interpolater
println(dsol(1e3))
println(sol2[end])
