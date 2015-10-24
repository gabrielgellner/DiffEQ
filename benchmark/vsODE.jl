import DiffEQ
import ODE

function f(t, y)
    ydot = similar(y)
    ydot[1] = -0.04*y[1] + 1.0e4*y[2]*y[3]
    ydot[3] = 3.0e7*y[2]*y[2]
    ydot[2] = -ydot[1] - ydot[3]
    ydot
end

t = linspace(0.0, 1e3, 1000)

# warm up both solvers
DiffEQ.ode45(f, [1.0, 0.0, 0.0], linspace(0, 1, 10))
ODE.ode45(f, [1.0, 0.0, 0.0], linspace(0, 1, 10))

# to the real benchmarks
@time sol = DiffEQ.ode45(f, [1.0, 0.0, 0.0], t; abstol=1e-8, reltol=1e-8, maxstep=1e11/10, minstep=1e11/1e18);
@time t2, sol2 = ODE.ode45(f, [1.0, 0.0, 0.0], t; abstol=1e-8, reltol=1e-8, maxstep=1e11/10, minstep=1e11/1e18, points = :specified);
