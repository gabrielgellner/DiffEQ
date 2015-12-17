using DiffEQ
using ODE

function f(t, y)
    global i
    i += 1
    ydot = zeros(3)
    ydot[1] = -0.04*y[1] + 1.0e4*y[2]*y[3]
    ydot[3] = 3.0e7*y[2]*y[2]
    ydot[2] = -ydot[1] - ydot[3]
    ydot
end

i = 0
# warm up both solvers
aode(Dopri54(f, [1.0, 0.0, 0.0]), linspace(0, 1, 10))
ode45(f, [1.0, 0.0, 0.0], linspace(0, 1, 10))

# to the real benchmarks
## Amazingly this takes around 5sec vs 72sec on matlab 2015a!
## Strangely the ode45 code is faster, yet for other problems aode is faster ... try to figure out how
tout = linspace(0.0, 1e3, 1000)
i = 0
@time sol = aode(Dopri54(f, [1.0, 0.0, 0.0]), tout; abstol=1e-8, reltol=1e-8, maxstep=1e11/10, minstep=1e11/1e18);
println("aode: $i")

i = 0
@time t2, sol2 = ode45(f, [1.0, 0.0, 0.0], tout; abstol=1e-8, reltol=1e-8, maxstep=1e11/10, minstep=1e11/1e18, points = :specified);
println("ode45: $i")

println("aode:  ", sol.y[end, :])
println("ode45: ", sol2[end])
