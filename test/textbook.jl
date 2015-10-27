using DiffEQ
## Some tests from an old Numerical Methods book

function test1(t, y)
    ydot = zeros(2)
    ydot[1] = 2.0*y[1] + 2.0*y[2]
    ydot[2] = 3.0*y[1] + y[2]
    return ydot
end

function actual_test1(t)
    y = zeros(2)
    y[1] = exp(4.0*t) - 2.0*exp(-t)
    y[2] = exp(4.0*t) + 3.0*exp(-t)
    return y
end

tout = linspace(0, 2, 100)
sol = ode45(test1, [-1.0, 4.0], tout)
##TODO what is the correct way to turn an array of arrays into a matrix?
actual_sol = hcat([actual_test1(t) for t in tout])
println(sol.y .- actual_sol)
