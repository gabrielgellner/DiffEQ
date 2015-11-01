using DiffEQ
using Gadfly

function linear(t, y, A::Array{Float64, 2})
    yprime = A*y
    return yprime
end

function solution_linear(t, y0, A)
    return exp(t*A)*y0
end

A1 = [2.0 2.0; 3.0 0.0]
y01 = [-1.0, 4.0]

tout = linspace(0, 2, 100)
approx = aode(Dopri5((t, y) -> linear(t, y, A1), y01), tout)
analytic = hcat([solution_linear(t, y01, A1) for t in tout]...)'

#println(analytic - approx.y)
println(analytic[end, :])
println(approx.y[end, :])

p_anal1 = plot(x = tout, y = analytic[:, 1])
draw(PNG("anal1.png", 12cm, 6cm), p_anal1)
p_approx1 = plot(x = approx.x, y = approx.y[:, 1])
draw(PNG("approx1.png", 12cm, 6cm), p_approx1)
