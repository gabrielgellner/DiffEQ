using DiffEQ
using Base.Test

function linear(t, y, A::Array{Float64, 2})
    yprime = A*y
    return yprime
end

function solution_linear(t, y0, A)
    return expm(A*t)*y0
end

A1 = [2.0 2.0; 3.0 1.0]
y01 = [-1.0, 4.0]

tout = linspace(0, 2, 100)
approx = aode(Dopri54((t, y) -> linear(t, y, A1), y01), tout)
analytic = hcat([solution_linear(t, y01, A1) for t in tout]...)'

# So this is about the same as Matlab. I need to figure out how to actually
# check for this working, since it can diverge a lot as in this case we are
# looking at exponentially growing solutions. It might make sense
# to keep this contrained to stable solutions
println(sum(abs(analytic[:, 1] - approx.y[:, 1]) .< 1e-3))

#println(analytic - approx.y)
println(analytic[end, :])
println(approx.y[end, :])

#@test abs(norm(analytic[:, 1] - approx.y[:, 1])) < 1e-3
