# Astronomy example from Hairer et al 1992 p.129 (see also Szebehely
# 1967): the orbit of a satellite in planet-moon system.  This system
# has know periodic solutions, so called “Arenstorf orbits”
# (Arenstorf, 1963); one of which is used as a test-case.
#
# This is a non-stiff problem.  Adaptive solvers will fare much better
# as there are three regions where the satellite experiences very
# strong acceleration interspersed with regions of very small
# acceleration
## Adapted from `IVPTestSuite.jl`
using DiffEQ
using Base.Test

function fn(t, y, μ, μprime)
    dydt = zeros(4)
    D1 = ((y[1] + μ)^2 + y[2]^2)^(3/2)
    D2 = ((y[1] - μprime)^2 + y[2]^2)^(3/2)
    dydt[1] = y[3]
    dydt[2] = y[4]
    dydt[3] = y[1] + 2*y[4] - μprime*(y[1] + μ)/D1 - μ*(y[1] - μprime)/D2
    dydt[4] = y[2] - 2*y[3] - μprime*y[2]/D1 - μ*y[2]/D2
    return dydt
end

y0 = [0.994, 0.0, 0.0, -2.00158510637908252240537862224]
tspan = [0, 17.0652165601579625588917206249]
refsol = y0 # reference solution at tspan[2]

μ = 0.012277471
μprime = 1 - μ
tout = linspace(tspan[1], tspan[2], 100)
sol = aode(Dopri54((t, y) -> fn(t, y, μ, μprime), y0), tout; reltol = 1e-8)

#println(sol.y[end, :])
#println(refsol)
#println(abs(refsol - sol.y[end, :]'))

##TODO: I need to set the reltol = 1e-8 to get to 1e-4 abstol (abstol = 1e-8 by default)
## I need to figure out what error is expected for the input tolerances
@test all(abs(refsol - sol.y[end, :]') .< 1e-4)
