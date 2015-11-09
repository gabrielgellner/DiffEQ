using DiffEQ
using Base.Test

function f(t, y)
    ydot = similar(y)
    ydot[1] = y[2]*y[3]
    ydot[2] = -y[1]*y[3]
    ydot[3] = -0.51*y[1]*y[2]
    ydot
end

# the matlab version has a vector of absol = [1e-4 1e-4 1e-5]
dsol = dode(Dopri54(f, [0.0, 1.0, 1.0]), [0.0, 12.0]; reltol = 1e-4, abstol = 1e-4)

#println(dsol(12.0))
#println(dsol(linspace(0.0, 12.0, 10)))

refsol = [-0.705810884230406, -0.708700436943760, 0.863899310414768]
t2 = 9.995799579957996
refsol2 = [0.879460191254681, -0.475972448567555, 0.778164753505625]

# is this correct? if my abstol = 1e-4 do I only have accuracy to 1e-3?
@test abs(norm(refsol - dsol(12.0))) < 1e-3
@test abs(norm(refsol2 - dsol(t2))) < 1e-3
