using ODE
using DiffEQ

function rigid(t, y)
    ydot = zeros(3)
    ydot[1] = y[2]*y[3]
    ydot[2] = -y[1]*y[3]
    ydot[3] = -0.51*y[1]*y[2]
    ydot
end

tout = linspace(0.0, 12.0, 100)
y0 = [0.0, 1.0, 1.0]

# the rigid system from matlab using atol = rtol = 1e-13
reft = 12.0
refsol = [-0.705397809522307, -0.708811632467256, 0.863846690370256]

tol = 1e-5
aode(Dopri54(rigid, y0), tout, abstol = tol, reltol = tol)
